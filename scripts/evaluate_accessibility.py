"""Does target-site accessibility predict knockdown?

The proposal names "accessible binding regions" as a parameter worth
learning. This tests it on real data before anything larger is built.

Design note that matters more than the result. Only a minority of ASO Atlas
oligos locate in the canonical cDNA: for DGAT2, 96 of 400. Orientation is
not the problem, since the reverse complement matches and the forward
sequence does not. The reason is that gapmers frequently target pre-mRNA,
and a cDNA contains no introns. Computing accessibility for the rest would
need genomic sequence, which is a far larger fetch and a far slower fold.

So this is deliberately the cheap version: evaluate accessibility on the
sites that do map, and only invest in pre-mRNA if it earns it. If
accessibility adds nothing where it is computable, it will not start adding
something on a bigger download.

Gated like everything else: target-gene-grouped folds, compared against the
same features without accessibility, so the question is whether it adds
anything rather than whether it correlates.

    python scripts/evaluate_accessibility.py

Writes reports/results/accessibility.csv and .json.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from oligotox.efficacy.potency import guide_features  # noqa: E402
from oligotox.features.accessibility import site_accessibility  # noqa: E402

ATLAS = ROOT / "data" / "external" / "aso_atlas.parquet"
TRANSCRIPTS = ROOT / "data" / "external" / "transcripts"
RESULTS = ROOT / "reports" / "results"
MIN_MAPPED_PER_GENE = 30
CLIP = (0.0, 100.0)


def _model():
    return make_pipeline(StandardScaler(), Ridge(alpha=1.0))


def load_transcripts() -> dict[str, str]:
    out = {}
    for path in sorted(TRANSCRIPTS.glob("*.json")):
        if path.name.startswith("_"):
            continue
        try:
            record = json.loads(path.read_text())
        except json.JSONDecodeError:
            continue
        sequence = (record.get("sequence") or "").upper().replace("U", "T")
        if sequence:
            out[record.get("symbol", path.stem)] = sequence
    return out


def build() -> tuple[pd.DataFrame, dict]:
    transcripts = load_transcripts()
    if not transcripts:
        raise SystemExit("no transcripts cached; run scripts/fetch_transcripts.py")

    raw = pd.read_parquet(ATLAS).dropna(
        subset=["aso_sequence_5_to_3", "inhibition_percent", "target_gene"]
    )
    raw["sequence"] = raw["aso_sequence_5_to_3"].astype(str).str.upper().str.strip()
    raw = raw[raw["sequence"].str.fullmatch(r"[ACGT]{12,30}")]
    raw = raw[raw["target_gene"].isin(transcripts)]
    raw["y"] = raw["inhibition_percent"].clip(*CLIP)
    agg = raw.groupby(["sequence", "target_gene"], as_index=False).agg(y=("y", "mean"))

    rows, mapped = [], 0
    for sequence, gene, y in agg[["sequence", "target_gene", "y"]].itertuples(index=False):
        access = site_accessibility(sequence, transcripts[gene])
        if not access.site_found:
            continue
        mapped += 1
        rows.append({
            "sequence": sequence, "target_gene": gene, "y": y,
            **guide_features(sequence), **access.as_features(),
        })

    coverage = {
        "genes_with_transcript": len(transcripts),
        "candidate_sequences": int(len(agg)),
        "sites_located": mapped,
        "coverage": round(mapped / max(len(agg), 1), 4),
        "note": "Unlocated sites are mostly intronic: a cDNA has no introns, and "
                "gapmers frequently target pre-mRNA.",
    }
    return pd.DataFrame(rows), coverage


def main() -> None:
    frame, coverage = build()
    print(f"{coverage['sites_located']:,} of {coverage['candidate_sequences']:,} "
          f"sequences located in a transcript ({coverage['coverage']:.1%})")

    counts = frame["target_gene"].value_counts()
    frame = frame[frame["target_gene"].isin(counts[counts >= MIN_MAPPED_PER_GENE].index)]
    if frame.empty or frame["target_gene"].nunique() < 3:
        raise SystemExit("not enough mapped sites per gene to evaluate; fetch more transcripts")

    print(f"evaluating {len(frame):,} sequences over "
          f"{frame['target_gene'].nunique()} genes\n")

    access_cols = ["site_accessibility", "site_log_accessibility",
                   "site_flank_accessibility", "site_relative_position"]
    design_cols = [c for c in frame.columns
                   if c not in {"sequence", "target_gene", "y", "site_found", *access_cols}
                   and pd.api.types.is_numeric_dtype(frame[c])]

    y = frame["y"].values
    genes = frame["target_gene"].values
    cv = GroupKFold(n_splits=min(5, frame["target_gene"].nunique()))

    rows = []
    for label, cols in {
        "Design features only": design_cols,
        "Design + accessibility": design_cols + access_cols,
        "Accessibility only": access_cols,
    }.items():
        pred = cross_val_predict(_model(), frame[cols].astype(float), y, cv=cv, groups=genes)
        rho = float(spearmanr(y, pred).statistic)
        rows.append({"model": label, "n_features": len(cols), "spearman": round(rho, 4),
                     "n": int(len(y)), "n_genes": int(frame["target_gene"].nunique())})
        print(f"  {label:26s} {len(cols):3d} feat  rho={rho:+.3f}")

    design = next(r for r in rows if r["model"] == "Design features only")["spearman"]
    combined = next(r for r in rows if r["model"] == "Design + accessibility")["spearman"]
    gain = round(combined - design, 4)
    verdict = ("accessibility adds signal" if gain > 0.01 else
               "accessibility adds nothing here")
    print(f"\n  gain from accessibility: {gain:+.4f} -> {verdict}")

    RESULTS.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(RESULTS / "accessibility.csv", index=False)
    (RESULTS / "accessibility.json").write_text(json.dumps({
        "coverage": coverage,
        "gain_from_accessibility": gain,
        "verdict": verdict,
        "method": "ViennaRNA pfl_fold_up (RNAplfold), window 80, max span 40. "
                  "Folds grouped by target gene.",
        "results": rows,
    }, indent=2) + "\n")
    print(f"\nwrote {(RESULTS / 'accessibility.csv').relative_to(ROOT)}")


if __name__ == "__main__":
    main()
