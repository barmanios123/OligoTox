"""Per-context potency models for ASO Atlas.

Pooled, ASO Atlas scores worse than the 3,867-sequence model it was meant to
replace (0.18 against 0.26) because percent inhibition is not comparable
across 81 cell lines and many dosages. Within one experimental context it
recovers: A431 at a single dose reaches 0.29 on 8% of the data.

So this trains one model per context rather than one model over everything,
and gates each the same way every other model here is gated. A context ships
only if it beats, within its own data:

    the target-gene-mean baseline   it must read the oligo, not the gene
    the GC-only baseline            it must beat one crude number
    a floor of 0.15 Spearman        below that it is not worth presenting

Folds are grouped by target gene inside each context. A context is only
considered at all with at least MIN_ROWS rows and MIN_GENES genes, because
grouped CV needs several genes to hold out.

What this does not claim: these are not a general potency predictor. Each
one predicts knockdown in its own assay, and the config records which. The
app labels them that way.

    python scripts/train_context_potency.py

Writes artifacts/models/potency_<context>.joblib,
artifacts/config/potency_contexts.json and reports/results/context_potency.csv.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import joblib
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

ATLAS = ROOT / "data" / "external" / "aso_atlas.parquet"
MODELS = ROOT / "artifacts" / "models"
CONFIG = ROOT / "artifacts" / "config"
RESULTS = ROOT / "reports" / "results"

MIN_ROWS = 2000
MIN_GENES = 5
FLOOR = 0.15
CLIP = (0.0, 100.0)
N_SPLITS = 5


def _model():
    return make_pipeline(StandardScaler(), Ridge(alpha=1.0))


def _normalise_cell_line(name: str) -> str:
    """A431 and A-431 are the same line written two ways."""
    return re.sub(r"[^A-Z0-9]", "", str(name).upper())


def _slug(cell_line: str, dose: float, method: str) -> str:
    dose_part = "anydose" if dose < 0 else f"{dose:g}"
    method_part = re.sub(r"[^a-z0-9]", "", str(method).lower())[:12] or "unknown"
    return f"{_normalise_cell_line(cell_line).lower()}_{dose_part}_{method_part}"


def load_contexts() -> list[tuple[str, pd.DataFrame, dict]]:
    raw = pd.read_parquet(ATLAS).dropna(
        subset=["aso_sequence_5_to_3", "inhibition_percent", "target_gene", "cell_line"]
    )
    raw["sequence"] = raw["aso_sequence_5_to_3"].astype(str).str.upper().str.strip()
    raw = raw[raw["sequence"].str.fullmatch(r"[ACGT]{12,30}")]
    raw["cell"] = raw["cell_line"].map(_normalise_cell_line)
    raw["dose"] = raw["dosage"].fillna(-1.0)
    raw["method"] = raw["transfection_method"].fillna("unknown").astype(str)
    raw["y"] = raw["inhibition_percent"].clip(*CLIP)

    contexts = []
    for (cell, dose, method), block in raw.groupby(["cell", "dose", "method"]):
        if len(block) < MIN_ROWS or block["target_gene"].nunique() < MIN_GENES:
            continue
        # One row per sequence-gene pair, so replicates never span a fold.
        agg = block.groupby(["sequence", "target_gene"], as_index=False).agg(y=("y", "mean"))
        if len(agg) < MIN_ROWS or agg["target_gene"].nunique() < MIN_GENES:
            continue
        meta = {"cell_line": cell, "dose": float(dose), "transfection_method": method}
        contexts.append((_slug(cell, dose, method), agg, meta))
    return sorted(contexts, key=lambda c: -len(c[1]))


def evaluate(frame: pd.DataFrame) -> dict:
    feats = pd.DataFrame([guide_features(s) for s in frame["sequence"]], index=frame.index)
    cols = list(feats.columns)
    y = frame["y"].values
    genes = frame["target_gene"].values
    cv = GroupKFold(n_splits=min(N_SPLITS, frame["target_gene"].nunique()))

    def rho(X):
        pred = cross_val_predict(_model(), X.astype(float), y, cv=cv, groups=genes)
        value = spearmanr(y, pred).statistic
        return float(value) if np.isfinite(value) else 0.0

    # Gene-mean baseline: no sequence information at all.
    gene_pred = np.full(len(frame), np.nan)
    for train_idx, test_idx in cv.split(frame, y, genes):
        gene_pred[test_idx] = y[train_idx].mean()
    gene_rho = spearmanr(y, gene_pred).statistic

    return {
        "model": rho(feats[cols]),
        "gc_only": rho(feats[["gc_content"]]),
        "gene_mean": float(gene_rho) if np.isfinite(gene_rho) else 0.0,
        "features": cols,
        "frame_features": feats,
    }


def main() -> None:
    if not ATLAS.exists():
        raise SystemExit(f"missing {ATLAS}; run scripts/ingest_aso_atlas.py first")

    contexts = load_contexts()
    print(f"{len(contexts)} candidate context(s) with >= {MIN_ROWS} rows "
          f"and >= {MIN_GENES} genes\n")

    rows, shipped = [], {}
    MODELS.mkdir(parents=True, exist_ok=True)

    for slug, frame, meta in contexts:
        scores = evaluate(frame)
        beats_gene = scores["model"] > scores["gene_mean"]
        beats_gc = scores["model"] > scores["gc_only"]
        above_floor = scores["model"] >= FLOOR
        passes = beats_gene and beats_gc and above_floor

        reasons = []
        if not beats_gene:
            reasons.append(f"gene-mean baseline scores {scores['gene_mean']:.3f}")
        if not beats_gc:
            reasons.append(f"GC alone scores {scores['gc_only']:.3f}")
        if not above_floor:
            reasons.append(f"below the {FLOOR} floor")

        rows.append({
            "context": slug, **meta, "sequences": int(len(frame)),
            "genes": int(frame["target_gene"].nunique()),
            "spearman": round(scores["model"], 4),
            "gc_only": round(scores["gc_only"], 4),
            "gene_mean": round(scores["gene_mean"], 4),
            "ships": passes, "blocked_by": "; ".join(reasons) or None,
        })

        mark = "ship " if passes else "skip "
        print(f"  {mark} {slug:34s} n={len(frame):6d} genes={rows[-1]['genes']:3d} "
              f"rho={scores['model']:+.3f} (gc {scores['gc_only']:+.3f}, "
              f"gene {scores['gene_mean']:+.3f})")
        if reasons:
            print(f"         blocked: {'; '.join(reasons)}")

        if passes:
            fitted = _model().fit(
                scores["frame_features"][scores["features"]].astype(float), frame["y"].values
            )
            joblib.dump(fitted, MODELS / f"potency_{slug}.joblib")
            shipped[slug] = {
                "model_file": f"potency_{slug}.joblib",
                "modality": "ASO",
                "readout": "percent knockdown",
                "units": f"percent knockdown in {meta['cell_line']}",
                "context": meta,
                "feature_columns": scores["features"],
                "sequences": int(len(frame)),
                "genes": int(frame["target_gene"].nunique()),
                "spearman": round(scores["model"], 4),
                "evidence": (
                    f"Spearman {scores['model']:.2f} under target-gene-grouped CV on "
                    f"{len(frame):,} sequences in {meta['cell_line']}, against "
                    f"{scores['gc_only']:.2f} for GC alone and {scores['gene_mean']:.2f} "
                    "for a target-gene-mean baseline that ignores the sequence."
                ),
                "caveat": (
                    "Predicts knockdown in this assay only. ASO Atlas pools percent "
                    "inhibition across 81 cell lines and many dosages, and the pooled "
                    "model scores 0.18 against this context's "
                    f"{scores['model']:.2f}. Do not read it as general potency."
                ),
            }

    (CONFIG / "potency_contexts.json").write_text(json.dumps({
        "source": "ASO Atlas, https://github.com/barneyhill/aso_atlas (Apache-2.0)",
        "why_per_context": "Pooled, the corpus scores below the model it would replace. "
                           "The sequence-potency relationship differs between assays, and "
                           "adding context as features made it worse, so each assay gets "
                           "its own model.",
        "gates": {"min_rows": MIN_ROWS, "min_genes": MIN_GENES, "floor_spearman": FLOOR,
                  "must_beat": ["target gene mean", "GC only"]},
        "contexts": shipped,
    }, indent=2) + "\n")

    RESULTS.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(RESULTS / "context_potency.csv", index=False)
    print(f"\nshipped {len(shipped)} of {len(contexts)} contexts")
    print(f"wrote {(CONFIG / 'potency_contexts.json').relative_to(ROOT)}")


if __name__ == "__main__":
    main()
