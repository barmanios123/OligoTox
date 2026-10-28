"""Does ASO Atlas improve on the shipped ASO activity model?

The shipped model is trained on McQuisten 2007: 3,867 sequences, Spearman
0.26 under target-gene-grouped CV, tiered `weak`. ASO Atlas is 190,927 rows
over 165,555 unique sequences and 343 genes, with chemistry annotation no
other source here carries.

More data is not automatically better, so this applies the same rules the
rest of the project uses:

  grouped by target gene   potency varies far more between genes than within
                           one, and a model that has seen other ASOs against
                           the same gene can score well knowing nothing about
                           design
  gene-mean baseline       predicts each gene's mean and ignores the sequence
  simple baselines         GC alone, length alone
  held-out study           scored against McQuisten, which shares no rows, to
                           check it transfers rather than just fits

Readout is percent inhibition, higher being more potent. Values outside
0-100 appear (down to -786 and up to 224) and are clipped, which is stated
rather than silent because clipping is a modelling choice.

Writes reports/results/aso_atlas.csv and .json.
"""
from __future__ import annotations

import json
import re
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
sys.path.insert(0, str(ROOT / "scripts"))

from oligotox.efficacy.potency import guide_features  # noqa: E402

ATLAS = ROOT / "data" / "external" / "aso_atlas.parquet"
RESULTS = ROOT / "reports" / "results"
CLIP = (0.0, 100.0)
MIN_PER_GENE = 20
SEED = 42
# Modifications named in the chemistry field, which arrives as a flattened stub.
CHEMISTRIES = ("MOE", "PS", "cEt", "LNA", "PO", "OMe")


def _model():
    return make_pipeline(StandardScaler(), Ridge(alpha=1.0))


def chemistry_flags(blob: str) -> dict:
    text = blob or ""
    flags = {f"chem_{c.lower()}": float(bool(re.search(rf"'{c}'", text))) for c in CHEMISTRIES}
    length = re.search(r"'length':\s*(\d+)", text)
    flags["chem_declared_length"] = float(length.group(1)) if length else 0.0
    flags["chem_n_modifications"] = float(text.count("'modification':"))
    return flags


def build() -> pd.DataFrame:
    raw = pd.read_parquet(ATLAS)
    raw = raw.dropna(subset=["aso_sequence_5_to_3", "inhibition_percent", "target_gene"])
    raw["sequence"] = raw["aso_sequence_5_to_3"].astype(str).str.upper().str.strip()
    raw = raw[raw["sequence"].str.fullmatch(r"[ACGT]{12,30}")]
    raw["y"] = raw["inhibition_percent"].clip(*CLIP)

    # One row per sequence-gene pair; replicates averaged so none spans a split.
    agg = raw.groupby(["sequence", "target_gene"], as_index=False).agg(
        y=("y", "mean"), n_replicates=("y", "size"),
        chemistry=("chemistry", "first"),
    )
    # Genes with too few sequences cannot support a held-out fold.
    counts = agg["target_gene"].value_counts()
    agg = agg[agg["target_gene"].isin(counts[counts >= MIN_PER_GENE].index)]

    feats = pd.DataFrame([guide_features(s) for s in agg["sequence"]], index=agg.index)
    chem = pd.DataFrame([chemistry_flags(c) for c in agg["chemistry"]], index=agg.index)
    out = pd.concat([agg[["sequence", "target_gene", "y"]], feats, chem], axis=1)
    return out.reset_index(drop=True)


def _features(df: pd.DataFrame) -> list[str]:
    drop = {"sequence", "target_gene", "y"}
    return [c for c in df.columns if c not in drop and pd.api.types.is_numeric_dtype(df[c])]


def _gene_mean_oof(df, cv, groups):
    pred = np.full(len(df), np.nan)
    y = df["y"].values
    for train_idx, test_idx in cv.split(df, y, groups):
        pred[test_idx] = y[train_idx].mean()
    return pred


def main() -> None:
    if not ATLAS.exists():
        raise SystemExit(
            f"missing {ATLAS}. Convert the pickle first:\n"
            "  python scripts/ingest_aso_atlas.py path/to/aso_atlas.pkl"
        )

    df = build()
    print(f"{len(df):,} sequence-gene pairs over {df.target_gene.nunique()} genes")

    y = df["y"].values
    genes = df["target_gene"].values
    cv = GroupKFold(n_splits=5)
    cols = _features(df)
    sequence_only = [c for c in cols if not c.startswith("chem_")]

    rows = []
    sets = {
        "All features (sequence + chemistry)": cols,
        "Sequence features only": sequence_only,
        "Chemistry only": [c for c in cols if c.startswith("chem_")],
        "GC only": ["gc_content"],
    }
    for label, feats in sets.items():
        if not feats:
            continue
        pred = cross_val_predict(_model(), df[feats].astype(float), y, cv=cv, groups=genes)
        rho = float(spearmanr(y, pred).statistic)
        rows.append({"dataset": "ASO Atlas", "evaluation": "gene-grouped CV",
                     "model": label, "n_features": len(feats), "spearman": rho,
                     "n": int(len(y)), "n_genes": int(df.target_gene.nunique())})
        print(f"  {label:38s} {len(feats):3d} feat  rho={rho:.3f}")

    pred = _gene_mean_oof(df, cv, genes)
    rho = spearmanr(y, pred).statistic
    rows.append({"dataset": "ASO Atlas", "evaluation": "gene-grouped CV",
                 "model": "Target gene mean (no sequence)", "n_features": 0,
                 "spearman": float(rho) if np.isfinite(rho) else 0.0,
                 "n": int(len(y)), "n_genes": int(df.target_gene.nunique())})
    print(f"  {'Target gene mean (no sequence)':38s}   0 feat  rho={rows[-1]['spearman']:.3f}")

    # Transfer to the study the shipped model was trained on.
    from evaluate_efficacy import ASO, build as build_mcquisten

    mcq = build_mcquisten(ASO)
    shared = sorted(set(sequence_only) & set(mcq.columns))
    fitted = _model().fit(df[shared].astype(float), y)
    transfer = spearmanr(mcq["y_raw"].values,
                         fitted.predict(mcq[shared].astype(float))).statistic
    rows.append({"dataset": "ASO Atlas -> McQuisten_2007", "evaluation": "cross-study transfer",
                 "model": "Sequence features only", "n_features": len(shared),
                 "spearman": float(transfer), "n": int(len(mcq)),
                 "n_genes": int(mcq.target_gene.nunique())})
    print(f"\n  transfer to McQuisten 2007: rho={transfer:.3f} "
          f"(shipped model scores 0.26 on its own data)")

    RESULTS.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(RESULTS / "aso_atlas.csv", index=False)
    (RESULTS / "aso_atlas.json").write_text(json.dumps({
        "source": "https://github.com/barneyhill/aso_atlas (Apache-2.0)",
        "clip": list(CLIP),
        "min_sequences_per_gene": MIN_PER_GENE,
        "note": "Percent inhibition, higher is more potent. Values outside 0-100 "
                "are clipped, which is a stated modelling choice.",
        "results": rows,
    }, indent=2) + "\n")
    print(f"\nwrote {(RESULTS / 'aso_atlas.csv').relative_to(ROOT)}")


if __name__ == "__main__":
    main()
