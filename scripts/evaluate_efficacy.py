"""Can OligoTox predict potency, and is it predicting the oligo or the target?

6,718 potency sequences sit in data/oligogym unused. This asks whether they
support an efficacy model, under the same rules the toxicity side now
follows: grouped splits, a baseline that could explain the result away, and
confidence intervals.

The baseline that matters here is the target gene. Potency varies far more
between genes than within one, so a model that has seen other siRNAs against
the same gene can score well while knowing nothing about oligo design. Two
guards:

  * folds are grouped by target gene, so a gene never spans a split;
  * a "target gene mean" baseline predicts each gene's mean potency from the
    training folds and ignores the sequence entirely.

If the sequence model cannot beat that baseline, it has learned which gene is
easy to silence, which is the potency equivalent of the data-source shortcut
found on the toxicity side.

Spearman is the primary metric, so no potency threshold has to be invented.

Datasets, from data/processed/external_efficacy_normalized.csv:
  Ichihara_2007_1  2,431 siRNA duplexes, 30 genes, percent inhibition
  Ichihara_2007_2    419 siRNA duplexes, 11 genes, percent inhibition
  McQuisten_2007   3,868 PS-ASOs,        85 genes, activity 0-1

Writes reports/results/efficacy.csv and efficacy.json.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from oligotox.features.candidate import g_motif_features, sequence_features

ROOT = Path(__file__).resolve().parent.parent
PROCESSED = ROOT / "data" / "processed"
RESULTS = ROOT / "reports" / "results"

SIRNA = ("Ichihara_2007_1", "Ichihara_2007_2")
ASO = "McQuisten_2007"
N_BOOT = 1000
SEED = 42


def _model():
    return make_pipeline(StandardScaler(), Ridge(alpha=1.0))


def guide_features(seq: str) -> dict:
    """Design features of the guide strand, plus siRNA-specific ones."""
    seq = seq.strip().upper().replace("U", "T")
    feats = sequence_features(seq)
    feats.pop("ps_x_len", None)
    feats.pop("lna_x_gc", None)
    feats.update(g_motif_features(seq))

    n = len(seq)
    # Thermodynamic asymmetry: an AU-rich 5' end of the guide favours loading
    # into RISC, which is the best known sequence determinant of siRNA potency.
    five = seq[:5]
    three = seq[-5:]
    feats["five_prime_au"] = sum(c in "AT" for c in five) / max(len(five), 1)
    feats["three_prime_au"] = sum(c in "AT" for c in three) / max(len(three), 1)
    feats["asymmetry"] = feats["five_prime_au"] - feats["three_prime_au"]
    for base in "ACGT":
        feats[f"frac_{base}"] = seq.count(base) / n if n else 0.0
    # Position 1 and the seed, both standard potency determinants.
    for base in "ACGT":
        feats[f"pos1_{base}"] = float(seq[:1] == base)
    return feats


def build(source: str) -> pd.DataFrame:
    df = pd.read_csv(PROCESSED / "external_efficacy_normalized.csv")
    df = df[df.source == source].dropna(subset=["target_gene", "y_raw"]).copy()
    feats = pd.DataFrame([guide_features(s) for s in df["guide_sequence"]], index=df.index)
    out = pd.concat([df[["guide_sequence", "y_raw", "target_gene"]], feats], axis=1)
    return out.reset_index(drop=True)


def _feature_cols(df: pd.DataFrame) -> list[str]:
    drop = {"guide_sequence", "y_raw", "target_gene"}
    return [c for c in df.columns if c not in drop and pd.api.types.is_numeric_dtype(df[c])]


def _gene_mean_oof(df: pd.DataFrame, cv, groups) -> np.ndarray:
    """Predict each gene's mean potency, learned only from training folds."""
    pred = np.full(len(df), np.nan)
    y = df["y_raw"].values
    for train_idx, test_idx in cv.split(df, y, groups):
        # A held-out gene is unseen by construction, so fall back to the
        # global training mean. That is the honest version of this baseline.
        pred[test_idx] = y[train_idx].mean()
    return pred


def _boot_spearman(y, p, rng, n=N_BOOT):
    idx = np.arange(len(y))
    vals = []
    for _ in range(n):
        s = rng.choice(idx, size=len(idx), replace=True)
        if np.std(p[s]) == 0 or np.std(y[s]) == 0:
            continue
        vals.append(spearmanr(y[s], p[s]).statistic)
    if not vals:
        return (float("nan"), float("nan"))
    return tuple(float(v) for v in np.percentile(vals, [2.5, 97.5]))


def evaluate(df: pd.DataFrame, name: str, rng) -> list[dict]:
    y = df["y_raw"].values
    genes = df["target_gene"].values
    n_splits = min(5, df["target_gene"].nunique())
    cv = GroupKFold(n_splits=n_splits)
    cols = _feature_cols(df)

    sets = {
        "All guide features": cols,
        "Asymmetry + seed only": ["five_prime_au", "three_prime_au", "asymmetry", "seed_gc"],
        "GC only": ["gc_content"],
    }
    rows = []
    for label, feats in sets.items():
        pred = cross_val_predict(_model(), df[feats].astype(float), y, cv=cv, groups=genes)
        rho = spearmanr(y, pred).statistic
        lo, hi = _boot_spearman(y, pred, rng)
        rows.append({"dataset": name, "evaluation": "gene-grouped CV", "model": label,
                     "spearman": float(rho), "ci_low": lo, "ci_high": hi,
                     "n": int(len(y)), "n_genes": int(df.target_gene.nunique())})

    pred = _gene_mean_oof(df, cv, genes)
    rho = spearmanr(y, pred).statistic
    rows.append({"dataset": name, "evaluation": "gene-grouped CV",
                 "model": "Target gene mean (no sequence)",
                 "spearman": float(rho) if np.isfinite(rho) else 0.0,
                 "ci_low": None, "ci_high": None,
                 "n": int(len(y)), "n_genes": int(df.target_gene.nunique())})

    # How much of the variance is between genes rather than within them?
    overall = y.mean()
    between = df.groupby("target_gene")["y_raw"].transform("mean")
    rows.append({"dataset": name, "evaluation": "variance decomposition",
                 "model": "fraction of variance between genes",
                 "spearman": float(((between - overall) ** 2).sum() / ((y - overall) ** 2).sum()),
                 "ci_low": None, "ci_high": None,
                 "n": int(len(y)), "n_genes": int(df.target_gene.nunique())})
    return rows


def transfer(train: pd.DataFrame, test: pd.DataFrame, a: str, b: str, rng) -> list[dict]:
    """Train on one siRNA study, predict the other. Ranks only, scales differ."""
    shared = sorted(set(_feature_cols(train)) & set(_feature_cols(test)))
    rows = []
    for label, feats in {"All guide features": shared, "GC only": ["gc_content"]}.items():
        model = _model().fit(train[feats].astype(float), train["y_raw"].values)
        pred = model.predict(test[feats].astype(float))
        y = test["y_raw"].values
        lo, hi = _boot_spearman(y, pred, rng)
        rows.append({"dataset": f"{a} -> {b}", "evaluation": "cross-study transfer",
                     "model": label, "spearman": float(spearmanr(y, pred).statistic),
                     "ci_low": lo, "ci_high": hi, "n": int(len(y)),
                     "n_genes": int(test.target_gene.nunique())})
    return rows


def main() -> None:
    rng = np.random.default_rng(SEED)
    frames = {s: build(s) for s in (*SIRNA, ASO)}

    rows = []
    for name, df in frames.items():
        rows += evaluate(df, name, rng)
    rows += transfer(frames[SIRNA[0]], frames[SIRNA[1]], *SIRNA, rng=rng)
    rows += transfer(frames[SIRNA[1]], frames[SIRNA[0]], SIRNA[1], SIRNA[0], rng=rng)

    df = pd.DataFrame(rows)
    RESULTS.mkdir(parents=True, exist_ok=True)
    df.to_csv(RESULTS / "efficacy.csv", index=False)
    (RESULTS / "efficacy.json").write_text(json.dumps({
        "note": "Spearman between measured potency and out-of-fold prediction. "
                "Folds are grouped by target gene. The target-gene-mean baseline "
                "uses no sequence information.",
        "datasets": {k: {"sequences": int(len(v)), "genes": int(v.target_gene.nunique())}
                     for k, v in frames.items()},
        "results": rows,
    }, indent=2) + "\n")

    print(df.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    print(f"\nwrote {(RESULTS / 'efficacy.csv').relative_to(ROOT)}")


if __name__ == "__main__":
    main()
