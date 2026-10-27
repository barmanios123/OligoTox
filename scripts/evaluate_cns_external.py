"""Build and test a CNS neurotoxicity endpoint from two independent assays.

The shipped cns_toxicity model rests on 2 positive sequences and is the one
endpoint training refuses to evaluate. Two external datasets in
data/oligogym/ measure CNS neurotoxicity on thousands of sequences:

  MOE_Neurotox_2024  2,437 rows / 1,685 sequences, in vivo, FOB score at 3h
  Hagedorn_2022      1,825 rows / 1,653 sequences, in vitro, calcium
                     oscillation score, LNA gapmers

They share no sequences, so training on one and testing on the other is a
real external test, not a resplit.

Label direction matters and is not symmetric between the assays:

  FOB score      higher = more functional abnormality = toxic.
  Calcium score  LOWER = toxic. Hagedorn et al. (Nucleic Acid Ther. 2022)
                 report that acute tolerability correlates with REDUCTIONS in
                 spontaneous calcium oscillations. The `y_binary` column in
                 data/oligogym/oligogym_tox_combined_fixed.csv has this
                 inverted and calls 83% of the compounds toxic.

Neither threshold is given by the data, so each is a stated choice and both
are swept in the sensitivity section.

Baselines include the two sequence drivers Hagedorn et al. name explicitly,
G count and the distance from the 3' end to the nearest G. If the full model
cannot beat those, the model has learned nothing beyond the published rule.

Writes reports/results/cns_external.json and cns_external.csv.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.calibration import CalibratedClassifierCV
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict

from oligotox.features.candidate import sequence_features
from oligotox.features.structure_features import build_structure_features

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data" / "oligogym"
RESULTS = ROOT / "reports" / "results"

COMBINED = DATA / "oligogym_tox_combined.csv"
MOE = "MOE_Neurotox_2024"
HAGEDORN = "Hagedorn_2022"

# Stated label thresholds. Swept in sensitivity().
FOB_TOXIC_AT_OR_ABOVE = 3.0
CALCIUM_TOXIC_BELOW = 50.0
N_BOOT = 1000
SEED = 42


def _model():
    base = LogisticRegression(max_iter=2000, class_weight="balanced", solver="liblinear")
    return CalibratedClassifierCV(base, method="sigmoid", cv=3)


def _normalize(series: pd.Series) -> pd.Series:
    return series.astype(str).str.upper().str.replace("U", "T", regex=False).str.strip()


def g_features(seq: str) -> dict:
    """The two drivers Hagedorn et al. name, plus the G run length."""
    n = len(seq)
    g_count = seq.count("G")
    # Distance from the 3' end to the nearest G; n when the sequence has none.
    nearest = next((i for i, ch in enumerate(reversed(seq)) if ch == "G"), n)
    longest_g = max((len(r) for r in "".join(c if c == "G" else " " for c in seq).split()), default=0)
    return {
        "g_count": float(g_count),
        "g_fraction": g_count / n if n else 0.0,
        "dist_3p_to_g": float(nearest),
        "longest_g_run": float(longest_g),
    }


def build(source: str) -> pd.DataFrame:
    """One row per unique sequence, with features and a binary label."""
    raw = pd.read_csv(COMBINED, low_memory=False)
    df = raw[raw.source == source].copy()
    df["seq"] = _normalize(df["sequence"])
    df["y_raw"] = pd.to_numeric(df["y_raw"], errors="coerce")
    df = df.dropna(subset=["y_raw"])

    # Replicates of a sequence are averaged, so no sequence appears twice.
    agg = df.groupby("seq", as_index=False)["y_raw"].mean()

    if source == MOE:
        agg["toxic"] = (agg["y_raw"] >= FOB_TOXIC_AT_OR_ABOVE).astype(int)
    else:
        agg["toxic"] = (agg["y_raw"] < CALCIUM_TOXIC_BELOW).astype(int)

    rows = []
    for seq in agg["seq"]:
        feats = sequence_features(seq)
        feats.pop("ps_x_len", None)   # chemistry is unknown for these sets
        feats.pop("lna_x_gc", None)
        feats.update(g_features(seq))
        rows.append(feats)
    feat = pd.DataFrame(rows, index=agg.index)

    struct = build_structure_features(pd.DataFrame({"sequence": agg["seq"]}))
    struct = struct.loc[:, [c for c in struct.columns if struct[c].nunique() > 1]]
    out = pd.concat([agg.reset_index(drop=True), feat.reset_index(drop=True),
                     struct.reset_index(drop=True)], axis=1)
    return out.loc[:, ~out.columns.duplicated()]


def _feature_cols(df: pd.DataFrame) -> list[str]:
    drop = {"seq", "y_raw", "toxic"}
    return [c for c in df.columns if c not in drop and pd.api.types.is_numeric_dtype(df[c])]


def _boot_auc(y, p, rng, n=N_BOOT):
    scores = []
    idx = np.arange(len(y))
    for _ in range(n):
        s = rng.choice(idx, size=len(idx), replace=True)
        if len(np.unique(y[s])) < 2:
            continue
        scores.append(roc_auc_score(y[s], p[s]))
    if not scores:
        return (float("nan"), float("nan"))
    return tuple(float(v) for v in np.percentile(scores, [2.5, 97.5]))


def within_assay(df: pd.DataFrame, name: str, rng) -> list[dict]:
    """Does either assay contain learnable sequence signal at all?"""
    y = df["toxic"].values
    cols = _feature_cols(df)
    sets = {
        "All sequence features": cols,
        "G features only (Hagedorn rule)": ["g_count", "g_fraction", "dist_3p_to_g", "longest_g_run"],
        "GC only": ["gc_content"],
        "Length only": ["seq_length"],
    }
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
    out = []
    for label, feats in sets.items():
        X = df[feats].astype(float)
        prob = cross_val_predict(_model(), X, y, cv=cv, method="predict_proba")[:, 1]
        lo, hi = _boot_auc(y, prob, rng)
        out.append({
            "assay": name, "evaluation": "within-assay 5-fold CV", "features": label,
            "roc_auc": float(roc_auc_score(y, prob)), "ci_low": lo, "ci_high": hi,
            "pr_auc": float(average_precision_score(y, prob)),
            "base_rate": float(y.mean()), "n": int(len(y)), "n_positive": int(y.sum()),
        })
    return out


def cross_assay(train: pd.DataFrame, test: pd.DataFrame, a: str, b: str, rng) -> list[dict]:
    """Train on one assay, predict the other. No shared sequences."""
    shared = sorted(set(_feature_cols(train)) & set(_feature_cols(test)))
    out = []
    sets = {
        "All sequence features": shared,
        "G features only (Hagedorn rule)": [c for c in
                                            ["g_count", "g_fraction", "dist_3p_to_g", "longest_g_run"]
                                            if c in shared],
        "GC only": ["gc_content"],
    }
    for label, feats in sets.items():
        model = _model().fit(train[feats].astype(float), train["toxic"].values)
        prob = model.predict_proba(test[feats].astype(float))[:, 1]
        y = test["toxic"].values
        lo, hi = _boot_auc(y, prob, rng)
        out.append({
            "assay": f"{a} -> {b}", "evaluation": "cross-assay transfer", "features": label,
            "roc_auc": float(roc_auc_score(y, prob)), "ci_low": lo, "ci_high": hi,
            "pr_auc": float(average_precision_score(y, prob)),
            "base_rate": float(y.mean()), "n": int(len(y)), "n_positive": int(y.sum()),
        })
    return out


def sensitivity(moe: pd.DataFrame, hag: pd.DataFrame) -> list[dict]:
    """Both thresholds are choices, so vary them and see what moves."""
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
    out = []
    for name, df, col, thresholds, direction in [
        (MOE, moe, "y_raw", [2.0, 2.5, 3.0, 3.5, 4.0], "above"),
        (HAGEDORN, hag, "y_raw", [25.0, 35.0, 50.0, 62.0, 75.0], "below"),
    ]:
        cols = _feature_cols(df)
        for t in thresholds:
            y = ((df[col] >= t) if direction == "above" else (df[col] < t)).astype(int).values
            if y.sum() < 20 or (1 - y).sum() < 20:
                continue
            prob = cross_val_predict(_model(), df[cols].astype(float), y, cv=cv,
                                     method="predict_proba")[:, 1]
            out.append({
                "assay": name, "evaluation": "threshold sensitivity",
                "features": f"toxic if y {'>=' if direction == 'above' else '<'} {t}",
                "roc_auc": float(roc_auc_score(y, prob)), "ci_low": None, "ci_high": None,
                "pr_auc": float(average_precision_score(y, prob)),
                "base_rate": float(y.mean()), "n": int(len(y)), "n_positive": int(y.sum()),
            })
    return out


def main() -> None:
    rng = np.random.default_rng(SEED)
    moe, hag = build(MOE), build(HAGEDORN)

    overlap = set(moe["seq"]) & set(hag["seq"])
    assert not overlap, f"assays share {len(overlap)} sequences; transfer would not be external"

    rows = []
    rows += within_assay(moe, MOE, rng)
    rows += within_assay(hag, HAGEDORN, rng)
    rows += cross_assay(moe, hag, MOE, HAGEDORN, rng)
    rows += cross_assay(hag, moe, HAGEDORN, MOE, rng)
    rows += sensitivity(moe, hag)

    df = pd.DataFrame(rows)
    RESULTS.mkdir(parents=True, exist_ok=True)
    df.to_csv(RESULTS / "cns_external.csv", index=False)

    summary = {
        "datasets": {
            MOE: {"sequences": int(len(moe)), "positives": int(moe["toxic"].sum()),
                  "readout": "FOB_score_3h", "assay": "in vivo",
                  "label": f"toxic if FOB >= {FOB_TOXIC_AT_OR_ABOVE}"},
            HAGEDORN: {"sequences": int(len(hag)), "positives": int(hag["toxic"].sum()),
                       "readout": "calcium_oscillation_score", "assay": "in vitro",
                       "label": f"toxic if calcium < {CALCIUM_TOXIC_BELOW}"},
        },
        "shared_sequences": 0,
        "results": rows,
    }
    (RESULTS / "cns_external.json").write_text(json.dumps(summary, indent=2) + "\n")

    show = df[df.evaluation != "threshold sensitivity"]
    print(show.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    print(f"\nwrote {(RESULTS / 'cns_external.csv').relative_to(ROOT)} and .json")


if __name__ == "__main__":
    main()
