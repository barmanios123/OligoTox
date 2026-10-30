"""Compare the shipped model against simple baselines, with bootstrap CIs.

Every model is scored with the same sequence-grouped folds used in training.
The "data source only" baseline sees nothing but which dataset a row came
from; if it matches the real model, the model is mostly recognising sources.

Writes reports/results/baselines.csv and reports/results/within_source.csv.
Needs the processed training data in data/processed/.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import average_precision_score
from sklearn.model_selection import StratifiedGroupKFold, cross_val_predict

from oligotox.training.train import _make_model

ROOT = Path(__file__).resolve().parent.parent
PROCESSED = ROOT / "data" / "processed"
RESULTS = ROOT / "reports" / "results"
MAIN_SOURCES = ("oligotox_main", "papargyri_2020")


def _label(n: int, kind: str) -> str:
    return f"{kind} ({n} features)"
N_BOOT = 1000


def _boosted():
    return HistGradientBoostingClassifier(
        max_iter=200, learning_rate=0.05, max_leaf_nodes=15,
        class_weight="balanced", random_state=42,
    )


def _oof(model, X, y, groups, n_splits):
    cv = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=42)
    return cross_val_predict(model, X, y, cv=cv.split(X, y, groups), method="predict_proba")[:, 1]


def _bootstrap_ci(y, prob, groups, rng):
    """95% CI for PR-AUC, resampling whole sequences rather than rows."""
    uniq = np.unique(groups)
    idx_by_group = {g: np.flatnonzero(groups == g) for g in uniq}
    scores = []
    for _ in range(N_BOOT):
        idx = np.concatenate([idx_by_group[g] for g in rng.choice(uniq, len(uniq))])
        if y[idx].sum() == 0 or y[idx].sum() == len(idx):
            continue
        scores.append(average_precision_score(y[idx], prob[idx]))
    return np.percentile(scores, [2.5, 97.5])


def main() -> None:
    X = pd.read_csv(PROCESSED / "feature_matrix.csv")
    Y = pd.read_csv(PROCESSED / "label_matrix.csv")
    raw = pd.read_csv(PROCESSED / "training_dataset.csv", low_memory=False)
    summary = json.loads((RESULTS / "training_summary_after.json").read_text())

    groups = raw["sequence"].astype(str).str.upper().str.replace("U", "T", regex=False).values
    source = raw["source"].where(raw["source"].isin(MAIN_SOURCES), "curated_other")
    source_X = pd.get_dummies(source, dtype=float)

    n_feat = X.shape[1]
    LOGISTIC, BOOSTED = _label(n_feat, "Logistic"), _label(n_feat, "Gradient boosting")
    models = {
        LOGISTIC: (_make_model, X),
        BOOSTED: (_boosted, X),
        "GC content only": (_make_model, X[["gc_content"]]),
        "Length only": (_make_model, X[["seq_length"]]),
        "Data source only": (_make_model, source_X),
    }

    rng = np.random.default_rng(42)
    rows, within = [], []
    for endpoint, s in summary.items():
        if s["status"] != "ok":
            continue
        # Only the rows measured for this endpoint, matching the trainer.
        keep = pd.to_numeric(Y[endpoint], errors="coerce").notna().values
        y = pd.to_numeric(Y[endpoint][keep], errors="coerce").astype(int).values
        groups_e, source_e = groups[keep], source[keep].reset_index(drop=True)
        # Every baseline is scored on the same rows as the model it is a
        # baseline for, or the comparison means nothing.
        models_e = {name: (factory, feats[keep].reset_index(drop=True))
                    for name, (factory, feats) in models.items()}
        n_splits = s["n_splits"]
        for name, (factory, feats) in models_e.items():
            prob = _oof(factory(), feats, y, groups_e, n_splits)
            lo, hi = _bootstrap_ci(y, prob, groups_e, rng)
            rows.append({
                "endpoint": endpoint, "model": name,
                "pr_auc": average_precision_score(y, prob), "ci_low": lo, "ci_high": hi,
                "base_rate": y.mean(), "positive_sequences": s["positive_sequences"],
            })

        # Within the source that holds most positives, the source shortcut is gone.
        top = source_e[y == 1].value_counts()
        mask = (source_e == top.index[0]).values
        ys, gs = y[mask], groups_e[mask]
        pos_seq = len(set(gs[ys == 1]))
        neg_seq = len(set(gs[ys == 0]))
        if pos_seq >= 3 and neg_seq >= 3:
            for name in (LOGISTIC, "GC content only", "Length only"):
                # models_e, not models: `mask` indexes this endpoint's rows.
                factory, feats = models_e[name]
                prob = _oof(factory(), feats[mask].reset_index(drop=True), ys, gs,
                            min(5, pos_seq, neg_seq))
                lo, hi = _bootstrap_ci(ys, prob, gs, rng)
                within.append({
                    "endpoint": endpoint, "source": top.index[0], "model": name,
                    "pr_auc": average_precision_score(ys, prob), "ci_low": lo, "ci_high": hi,
                    "base_rate": ys.mean(), "positive_sequences": pos_seq,
                    "negative_sequences": neg_seq,
                })

    pd.DataFrame(rows).to_csv(RESULTS / "baselines.csv", index=False, float_format="%.4f")
    pd.DataFrame(within).to_csv(RESULTS / "within_source.csv", index=False, float_format="%.4f")
    print(pd.DataFrame(rows).pivot(index="endpoint", columns="model", values="pr_auc").round(2).to_string())
    print()
    print(pd.DataFrame(within).round(2).to_string(index=False))


if __name__ == "__main__":
    main()
