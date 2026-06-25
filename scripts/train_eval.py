from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    average_precision_score,
    brier_score_loss,
    roc_auc_score,
)
from sklearn.calibration import CalibratedClassifierCV
from sklearn.model_selection import StratifiedShuffleSplit, GroupShuffleSplit
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data" / "processed"
OUT_DIR = ROOT / "artifacts" / "eval"
OUT_DIR.mkdir(parents=True, exist_ok=True)

TRAIN_PATH = DATA_DIR / "training_dataset.csv"
FEATURE_PATH = DATA_DIR / "feature_matrix.csv"
LABEL_PATH = DATA_DIR / "label_matrix.csv"
META_PATH = DATA_DIR / "dataset_metadata.json"


def load_inputs():
    X = pd.read_csv(FEATURE_PATH)
    Y = pd.read_csv(LABEL_PATH)
    train = pd.read_csv(TRAIN_PATH)
    meta = json.loads(META_PATH.read_text())
    return X, Y, train, meta


def build_groups(train: pd.DataFrame) -> pd.Series:
    source = train.get("source", pd.Series(index=train.index, data="unknown")).astype(str).fillna("unknown")
    series = train.get("series_id", pd.Series(index=train.index, data="no_series")).astype(str).fillna("no_series")
    seq_len = train.get("seq_length", pd.Series(index=train.index, data=-1)).fillna(-1).astype(int)
    ps_frac = train.get("ps_fraction", pd.Series(index=train.index, data=-1.0)).fillna(-1.0).astype(float)
    ps_bin = np.round(ps_frac, 2)

    groups = (
        source
        + "||"
        + series
        + "||"
        + seq_len.astype(str)
        + "||"
        + ps_bin.astype(str)
    )
    groups = groups.fillna("unknown_group")
    groups = groups.replace({"nan": "unknown_group"})
    return groups


def make_model():
    return Pipeline(
        steps=[
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
            ("clf", LogisticRegression(max_iter=2000, class_weight="balanced")),
        ]
    )


def fit_calibrated(X_train, y_train):
    base = make_model()
    calib = CalibratedClassifierCV(base, method="sigmoid", cv=3)
    calib.fit(X_train, y_train)
    return calib


def safe_metric(fn, y_true, y_score):
    try:
        if len(np.unique(y_true)) < 2:
            return np.nan
        return fn(y_true, y_score)
    except Exception:
        return np.nan


def evaluate_endpoint(name, X, y, train, groups):
    rows = []

    y = pd.Series(y).fillna(0).astype(int)

    # Global support metrics for this endpoint
    total_pos = int(y.sum())
    total_neg = int(len(y) - y.sum())
    total_groups = int(pd.Series(groups).nunique())
    groups_with_pos = int(pd.Series(groups)[y == 1].nunique()) if total_pos > 0 else 0

    if total_pos < 5 or total_neg < 5:
        for split_type in ["random", "grouped"]:
            rows.append(
                {
                    "endpoint": name,
                    "split_type": split_type,
                    "n_train": 0,
                    "n_test": 0,
                    "positive_rate_test": np.nan,
                    "pr_auc": np.nan,
                    "roc_auc": np.nan,
                    "brier": np.nan,
                    "total_positives": total_pos,
                    "total_negatives": total_neg,
                    "total_groups": total_groups,
                    "groups_with_positives": groups_with_pos,
                    "status": "skipped_low_support",
                }
            )
        return rows

    splitters = {
        "random": StratifiedShuffleSplit(n_splits=1, test_size=0.2, random_state=42),
        "grouped": GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=42),
    }

    for split_type, splitter in splitters.items():
        if split_type == "random":
            train_idx, test_idx = next(splitter.split(X, y))
        else:
            train_idx, test_idx = next(splitter.split(X, y, groups=groups))

        X_train = X.iloc[train_idx]
        X_test = X.iloc[test_idx]
        y_train = y.iloc[train_idx]
        y_test = y.iloc[test_idx]

        if y_train.nunique() < 2 or y_test.nunique() < 2:
            rows.append(
                {
                    "endpoint": name,
                    "split_type": split_type,
                    "n_train": len(train_idx),
                    "n_test": len(test_idx),
                    "positive_rate_test": float(y_test.mean()) if len(y_test) else np.nan,
                    "pr_auc": np.nan,
                    "roc_auc": np.nan,
                    "brier": np.nan,
                    "total_positives": total_pos,
                    "total_negatives": total_neg,
                    "total_groups": total_groups,
                    "groups_with_positives": groups_with_pos,
                    "status": "skipped_single_class_fold",
                }
            )
            continue

        model = fit_calibrated(X_train, y_train)
        y_score = model.predict_proba(X_test)[:, 1]

        rows.append(
            {
                "endpoint": name,
                "split_type": split_type,
                "n_train": len(train_idx),
                "n_test": len(test_idx),
                "positive_rate_test": float(y_test.mean()),
                "pr_auc": safe_metric(average_precision_score, y_test, y_score),
                "roc_auc": safe_metric(roc_auc_score, y_test, y_score),
                "brier": safe_metric(brier_score_loss, y_test, y_score),
                "total_positives": total_pos,
                "total_negatives": total_neg,
                "total_groups": total_groups,
                "groups_with_positives": groups_with_pos,
                "status": "ok",
            }
        )

    return rows


def main():
    X, Y, train, meta = load_inputs()
    groups = build_groups(train)

    common_n = min(len(X), len(Y), len(train), len(groups))
    X = X.iloc[:common_n].copy()
    Y = Y.iloc[:common_n].copy()
    train = train.iloc[:common_n].copy()
    groups = groups.iloc[:common_n].copy()

    assert not pd.isna(groups).any(), "groups contains NaN after build_groups"

    results = []
    for endpoint in Y.columns:
        results.extend(evaluate_endpoint(endpoint, X, Y[endpoint], train, groups))

    df = pd.DataFrame(results)
    df.to_csv(OUT_DIR / "train_eval_report.csv", index=False)

    pivot = df.pivot_table(
        index="endpoint",
        columns="split_type",
        values=["pr_auc", "roc_auc", "brier"],
        aggfunc="first",
    )
    pivot.columns = ["_".join(col).strip() for col in pivot.columns.values]
    pivot = pivot.reset_index()

    if "pr_auc_random" in pivot.columns and "pr_auc_grouped" in pivot.columns:
        pivot["delta_pr_auc_grouped_minus_random"] = pivot["pr_auc_grouped"] - pivot["pr_auc_random"]
    if "roc_auc_random" in pivot.columns and "roc_auc_grouped" in pivot.columns:
        pivot["delta_roc_auc_grouped_minus_random"] = pivot["roc_auc_grouped"] - pivot["roc_auc_random"]
    if "brier_random" in pivot.columns and "brier_grouped" in pivot.columns:
        pivot["delta_brier_grouped_minus_random"] = pivot["brier_grouped"] - pivot["brier_random"]

    pivot.to_csv(OUT_DIR / "train_eval_comparison.csv", index=False)

    summary = {
        "n_rows_evaluated": int(common_n),
        "n_features": int(X.shape[1]),
        "n_endpoints": int(Y.shape[1]),
        "group_definition": "source||series_id||seq_length||ps_fraction_rounded",
        "files": {
            "report_csv": str(OUT_DIR / "train_eval_report.csv"),
            "comparison_csv": str(OUT_DIR / "train_eval_comparison.csv"),
        },
    }
    (OUT_DIR / "train_eval_summary.json").write_text(json.dumps(summary, indent=2))

    print("Saved:", OUT_DIR / "train_eval_report.csv")
    print("Saved:", OUT_DIR / "train_eval_comparison.csv")
    print("Saved:", OUT_DIR / "train_eval_summary.json")
    print()
    print("Top comparison rows:")
    print(pivot.head(15).to_string(index=False))


if __name__ == "__main__":
    main()
