from __future__ import annotations

import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


ROOT = Path(".")
PROCESSED = ROOT / "data" / "processed"
ARTIFACTS = ROOT / "artifacts"
MODEL_DROP = ARTIFACTS / "models" / "latest_train_baseline"
MODEL_DROP.mkdir(parents=True, exist_ok=True)

X_CANDIDATES = [
    PROCESSED / "X_full.csv",
    PROCESSED / "features.csv",
    PROCESSED / "X.csv",
    PROCESSED / "feature_matrix.csv",
]
Y_CANDIDATES = [
    PROCESSED / "y_full.csv",
    PROCESSED / "labels.csv",
    PROCESSED / "y.csv",
    PROCESSED / "label_matrix.csv",
]

def find_first(paths):
    for p in paths:
        if p.exists():
            return p
    raise FileNotFoundError(f"None of these files exist: {paths}")

x_path = find_first(X_CANDIDATES)
y_path = find_first(Y_CANDIDATES)

X = pd.read_csv(x_path, low_memory=False)
y = pd.read_csv(y_path, low_memory=False)

non_feature_cols = [c for c in ["candidate_id", "sequence", "oligo_id"] if c in X.columns]
X_num = X.drop(columns=non_feature_cols, errors="ignore").copy()

for col in X_num.columns:
    X_num[col] = pd.to_numeric(X_num[col], errors="coerce")

all_nan_cols = [c for c in X_num.columns if X_num[c].isna().all()]
constant_cols = [c for c in X_num.columns if X_num[c].nunique(dropna=True) <= 1]
drop_cols = sorted(set(all_nan_cols + constant_cols))
X_model = X_num.drop(columns=drop_cols, errors="ignore")

usable_labels = []
label_summary = []

for col in y.columns:
    yc = pd.to_numeric(y[col], errors="coerce")
    mask = yc.notna()
    yv = yc[mask]
    if len(yv) == 0:
        label_summary.append({"endpoint": col, "status": "all_nan"})
        continue
    if yv.nunique() < 2:
        label_summary.append({"endpoint": col, "status": "single_class", "n_rows": int(mask.sum())})
        continue
    usable_labels.append(col)

cv = StratifiedKFold(n_splits=3, shuffle=True, random_state=42)
manifest = {
    "x_path": str(x_path),
    "y_path": str(y_path),
    "n_rows_X": int(len(X_model)),
    "n_features_retained": int(X_model.shape[1]),
    "n_features_dropped_all_nan": int(len(all_nan_cols)),
    "n_features_dropped_constant": int(len(constant_cols)),
    "usable_labels": [],
    "models": {},
    "summary": [],
}

summary_rows = []

for endpoint in usable_labels:
    yc = pd.to_numeric(y[endpoint], errors="coerce")
    mask = yc.notna()

    X_ep = X_model.loc[mask].copy()
    y_ep = yc.loc[mask].astype(int).copy()

    if y_ep.nunique() < 2 or len(y_ep) < 8:
        summary_rows.append({
            "endpoint": endpoint,
            "status": "skipped_insufficient_rows",
            "n_rows": int(len(y_ep)),
        })
        continue

    pipe = Pipeline(
        steps=[
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
            ("model", LogisticRegression(max_iter=2000, class_weight="balanced")),
        ]
    )

    try:
        proba = cross_val_predict(pipe, X_ep, y_ep, cv=cv, method="predict_proba")[:, 1]
        roc_auc = float(roc_auc_score(y_ep, proba))
        pr_auc = float(average_precision_score(y_ep, proba))
    except Exception:
        roc_auc = None
        pr_auc = None

    pipe.fit(X_ep, y_ep)

    model_path = MODEL_DROP / f"{endpoint}.joblib"
    joblib.dump(pipe, model_path)

    manifest["usable_labels"].append(endpoint)
    manifest["models"][endpoint] = str(model_path)
    summary_rows.append({
        "endpoint": endpoint,
        "status": "trained",
        "n_rows": int(len(y_ep)),
        "positive_rate": float(y_ep.mean()),
        "roc_auc_cv": roc_auc,
        "pr_auc_cv": pr_auc,
        "model_path": str(model_path),
    })

summary_df = pd.DataFrame(summary_rows).sort_values(["status", "endpoint"]).reset_index(drop=True)
summary_path = MODEL_DROP / "train_eval_report.csv"
summary_df.to_csv(summary_path, index=False)

manifest["summary"] = summary_rows
manifest_path = MODEL_DROP / "model_manifest.json"
manifest_path.write_text(json.dumps(manifest, indent=2))

meta = {
    "model_dir": str(MODEL_DROP),
    "summary_report": str(summary_path),
    "manifest": str(manifest_path),
    "n_models_saved": int(sum(1 for r in summary_rows if r["status"] == "trained")),
}
(MODEL_DROP / "train_metadata.json").write_text(json.dumps(meta, indent=2))

print(json.dumps(meta, indent=2))
