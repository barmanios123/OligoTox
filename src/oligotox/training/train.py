from __future__ import annotations

import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from oligotox.features.v03_bridge import build_v03_feature_df

from sklearn.calibration import CalibratedClassifierCV
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.dummy import DummyClassifier
from sklearn.metrics import average_precision_score, brier_score_loss
from sklearn.model_selection import StratifiedKFold, GroupKFold, cross_val_predict

from oligotox.config.paths import (
    PROCESSED_DATA_DIR,
    MODEL_ARTIFACTS_DIR,
    METRICS_ARTIFACTS_DIR,
    CONFIG_ARTIFACTS_DIR,
)
from oligotox.config.risk_config import (
    ENDPOINT_SEVERITY_WEIGHTS,
    DEFAULT_THRESHOLD_POLICY,
    RISK_BANDS,
)
from oligotox.pipelines.build_training_dataset import run as build_training_dataset
from oligotox.pipelines.thresholds import batch_select_thresholds

from src.oligotox.integrations_v03.datasail_splits import (
    scaffold_sequence_split,
    datasail_split,
    check_leakage,
    evaluate_endpoint,
)


def _json_dump(path: Path, obj: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2))


def _make_model() -> CalibratedClassifierCV:
    base = LogisticRegression(
        max_iter=2000,
        class_weight="balanced",
        solver="liblinear",
    )
    return CalibratedClassifierCV(base, method="sigmoid", cv=2)


def _prepare_training_matrix(X: pd.DataFrame) -> pd.DataFrame:
    X = X.copy()

    for col in X.columns:
        X[col] = pd.to_numeric(X[col], errors="coerce")

    imputer = SimpleImputer(
        strategy="constant",
        fill_value=0.0,
        keep_empty_features=True,
    )
    arr = imputer.fit_transform(X)
    return pd.DataFrame(arr, columns=X.columns, index=X.index)


def _make_constant_model() -> DummyClassifier:
    return DummyClassifier(strategy="prior")


def _positive_class_proba(model, X):
    proba = model.predict_proba(X)
    classes = list(getattr(model, "classes_", []))
    if proba.ndim != 2:
        raise ValueError("predict_proba returned unexpected shape")
    if len(classes) == 2:
        if 1 in classes:
            pos_idx = classes.index(1)
            return proba[:, pos_idx]
        return proba[:, -1]
    if len(classes) == 1:
        only = classes[0]
        if only == 1:
            return proba[:, 0]
        return proba[:, 0] * 0.0
    return proba[:, -1]


def _get_groups_for_cv(X: pd.DataFrame) -> "np.ndarray | None":
    """Return grouping labels for leakage-aware CV if available.

    Priority:
      1) gapmer_id series-level grouping
      2) oligo_id grouping
      3) None (fall back to standard StratifiedKFold)
    """
    import numpy as np

    for col in ["gapmer_id", "oligo_id"]:
        if col in X.columns:
            g = X[col]
            if g.nunique(dropna=True) > 1:
                return g.astype(str).fillna("NA").values
    return None


def _run_datasail_diagnostic_split(
    X: pd.DataFrame,
    Y: pd.DataFrame,
    meta: dict,
) -> dict:
    """
    Run a leakage-aware split diagnostic (DataSAIL → scaffold fallback).

    This does NOT affect training; it only writes a leakage report to metrics.
    """
    # Choose a sequence column if available.
    seq_col_candidates = ["sequence", "antisense_sequence", "oligo_sequence"]
    feature_columns = meta.get("feature_columns", list(X.columns))
    sequence_col = None
    for c in seq_col_candidates:
        if c in feature_columns or c in X.columns:
            sequence_col = c
            break
    if sequence_col is None:
        return {"enabled": False, "reason": "no_sequence_column"}

    # Ensure we have the raw sequence values, not imputed numerics.
    if sequence_col not in X.columns:
        try:
            raw_X = pd.read_csv(PROCESSED_DATA_DIR / "feature_matrix.csv")
        except Exception:
            return {"enabled": False, "reason": "sequence_column_missing_and_raw_load_failed"}
        if sequence_col not in raw_X.columns:
            return {"enabled": False, "reason": f"sequence_column_{sequence_col}_not_found_in_raw"}
        seqs = raw_X[sequence_col].fillna("").astype(str).tolist()
    else:
        seqs = X[sequence_col].fillna("").astype(str).tolist()

    # Simple binary label proxy: union of any positive endpoint.
    y_any = (Y.fillna(0) > 0).any(axis=1).astype(int).values

    # Try DataSAIL; fall back to scaffold split.
    try:
        train_idx, val_idx, test_idx = datasail_split(
            sequences=seqs,
            y=y_any,
            test_size=0.2,
            val_size=0.1,
            similarity_threshold=0.6,
            technique="I1",
        )
        strategy = "datasail_I1"
    except Exception:
        train_idx, val_idx, test_idx = scaffold_sequence_split(
            sequences=seqs,
            y=y_any,
            test_size=0.2,
            val_size=0.1,
        )
        strategy = "scaffold"

    train_seqs = [seqs[i] for i in train_idx]
    test_seqs = [seqs[i] for i in test_idx]
    leakage_report = check_leakage(train_seqs, test_seqs)

    return {
        "enabled": True,
        "strategy": strategy,
        "sequence_column": sequence_col,
        "train_size": int(len(train_idx)),
        "val_size": int(len(val_idx)),
        "test_size": int(len(test_idx)),
        "leakage_report": leakage_report,
    }


def main() -> None:
    print("=" * 80)
    print("OligoTox training run")
    print("=" * 80)

    build_info = build_training_dataset(
        str(PROCESSED_DATA_DIR / "training_dataset.csv"),
        str(PROCESSED_DATA_DIR / "feature_matrix.csv"),
        str(PROCESSED_DATA_DIR / "label_matrix.csv"),
        str(PROCESSED_DATA_DIR / "dataset_metadata.json"),
    )

    X = pd.read_csv(PROCESSED_DATA_DIR / "feature_matrix.csv")
    Y = pd.read_csv(PROCESSED_DATA_DIR / "label_matrix.csv")
    meta = json.loads((PROCESSED_DATA_DIR / "dataset_metadata.json").read_text())

    X = _prepare_training_matrix(X)

    print(f"Rows: {len(X)}")
    print(f"Feature columns: {X.shape[1]}")
    print(f"Endpoint labels: {Y.shape[1]}")
    print()

    summary: dict[str, dict] = {}
    oof_true: dict[str, np.ndarray] = {}
    oof_prob: dict[str, np.ndarray] = {}

    v03_cols_existing = [
        c for c in X.columns if c.startswith(("chem_", "d2_", "struct_", "ctx_"))
    ]
    if not v03_cols_existing:
        v03 = build_v03_feature_df(X, use_viennarna=False, use_embeddings=False)
        overlap = [c for c in v03.columns if c in X.columns]
        if overlap:
            X = X.drop(columns=overlap)
        X = pd.concat([X.reset_index(drop=True), v03.reset_index(drop=True)], axis=1)
        print(f"Re-attached v0.3 mechanistic features: {v03.shape[1]} columns")
    else:
        print(f"v0.3 mechanistic features already present: {len(v03_cols_existing)} columns")

    print("=" * 80)
    print("Training summary")
    print("=" * 80)

    for endpoint in Y.columns:
        y = pd.to_numeric(Y[endpoint], errors="coerce").fillna(0).astype(int).values
        positives = int(y.sum())
        negatives = int((y == 0).sum())
        positive_rate = float(np.mean(y))
        samples = int(len(y))
        min_class = min(positives, negatives)

        if positives == 0 or negatives == 0:
            status = "fallback_single_class_fit"
            model = _make_constant_model()
            model.fit(X, y)
            y_prob = _positive_class_proba(model, X)
            pr_auc = None
            brier = None
            n_splits = None
            endpoint_metrics = None
        elif min_class < 3:
            status = "fallback_low_count_fit"
            model = _make_model()
            model.fit(X, y)
            y_prob = _positive_class_proba(model, X)
            pr_auc = None
            brier = None
            n_splits = None
            endpoint_metrics = None
        else:
            status = "ok"
            n_splits = min(3, min_class)
            groups = _get_groups_for_cv(X)
            if groups is not None:
                unique_groups = len(set(map(str, groups)))
                if unique_groups >= n_splits:
                    cv = GroupKFold(n_splits=n_splits)
                    cv_iter = cv.split(X, y, groups=groups)
                else:
                    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)
                    cv_iter = cv.split(X, y)
            else:
                cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)
                cv_iter = cv.split(X, y)

            model = _make_model()

            y_prob = cross_val_predict(
                model,
                X,
                y,
                cv=cv_iter,
                method="predict_proba",
            )[:, 1]

            pr_auc = float(average_precision_score(y, y_prob))
            brier = float(brier_score_loss(y, y_prob))
            try:
                endpoint_metrics = evaluate_endpoint(y, y_prob, endpoint_name=endpoint)
            except Exception:
                endpoint_metrics = None

        model.fit(X, y)
        joblib.dump(model, MODEL_ARTIFACTS_DIR / f"{endpoint}.joblib")

        oof_true[endpoint] = y
        oof_prob[endpoint] = y_prob

        summary[endpoint] = {
            "samples": samples,
            "positive_rate": positive_rate,
            "status": status,
            "pr_auc": pr_auc,
            "brier": brier,
            "positives": positives,
            "negatives": negatives,
            "n_splits": n_splits,
            "datasail_eval": endpoint_metrics,
        }

        pr_auc_str = f"{pr_auc:0.3f}" if pr_auc is not None else "NA"
        brier_str = f"{brier:0.3f}" if brier is not None else "NA"
        print(
            f"{endpoint:28s} "
            f"samples={samples:<5d} "
            f"positive_rate={positive_rate:0.3f} "
            f"status={status:22s} "
            f"pr_auc={pr_auc_str:<6s} "
            f"brier={brier_str}"
        )

    thresholds = batch_select_thresholds(
        y_true_by_endpoint=oof_true,
        y_prob_by_endpoint=oof_prob,
        method=DEFAULT_THRESHOLD_POLICY["method"],
        min_threshold=DEFAULT_THRESHOLD_POLICY["min_threshold"],
        max_threshold=DEFAULT_THRESHOLD_POLICY["max_threshold"],
        fallback_threshold=DEFAULT_THRESHOLD_POLICY["fallback_threshold"],
    )

    composite_cfg = {
        "endpoint_weights": ENDPOINT_SEVERITY_WEIGHTS,
        "risk_bands": RISK_BANDS,
        "trained_feature_columns": list(X.columns),
    }

    _json_dump(METRICS_ARTIFACTS_DIR / "training_summary.json", summary)
    _json_dump(CONFIG_ARTIFACTS_DIR / "endpoint_thresholds.json", thresholds)
    _json_dump(CONFIG_ARTIFACTS_DIR / "composite_config.json", composite_cfg)
    _json_dump(METRICS_ARTIFACTS_DIR / "build_info.json", build_info)

    feature_matrix_path = PROCESSED_DATA_DIR / "feature_matrix.csv"
    label_matrix_path = PROCESSED_DATA_DIR / "label_matrix.csv"

    X.to_csv(feature_matrix_path, index=False)
    Y.to_csv(label_matrix_path, index=False)

    feature_matrix_final = pd.read_csv(feature_matrix_path, low_memory=False)
    label_matrix_final = pd.read_csv(label_matrix_path, low_memory=False)

    meta["rows"] = int(feature_matrix_final.shape[0])
    meta["n_features"] = int(feature_matrix_final.shape[1])
    meta["n_labels"] = int(label_matrix_final.shape[1])
    meta["feature_columns"] = list(feature_matrix_final.columns)
    meta["trained_feature_columns"] = list(feature_matrix_final.columns)
    meta["label_columns"] = list(label_matrix_final.columns)

    _json_dump(PROCESSED_DATA_DIR / "dataset_metadata.json", meta)

    # DataSAIL/scaffold leakage diagnostic
    datasail_info = _run_datasail_diagnostic_split(feature_matrix_final, label_matrix_final, meta)
    _json_dump(METRICS_ARTIFACTS_DIR / "datasail_split_diagnostic.json", datasail_info)

    print()
    print(f"Saved models to: {MODEL_ARTIFACTS_DIR}")
    print(f"Saved metrics to: {METRICS_ARTIFACTS_DIR}")
    print(f"Saved config to: {CONFIG_ARTIFACTS_DIR}")
    print(f"Saved feature matrix to: {feature_matrix_path}")
    print(f"Saved label matrix to: {label_matrix_path}")
    print(f"Saved dataset metadata to: {PROCESSED_DATA_DIR / 'dataset_metadata.json'}")


if __name__ == "__main__":
    main()