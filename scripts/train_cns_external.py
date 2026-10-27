"""Train the shipped cns_toxicity model on the two external CNS assays.

The previous cns_toxicity model was fit on 2 positive sequences and training
marked it `fallback_low_count_fit`. This replaces it with a model trained on
3,338 sequences from two independent assays, validated by
scripts/evaluate_cns_external.py by training on one assay and testing on the
other (ROC-AUC 0.83 and 0.78, no shared sequences).

The external data carries no chemistry annotation, so this model uses
sequence-derived features only. Its feature list is recorded separately from
the shared 55-column list in composite_config.json, and score_candidate reads
the per-endpoint list when one exists.

Writes artifacts/models/cns_toxicity.joblib,
artifacts/config/endpoint_features.json and the threshold into
endpoint_thresholds.json.
"""
from __future__ import annotations

import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.metrics import f1_score, roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict

from evaluate_cns_external import (  # same directory
    HAGEDORN,
    MOE,
    SEED,
    _feature_cols,
    _model,
    build,
)

ROOT = Path(__file__).resolve().parent.parent
MODELS = ROOT / "artifacts" / "models"
CONFIG = ROOT / "artifacts" / "config"
RESULTS = ROOT / "reports" / "results"
ENDPOINT = "cns_toxicity"


def main() -> None:
    moe, hag = build(MOE), build(HAGEDORN)
    shared = sorted(set(_feature_cols(moe)) & set(_feature_cols(hag)))

    moe["assay"], hag["assay"] = MOE, HAGEDORN
    full = pd.concat([moe, hag], ignore_index=True)
    X, y = full[shared].astype(float), full["toxic"].values

    # Out-of-fold probabilities on the pooled set, only to pick a threshold.
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
    oof = cross_val_predict(_model(), X, y, cv=cv, method="predict_proba")[:, 1]
    grid = np.linspace(0.05, 0.95, 91)
    threshold = float(max(grid, key=lambda t: f1_score(y, (oof >= t).astype(int), zero_division=0)))

    model = _model().fit(X, y)
    MODELS.mkdir(parents=True, exist_ok=True)
    joblib.dump(model, MODELS / f"{ENDPOINT}.joblib")

    # Per-endpoint feature columns, since this model does not use the shared 55.
    features_path = CONFIG / "endpoint_features.json"
    features = json.loads(features_path.read_text()) if features_path.exists() else {}
    features[ENDPOINT] = {
        "feature_columns": shared,
        "source": "external CNS assays; sequence-derived features only",
        "note": "The training data carries no chemistry annotation, so this model "
                "uses no chemistry or modality feature.",
    }
    features_path.write_text(json.dumps(features, indent=2) + "\n")

    thresholds_path = CONFIG / "endpoint_thresholds.json"
    thresholds = json.loads(thresholds_path.read_text())
    thresholds[ENDPOINT] = {"threshold": threshold, "source": "external CNS assays, F1 on pooled OOF"}
    thresholds_path.write_text(json.dumps(thresholds, indent=2) + "\n")

    summary = {
        "endpoint": ENDPOINT,
        "sequences": int(len(full)),
        "positives": int(y.sum()),
        "base_rate": float(y.mean()),
        "assays": {MOE: int(len(moe)), HAGEDORN: int(len(hag))},
        "n_features": len(shared),
        "feature_columns": shared,
        "pooled_oof_roc_auc": float(roc_auc_score(y, oof)),
        "threshold": threshold,
        "validation": "cross-assay, see reports/results/cns_external.json",
    }
    (RESULTS / "cns_training_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    print(f"trained {ENDPOINT} on {len(full)} sequences ({y.sum()} positive) "
          f"from {len(summary['assays'])} assays, {len(shared)} features")
    print(f"pooled out-of-fold ROC-AUC {summary['pooled_oof_roc_auc']:.3f}, threshold {threshold:.2f}")


if __name__ == "__main__":
    main()
