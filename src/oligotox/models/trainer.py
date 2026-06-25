from pathlib import Path
import numpy as np
import pandas as pd

from sklearn.calibration import CalibratedClassifierCV
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.pipeline import Pipeline

from oligotox.models.persistence import save_model, save_json

LABEL_COLUMNS = [
    "hepatocellular_injury",
    "hepatic_vacuolation",
    "renal_tubular_injury",
    "renal_accumulation",
    "thrombocytopenia",
    "coagulation_effect",
    "complement_activation",
    "innate_immune_activation",
    "cytokine_release",
    "kupffer_cell_activation",
    "cns_toxicity",
    "offtarget_gene_silencing",
    "injection_site_reaction",
]

RANDOM_STATE = 42
DEFAULT_N_SPLITS = 3
MIN_POSITIVES_FOR_TRAIN = 2
C_GRID = [0.01, 0.1, 1.0, 10.0]


def build_base_estimator():
    return Pipeline(
        steps=[
            ("imputer", SimpleImputer(strategy="median")),
            (
                "clf",
                LogisticRegression(
                    class_weight="balanced",
                    solver="liblinear",
                    penalty="l2",
                    max_iter=2000,
                    random_state=RANDOM_STATE,
                ),
            ),
        ]
    )


def choose_cv(y):
    positives = int(np.sum(y == 1))
    negatives = int(np.sum(y == 0))
    n_splits = min(DEFAULT_N_SPLITS, positives, negatives)
    if n_splits < 2:
        return None
    return StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=RANDOM_STATE)


def fit_one_endpoint(X, y, endpoint_name):
    y = np.asarray(y).astype(int)
    positives = int(np.sum(y == 1))
    negatives = int(np.sum(y == 0))

    metrics = {
        "endpoint": endpoint_name,
        "samples": int(len(y)),
        "positives": positives,
        "negatives": negatives,
        "positive_rate": float(np.mean(y)),
    }

    if positives < MIN_POSITIVES_FOR_TRAIN or negatives < MIN_POSITIVES_FOR_TRAIN:
        metrics["status"] = "skipped_too_few_examples"
        return None, metrics

    cv = choose_cv(y)
    base = build_base_estimator()

    if cv is None:
        base.fit(X, y)
        prob = base.predict_proba(X)[:, 1]
        metrics.update(
            {
                "status": "fit_no_cv",
                "best_C": None,
                "roc_auc_train": float(roc_auc_score(y, prob)) if len(np.unique(y)) > 1 else None,
                "pr_auc_train": float(average_precision_score(y, prob)),
                "brier_train": float(brier_score_loss(y, prob)),
            }
        )
        return base, metrics

    search = GridSearchCV(
        estimator=base,
        param_grid={"clf__C": C_GRID},
        scoring="average_precision",
        cv=cv,
        refit=True,
    )
    search.fit(X, y)

    calibrated = CalibratedClassifierCV(
        estimator=search.best_estimator_,
        method="sigmoid",
        cv=cv,
    )
    calibrated.fit(X, y)
    prob = calibrated.predict_proba(X)[:, 1]

    metrics.update(
        {
            "status": "ok",
            "best_C": float(search.best_params_["clf__C"]),
            "roc_auc_train": float(roc_auc_score(y, prob)) if len(np.unique(y)) > 1 else None,
            "pr_auc_train": float(average_precision_score(y, prob)),
            "brier_train": float(brier_score_loss(y, prob)),
        }
    )
    return calibrated, metrics


def train_all(feature_csv, label_csv, model_dir, metrics_dir):
    X = pd.read_csv(feature_csv)
    Y = pd.read_csv(label_csv)

    model_dir = Path(model_dir)
    metrics_dir = Path(metrics_dir)
    model_dir.mkdir(parents=True, exist_ok=True)
    metrics_dir.mkdir(parents=True, exist_ok=True)

    label_cols = [c for c in LABEL_COLUMNS if c in Y.columns]
    if not label_cols:
        raise ValueError("No known endpoint labels found in label matrix.")

    summary = []
    for endpoint in label_cols:
        model, metrics = fit_one_endpoint(X, Y[endpoint], endpoint)
        summary.append(metrics)
        save_json(metrics_dir / f"{endpoint}.json", metrics)
        if model is not None:
            save_model(model_dir / f"{endpoint}.joblib", model)

    save_json(metrics_dir / "training_summary.json", summary)
    return summary