"""The shipped estimator must be scaled and calibrated on enough folds.

`CalibratedClassifierCV(base, cv=2)` on unscaled features was costing
hepatocellular injury 0.25 PR-AUC and all of its calibration: each
calibrated sub-model saw half the rows, and `liblinear` applied one penalty
across coefficients whose natural scales run from 0.5 to 20.

Neither defect is visible in any output. The model trains, predicts, and
reports a number; the number is just much worse than it needs to be. These
assert the shape of the estimator, because that is where the bug lived.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest
import numpy as np
from sklearn.base import clone
from sklearn.calibration import CalibratedClassifierCV
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from oligotox.models.estimator import AdaptiveCalibrated
from oligotox.training.train import _make_model

ROOT = Path(__file__).resolve().parent.parent


def _fit(positives: int, rows: int = 60):
    X = np.random.RandomState(0).rand(rows, 3)
    y = np.array([1] * positives + [0] * (rows - positives))
    model = clone(_make_model())
    model.fit(X, y)
    return model


def _calibrator(model):
    inner = model.pipeline_
    assert model.folds_, "this model was fitted without calibration"
    assert isinstance(inner, Pipeline), "the estimator must scale its features"
    assert isinstance(inner.steps[0][1], StandardScaler)
    calibrator = inner.steps[-1][1]
    assert isinstance(calibrator, CalibratedClassifierCV)
    return calibrator


def test_the_features_are_scaled_before_a_penalised_fit():
    _calibrator(_fit(20))


def test_calibration_does_not_train_on_half_the_data():
    """cv=2 is the defect; five folds leave 80% in each fit."""
    assert _fit(20).folds_ == 5


def test_a_sparse_endpoint_gets_as_many_folds_as_it_can_support():
    """An endpoint with 3 positives cannot be calibrated on 5 folds."""
    assert _fit(4).folds_ == 4
    assert _fit(3).folds_ == 3
    assert _fit(2).folds_ == 2, "never fewer than two"
    assert _fit(1).folds_ == 0, "a class of one cannot be split, so no calibration"
    assert _fit(50, rows=200).folds_ == 5, "never more than five"


def test_the_fold_count_comes_from_the_rows_the_fit_was_given():
    """Deciding it from the endpoint's totals crashes inside an outer fold.

    An endpoint with 4 positives reaches a 5-fold outer split with 3, and a
    calibrator configured for 4 folds raises rather than adapting.
    """
    assert _fit(3).folds_ == 3


def test_the_fix_is_visible_in_the_endpoint_that_has_the_data():
    """Hepatocellular injury is the endpoint with enough rows to show it.

    Guards the direction rather than an exact value: the point is that the
    model is no longer worse than quoting the base rate.
    """
    summary = json.loads(
        (ROOT / "reports/results/training_summary_after.json").read_text())
    spec = summary.get("hepatocellular_injury", {})
    if spec.get("pr_auc") is None:
        pytest.skip("endpoint not trained in this build")
    assert spec["pr_auc"] > 0.5, (
        f"PR-AUC {spec['pr_auc']:.3f}; it was 0.371 with the cv=2 unscaled "
        "estimator and about 0.62 once scaled and calibrated on 5 folds"
    )


def test_the_model_can_be_loaded_outside_the_trainer():
    """The trainer runs as `python -m`, so its module is `__main__`.

    A class defined beside the training loop pickles as
    `__main__.AdaptiveCalibrated` and raises AttributeError in the app, the
    scorer and the external check. The shipped models were briefly exactly
    that broken. This asserts the import path the pickle will carry.
    """
    model = _make_model()
    assert type(model) is AdaptiveCalibrated
    assert type(model).__module__ == "oligotox.models.estimator", (
        "the estimator must come from an importable module, not from "
        f"{type(model).__module__}"
    )


def test_every_shipped_model_actually_opens():
    """Loading is the only proof that a pickle is loadable."""
    import joblib

    models = sorted((ROOT / "artifacts" / "models").glob("*.joblib"))
    assert models, "no models in this build"
    for path in models:
        joblib.load(path)


def test_the_model_reports_the_features_it_was_fitted_on():
    """The How it works page lists inputs by reading them off the model.

    A hand-written estimator never calls sklearn's `_validate_data`, so
    `feature_names_in_` is not set for free, and the page would show every
    toxicity model consuming nothing at all.
    """
    import pandas as pd

    frame = pd.DataFrame(np.random.RandomState(0).rand(40, 3),
                         columns=["gc_content", "seq_length", "tlr7_score"])
    y = np.array([1] * 10 + [0] * 30)
    model = clone(_make_model())
    model.fit(frame, y)
    assert list(model.feature_names_in_) == list(frame.columns)
