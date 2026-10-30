"""The shipped toxicity estimator.

It lives in its own module rather than beside the training loop because the
trainer runs as `python -m oligotox.training.train`, which makes its module
`__main__`. A class defined there pickles as `__main__._AdaptiveCalibrated`
and cannot be loaded by anything else, so every model this project ships
would have failed to open in the app, in the scorer, and in the external
check. Defining it here gives the pickle a stable import path.
"""
from __future__ import annotations

import numpy as np
from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.calibration import CalibratedClassifierCV
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


class AdaptiveCalibrated(BaseEstimator, ClassifierMixin):
    """Logistic regression, scaled, calibrated on as many folds as fit.

    Two defects lived in the three lines this replaces, and together they
    cost hepatocellular injury 0.25 PR-AUC and all of its calibration:

      - `cv=2` meant every calibrated sub-model was fitted on half the rows,
        and the endpoint with the most data is the one that suffers most
        from being given half of it. Raising it to 5 is worth +0.18 PR-AUC
        on that endpoint; paired bootstrap over sequences puts the whole fix
        at +0.195 [+0.024, +0.342], P(better) 0.994.
      - The features were never scaled, while `liblinear` applies one
        penalty to coefficients whose natural scales run from 0.5
        (gc_content) to 20 (seq_length), so regularisation fell unevenly
        across them. Worth a further +0.07.

    Brier skill on hepatocellular injury goes from -0.054, which is worse
    than quoting the base rate, to +0.128. A gradient boosting model scores
    higher still (0.679 against 0.619) but not distinguishably so
    (P(better) 0.90), and it would give up the calibrated probabilities the
    rest of the platform is built on, so the estimator was fixed rather than
    replaced. `scripts/evaluate_hepatotox_model.py` has the comparison.

    The fold count is chosen inside `fit`, from the labels actually handed
    to it. Deciding it from the endpoint's total class sizes is wrong,
    because each outer CV fold holds out part of the smaller class: an
    endpoint with 4 positives arrives here with 3.
    """

    def __init__(self, max_folds: int = 5):
        self.max_folds = max_folds

    def fit(self, X, y):
        # sklearn sets this inside `_validate_data`, which a hand-written
        # estimator never calls. The app's How it works page lists each
        # model's real inputs by reading it off the model, so without this
        # every toxicity model would appear to consume nothing.
        columns = getattr(X, "columns", None)
        if columns is not None:
            self.feature_names_in_ = np.asarray(columns, dtype=object)
            self.n_features_in_ = len(columns)
        y = np.asarray(y)
        counts = np.bincount(y.astype(int)) if len(y) else np.array([0])
        smallest = int(counts[counts > 0].min()) if counts.any() else 0
        base = LogisticRegression(
            max_iter=2000, class_weight="balanced", solver="liblinear",
        )
        # A class of one cannot be split at all, so there is nothing to
        # calibrate on. The trainer already marks those endpoints
        # `fallback_low_count_fit`; this keeps the fit from raising.
        folds = 0 if smallest < 2 else max(2, min(self.max_folds, smallest))
        self.folds_ = folds
        self.pipeline_ = make_pipeline(
            StandardScaler(),
            base if folds == 0
            else CalibratedClassifierCV(base, method="sigmoid", cv=folds),
        )
        self.pipeline_.fit(X, y)
        self.classes_ = self.pipeline_.classes_
        return self

    def predict_proba(self, X):
        return self.pipeline_.predict_proba(X)

    def predict(self, X):
        return self.pipeline_.predict(X)
