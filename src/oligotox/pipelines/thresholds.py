# pipelines/thresholds.py

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable
import numpy as np
from sklearn.metrics import precision_score, recall_score, f1_score


@dataclass
class ThresholdResult:
    endpoint: str
    threshold: float
    precision: float
    recall: float
    f1: float
    support: int
    method: str


def _safe_metric(func, y_true, y_pred) -> float:
    try:
        return float(func(y_true, y_pred, zero_division=0))
    except TypeError:
        return float(func(y_true, y_pred))


def select_best_threshold(
    y_true: np.ndarray,
    y_prob: np.ndarray,
    endpoint: str,
    method: str = "f1",
    min_threshold: float = 0.20,
    max_threshold: float = 0.80,
    fallback_threshold: float = 0.50,
    n_grid: int = 61,
) -> ThresholdResult:
    y_true = np.asarray(y_true).astype(int)
    y_prob = np.asarray(y_prob).astype(float)

    positives = int(y_true.sum())
    if positives < 2:
        y_pred = (y_prob >= fallback_threshold).astype(int)
        return ThresholdResult(
            endpoint=endpoint,
            threshold=float(fallback_threshold),
            precision=_safe_metric(precision_score, y_true, y_pred),
            recall=_safe_metric(recall_score, y_true, y_pred),
            f1=_safe_metric(f1_score, y_true, y_pred),
            support=positives,
            method=f"{method}_fallback",
        )

    grid = np.linspace(min_threshold, max_threshold, n_grid)
    best = None

    for thr in grid:
        y_pred = (y_prob >= thr).astype(int)
        precision = _safe_metric(precision_score, y_true, y_pred)
        recall = _safe_metric(recall_score, y_true, y_pred)
        f1 = _safe_metric(f1_score, y_true, y_pred)

        if method == "recall_at_precision_050":
            score = recall if precision >= 0.50 else -1.0
        elif method == "precision_at_recall_050":
            score = precision if recall >= 0.50 else -1.0
        else:
            score = f1

        row = (score, f1, precision, recall, -abs(thr - 0.50), thr)
        if best is None or row > best[0]:
            best = (
                row,
                ThresholdResult(
                    endpoint=endpoint,
                    threshold=float(thr),
                    precision=float(precision),
                    recall=float(recall),
                    f1=float(f1),
                    support=positives,
                    method=method,
                ),
            )

    return best[1]


def batch_select_thresholds(
    y_true_by_endpoint: dict[str, np.ndarray],
    y_prob_by_endpoint: dict[str, np.ndarray],
    method: str = "f1",
    min_threshold: float = 0.20,
    max_threshold: float = 0.80,
    fallback_threshold: float = 0.50,
) -> dict[str, dict]:
    out = {}
    for endpoint, y_true in y_true_by_endpoint.items():
        if endpoint not in y_prob_by_endpoint:
            continue
        res = select_best_threshold(
            y_true=y_true,
            y_prob=y_prob_by_endpoint[endpoint],
            endpoint=endpoint,
            method=method,
            min_threshold=min_threshold,
            max_threshold=max_threshold,
            fallback_threshold=fallback_threshold,
        )
        out[endpoint] = {
            "threshold": res.threshold,
            "precision": res.precision,
            "recall": res.recall,
            "f1": res.f1,
            "support": res.support,
            "method": res.method,
        }
    return out