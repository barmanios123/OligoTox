# pipelines/composite.py

from __future__ import annotations

from typing import Dict, Any
import numpy as np


def probability_to_band(p: float, bands: dict[str, float]) -> str:
    if p >= bands["high"]:
        return "high"
    if p >= bands["medium"]:
        return "medium"
    if p >= bands["low"]:
        return "elevated"
    return "low"


def weighted_composite_score(
    endpoint_probs: Dict[str, float],
    endpoint_weights: Dict[str, float],
) -> float:
    pairs = [(endpoint_probs[e], endpoint_weights[e]) for e in endpoint_probs if e in endpoint_weights]
    if not pairs:
        return 0.0
    probs = np.array([p for p, _ in pairs], dtype=float)
    weights = np.array([w for _, w in pairs], dtype=float)
    if weights.sum() <= 0:
        return float(np.mean(probs))
    return float(np.average(probs, weights=weights))


def top_risk_endpoint(endpoint_probs: Dict[str, float]) -> str | None:
    if not endpoint_probs:
        return None
    return max(endpoint_probs, key=endpoint_probs.get)


def summarize_endpoint(
    endpoint: str,
    prob: float,
    threshold: float | None,
    bands: dict[str, float],
) -> Dict[str, Any]:
    called_positive = bool(prob >= threshold) if threshold is not None else False
    return {
        "prob": float(prob),
        "threshold": float(threshold) if threshold is not None else None,
        "predicted_positive": called_positive,
        "risk_band": probability_to_band(float(prob), bands),
    }