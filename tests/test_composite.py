from oligotox.pipelines.composite import (
    probability_to_band,
    summarize_endpoint,
    top_risk_endpoint,
    weighted_composite_score,
)

BANDS = {"low": 0.2, "medium": 0.4, "high": 0.7}


def test_probability_to_band_edges():
    assert probability_to_band(0.1, BANDS) == "low"
    assert probability_to_band(0.2, BANDS) == "elevated"
    assert probability_to_band(0.4, BANDS) == "medium"
    assert probability_to_band(0.7, BANDS) == "high"


def test_weighted_composite_ignores_unweighted_endpoints():
    probs = {"a": 1.0, "b": 0.0, "c": 0.9}
    assert weighted_composite_score(probs, {"a": 3.0, "b": 1.0}) == 0.75


def test_weighted_composite_empty_is_zero():
    assert weighted_composite_score({}, {"a": 1.0}) == 0.0


def test_top_risk_endpoint():
    assert top_risk_endpoint({"a": 0.1, "b": 0.6}) == "b"
    assert top_risk_endpoint({}) is None


def test_summarize_endpoint_without_threshold_is_not_positive():
    out = summarize_endpoint("a", 0.9, None, BANDS)
    assert out["predicted_positive"] is False
    assert out["risk_band"] == "high"
