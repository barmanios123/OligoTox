"""The app must not present an endpoint the evaluation doesn't support.

These tests tie artifacts/config/endpoint_evidence.json to the numbers in
reports/results/, so the shipped tiers can't drift from the evaluation, and
check that the gate actually holds in scoring and in the report.
"""
import json

import pandas as pd
import pytest

from oligotox.features.candidate import candidate_features
from oligotox.models.loader import load_endpoint_evidence
from oligotox.pipelines.score_candidate import score_candidate
from oligotox.reporting.oswg_report import OSWGSafetyReport

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
from derive_evidence_tiers import (  # noqa: E402
    SIMPLE_BASELINES,
    SOURCE_BASELINE,
    _model_label,
    classify,
)


@pytest.fixture
def evidence(root):
    return load_endpoint_evidence(root / "artifacts/config")


@pytest.fixture
def scored(root, tmp_path):
    return score_candidate(
        candidate_features=candidate_features("AUGCUAGCUAGCUAGCUAGC", "siRNA", "PS", ps_fraction=0.9),
        processed_dir=tmp_path / "no_training_data_needed",
        model_dir=root / "artifacts/models",
        config_dir=root / "artifacts/config",
    )


def test_shipped_tiers_match_the_baseline_results(root, evidence):
    """Re-derive every tier from the CSVs; the shipped config must agree."""
    baselines = pd.read_csv(root / "reports/results/baselines.csv")
    within = pd.read_csv(root / "reports/results/within_source.csv")

    for endpoint, sub in baselines.groupby("endpoint"):
        by_model = sub.set_index("model")
        model = by_model.loc[_model_label(sub)]
        rows = within[(within.endpoint == endpoint)
                      & (within.model == _model_label(sub))]
        tier, _, _ = classify(
            {
                "pr_auc": float(model["pr_auc"]),
                "ci_low": float(model["ci_low"]),
                "ci_high": float(model["ci_high"]),
            },
            float(by_model.loc[SOURCE_BASELINE, "pr_auc"]),
            max(float(by_model.loc[name, "pr_auc"]) for name in SIMPLE_BASELINES),
            float(model["base_rate"]),
            int(model["positive_sequences"]),
            bool((rows.pr_auc > rows.base_rate).any()) if len(rows) else False,
        )
        assert evidence["endpoints"][endpoint]["tier"] == tier, endpoint


def test_composite_endpoints_are_exactly_the_validated_ones(evidence):
    validated = sorted(e for e, v in evidence["endpoints"].items() if v["tier"] == "validated")
    assert sorted(evidence["composite_endpoints"]) == validated
    assert validated, "at least one endpoint must be presentable, or the app has nothing to show"


def test_a_source_confounded_endpoint_loses_to_the_source_baseline(evidence):
    """Guards the claim the app makes about this tier."""
    for endpoint, info in evidence["endpoints"].items():
        if info["tier"] == "source_confounded":
            assert info["source_only_pr_auc"] >= info["pr_auc"], endpoint


def test_every_scored_endpoint_carries_a_tier(scored, evidence):
    for endpoint, info in scored["per_endpoint"].items():
        assert info["evidence_tier"] in evidence["tiers"], endpoint
        assert info["evidence_note"], endpoint


def test_composite_only_averages_validated_endpoints(scored):
    validated = [e for e, i in scored["per_endpoint"].items() if i["evidence_tier"] == "validated"]
    assert scored["composite_basis"] == validated
    assert len(scored["composite_basis"]) < len(scored["per_endpoint"])
    assert scored["dominant_endpoint"] in validated


def test_unsupported_endpoints_are_not_called_positive(scored):
    for endpoint, info in scored["per_endpoint"].items():
        if info["evidence_tier"] != "validated":
            assert info["predicted_positive"] is None, endpoint


def test_report_raises_no_recommendation_from_an_unsupported_endpoint():
    """A certain-looking score on an unsupported endpoint must stay silent."""
    report = OSWGSafetyReport(
        compound_id="TEST-GATE",
        sequence="AUGC" * 5,
        assessment_date="2026-10-27",
        toxicity_predictions={"complement_activation": 1.0, "renal_accumulation": 1.0},
        pkpd_summary={},
        evidence_tiers={"complement_activation": "insufficient_data", "renal_accumulation": "weak"},
    )
    notes = report._generate_regulatory_notes()
    assert not any("Immunogenicity" in n or "Renal function" in n for n in notes)
    assert report._calculate_risk_score() == 0.0
    assert report._classify_safety() == "UNCLASSIFIED_NO_VALIDATED_ENDPOINT"

    payload = json.loads(report.generate_report(format="json"))
    assert payload["evidence"]["endpoints_supporting_assessment"] == []
    # The raw outputs stay visible; they just carry no conclusion.
    assert payload["safety_profile"]["toxicity_endpoints"]["complement_activation"] == 1.0


def test_report_still_acts_on_a_validated_endpoint():
    report = OSWGSafetyReport(
        compound_id="TEST-ACT",
        sequence="AUGC" * 5,
        assessment_date="2026-10-27",
        toxicity_predictions={"hepatocellular_injury": 0.9, "thrombocytopenia": 1.0},
        pkpd_summary={},
        evidence_tiers={"hepatocellular_injury": "validated", "thrombocytopenia": "insufficient_data"},
    )
    notes = report._generate_regulatory_notes()
    assert any("Liver function monitoring" in n for n in notes)
    # thrombocytopenia at 1.0 must not move the score.
    assert report._calculate_risk_score() == 0.9
    assert report._classify_safety() == "HIGH_RISK"
