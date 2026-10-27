import json
import math

from oligotox.pkpd.bridge import compute_pkpd_enriched_features
from oligotox.reporting.oswg_report import OSWGSafetyReport


def test_pkpd_outputs_are_finite_and_positive():
    out = compute_pkpd_enriched_features("AUGCUAGCUAGCUAGCUAGC", "siRNA", "PS", 3.0, "weekly")
    for key in ("liver_auc_pmol_g_h", "liver_cmax_pmol_g", "liver_half_life_h"):
        assert math.isfinite(out[key]) and out[key] > 0
    assert 0 <= out["exposure_score"] <= 4.99


def test_pkpd_exposure_rises_with_dose():
    low = compute_pkpd_enriched_features("AUGC" * 5, "siRNA", "PS", 1.0, "weekly")
    high = compute_pkpd_enriched_features("AUGC" * 5, "siRNA", "PS", 10.0, "weekly")
    assert high["liver_auc_pmol_g_h"] > low["liver_auc_pmol_g_h"]


def test_report_generates_json():
    report = OSWGSafetyReport(
        compound_id="TEST-1",
        sequence="AUGCUAGCUAGCUAGCUAGC",
        assessment_date="2026-10-27",
        toxicity_predictions={"hepatocellular_injury": 0.3, "thrombocytopenia": 0.05},
        pkpd_summary={"exposure_score": 1.2},
    )
    payload = json.loads(report.generate_report(format="json"))
    assert payload
    assert 0.0 <= report._calculate_risk_score() <= 1.0
    assert isinstance(report._classify_safety(), str)


def test_report_omits_efficacy_when_not_scored():
    report = OSWGSafetyReport(
        compound_id="TEST-2",
        sequence="AUGC" * 5,
        assessment_date="2026-10-27",
        toxicity_predictions={},
        pkpd_summary={},
    )
    assert "efficacy" not in json.loads(report.generate_report(format="json"))
