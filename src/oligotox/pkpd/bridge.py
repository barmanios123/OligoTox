"""
PK/PD Bridge to OligoTox Toxicity Predictions
Phase 2B: Calibrated hepatotox adjustment logic
"""
from typing import Dict, Any
from .simulator import GalNAcPKPDSimulator


def _hepatotox_pkpd_adjustment(exposure_score: float) -> float:
    """
    Hepatotoxicity adjustment based on exposure.
    High exposure ALWAYS increases risk relative to baseline.
    """
    if exposure_score <= 1.0:
        return 0.00  # No adjustment at baseline
    if exposure_score <= 2.0:
        return 0.02
    if exposure_score <= 3.0:
        return 0.05
    if exposure_score <= 4.0:
        return 0.10
    return 0.18  # High exposure significantly increases risk


def compute_pkpd_enriched_features(
    sequence: str,
    modality: str,
    chemistry: str,
    dose_mg_kg: float,
    schedule: str
) -> Dict[str, float]:
    """
    Compute PK/PD-enriched features for GalNAc-siRNA.
    """
    sim = GalNAcPKPDSimulator(species="human", body_weight_kg=70.0)
    pk_metrics = sim.simulate(dose_mg_kg=dose_mg_kg, schedule=schedule)

    chem_modifiers = {
        "PS": 1.0,
        "PO": 0.4,
        "2'-OMe": 0.9,
        "2'-MOE": 0.95,
        "LNA": 0.85,
        "Mixed": 0.8,
    }
    stability_factor = chem_modifiers.get(chemistry, 0.8)

    auc_adj = pk_metrics["liver_auc_pmol_g_h"] * stability_factor
    cmax_adj = pk_metrics["liver_cmax_pmol_g"] * stability_factor  
    liver_half_life_h = pk_metrics["liver_half_life_h"]

    exposure_score = min(
        4.99,
        ((auc_adj * cmax_adj) / (liver_half_life_h + 1.0)) / 1e7
    )

    return {
        "liver_auc_pmol_g_h": auc_adj,
        "liver_cmax_pmol_g": cmax_adj,
        "liver_half_life_h": liver_half_life_h,
        "peak_knockdown_pct": pk_metrics["peak_knockdown_pct"],
        "exposure_score": exposure_score,
    }


def adjust_toxicity_risk_with_pkpd(
    oligotox_result: Dict[str, Any],
    sequence: str,
    modality: str,
    chemistry: str,
    dose_mg_kg: float,
    schedule: str
) -> Dict[str, Any]:
    """
    Adjust OligoTox risk predictions using PK/PD exposure metrics.
    """
    pkpd_features = compute_pkpd_enriched_features(
        sequence=sequence,
        modality=modality,
        chemistry=chemistry,
        dose_mg_kg=dose_mg_kg,
        schedule=schedule,
    )

    baseline_risks = oligotox_result["risk_scores"].copy()
    exposure_score = pkpd_features["exposure_score"]

    adjusted_risks = {}
    for endpoint, baseline_risk in baseline_risks.items():
        if endpoint == "hepatotoxicity":
            adjusted_risk = baseline_risk + _hepatotox_pkpd_adjustment(exposure_score)
        else:
            sensitivity = {
                "nephrotoxicity": 0.10,
                "thrombocytopenia": 0.15,
                "immune_activation": 0.12,
            }.get(endpoint, 0.10)
            adjustment_factor = 1.0 + sensitivity * (exposure_score - 1.0)
            adjusted_risk = baseline_risk * adjustment_factor

        adjusted_risks[endpoint] = max(0.0, min(1.0, adjusted_risk))

    result_adjusted = oligotox_result.copy()
    result_adjusted["risk_scores"] = adjusted_risks
    result_adjusted["pkpd_metrics"] = pkpd_features
    result_adjusted["dominant_endpoint"] = max(adjusted_risks.items(), key=lambda x: x[1])[0]
    return result_adjusted
