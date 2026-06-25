"""
PK/PD Bridge to OligoTox Toxicity Predictions
Phase 2A: Multiplicative risk adjustment based on liver exposure
"""
from typing import Dict, Any
from .simulator import GalNAcPKPDSimulator


def compute_pkpd_enriched_features(
    sequence: str,
    modality: str,
    chemistry: str,
    dose_mg_kg: float,
    schedule: str
) -> Dict[str, float]:
    """
    Compute PK/PD-enriched features for GalNAc-siRNA.
    
    Args:
        sequence: Oligonucleotide sequence
        modality: e.g., "siRNA"
        chemistry: e.g., "PS", "PO", "2'-OMe"
        dose_mg_kg: Dose in mg/kg
        schedule: Dosing schedule (e.g., "Q4W")
    
    Returns:
        Dict with PK/PD metrics and exposure score
    """
    # Initialize simulator
    sim = GalNAcPKPDSimulator(species="human", body_weight_kg=70.0)
    
    # Run PK/PD simulation
    pk_metrics = sim.simulate(dose_mg_kg=dose_mg_kg, schedule=schedule)
    
    # Chemistry stability modifiers
    chem_modifiers = {
        "PS": 1.0,     # Fully stabilized
        "PO": 0.4,     # Unstable
        "2'-OMe": 0.9,
        "2'-MOE": 0.95,
        "LNA": 0.85,
        "Mixed": 0.8
    }
    stability_factor = chem_modifiers.get(chemistry, 0.8)
    
    # Adjust exposure for chemistry
    auc_adj = pk_metrics["liver_auc_pmol_g_h"] * stability_factor
    cmax_adj = pk_metrics["liver_cmax_pmol_g"] * stability_factor
    
    # Compute composite exposure score (weighted AUC + Cmax)
    exposure_score = (0.7 * (auc_adj / 10000)) + (0.3 * (cmax_adj / 500))
    
    return {
        "liver_auc_pmol_g_h": auc_adj,
        "liver_cmax_pmol_g": cmax_adj,
        "liver_half_life_h": pk_metrics["liver_half_life_h"],
        "peak_knockdown_pct": pk_metrics["peak_knockdown_pct"],
        "exposure_score": exposure_score
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
    
    Phase 2A approach:
    - Multiplicative adjustment: baseline_risk * exposure_factor
    - Primarily affects hepatotoxicity (liver-targeted)
    - Other endpoints adjusted more conservatively
    
    Args:
        oligotox_result: Result dict from run_oligotox
        sequence: Oligo sequence
        modality: "siRNA", "ASO", etc.
        chemistry: Backbone chemistry
        dose_mg_kg: Dose
        schedule: Dosing schedule
    
    Returns:
        Updated result dict with PK/PD-adjusted risks and metrics
    """
    # Compute PK/PD metrics
    pkpd_features = compute_pkpd_enriched_features(
        sequence=sequence,
        modality=modality,
        chemistry=chemistry,
        dose_mg_kg=dose_mg_kg,
        schedule=schedule
    )
    
    # Extract baseline risk scores
    baseline_risks = oligotox_result["risk_scores"].copy()
    
    # Compute exposure adjustment factors
    exposure_score = pkpd_features["exposure_score"]
    
    # Define endpoint-specific sensitivity to liver exposure
    endpoint_sensitivity = {
        "hepatotoxicity": 0.8,      # Highly sensitive to liver exposure
        "nephrotoxicity": 0.2,      # Less sensitive (not liver-targeted)
        "thrombocytopenia": 0.3,    # Moderate (systemic exposure)
        "immune_activation": 0.25   # Moderate (systemic)
    }
    
    # Adjust each risk score
    adjusted_risks = {}
    for endpoint, baseline_risk in baseline_risks.items():
        sensitivity = endpoint_sensitivity.get(endpoint, 0.2)
        
        # Multiplicative adjustment: risk * (1 + sensitivity * exposure_score)
        adjustment_factor = 1.0 + (sensitivity * (exposure_score - 1.0))
        adjusted_risk = baseline_risk * adjustment_factor
        
        # Clamp to [0, 1]
        adjusted_risks[endpoint] = max(0.0, min(1.0, adjusted_risk))
    
    # Update result
    result_adjusted = oligotox_result.copy()
    result_adjusted["risk_scores"] = adjusted_risks
    result_adjusted["pkpd_metrics"] = pkpd_features
    
    # Recalculate dominant endpoint if needed
    dominant_endpoint = max(adjusted_risks.items(), key=lambda x: x[1])[0]
    result_adjusted["dominant_endpoint"] = dominant_endpoint
    
    return result_adjusted
