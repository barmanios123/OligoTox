from __future__ import annotations

from typing import Any, Dict

from oligotox.config.paths import (
    PROCESSED_DATA_DIR,
    MODEL_ARTIFACTS_DIR,
    CONFIG_ARTIFACTS_DIR,
)
from oligotox.pipelines.score_candidate import score_candidate
from oligotox.omics import (
    parse_omics_payload,
    apply_adjustments,
    recommend_assays,
    build_omics_explanation,
)


def _coerce_bool(v: Any) -> int:
    return int(bool(v))


def _safe_get(obj: Any, name: str, default=None):
    if isinstance(obj, dict):
        return obj.get(name, default)
    return getattr(obj, name, default)


def build_feature_dict_from_input(inp: Any) -> Dict[str, Any]:
    seq = _safe_get(inp, "sequence", "") or ""
    seq_upper = seq.upper()

    g = seq_upper.count("G")
    c = seq_upper.count("C")
    seq_len = max(len(seq_upper), 1)
    gc_content = (g + c) / seq_len

    route = str(_safe_get(inp, "route", "") or "").upper()
    delivery_vehicle = str(_safe_get(inp, "delivery_vehicle", "") or "").lower()

    return {
        "seq_length": len(seq_upper),
        "gc_content": gc_content,
        "cpg_density": seq_upper.count("CG") / seq_len,
        "poly_u_density": seq_upper.count("UUUU") / seq_len,
        "g4_motif_density": seq_upper.count("GGG") / seq_len,
        "self_complementarity": _safe_get(inp, "self_complementarity", 0.0),
        "ps_fraction": _safe_get(inp, "ps_fraction", 0.0),
        "lna_fraction": _safe_get(inp, "lna_fraction", 0.0),
        "lna_5flank": _safe_get(inp, "lna_5flank", 0.0),
        "ps_entropy": _safe_get(inp, "ps_entropy", 0.0),
        "is_lnp": _coerce_bool(_safe_get(inp, "is_lnp", False) or delivery_vehicle == "lnp"),
        "is_galnac": _coerce_bool(_safe_get(inp, "is_galnac", False) or delivery_vehicle == "galnac"),
        "is_it": _coerce_bool(_safe_get(inp, "is_it", False) or route == "IT"),
        "is_iv": _coerce_bool(_safe_get(inp, "is_iv", False) or route == "IV"),
        "peg_pct": _safe_get(inp, "peg_pct", 0.0),
        "particle_size": _safe_get(inp, "particle_size", 0.0),
        "pdi": _safe_get(inp, "pdi", 0.0),
        "galnac_valency": _safe_get(inp, "galnac_valency", 0.0),
        "dose_mg_kg": _safe_get(inp, "dose_mg_kg", 0.0),
        "repeat_dose": _coerce_bool(_safe_get(inp, "repeat_dose", False)),
        "platelet_drop_fraction": _safe_get(inp, "platelet_drop_fraction", 0.0),
        "pf4_od": _safe_get(inp, "pf4_od", 0.0),
        "c3a": _safe_get(inp, "c3a", 0.0),
        "c5a": _safe_get(inp, "c5a", 0.0),
        "sc5b9": _safe_get(inp, "sc5b9", 0.0),
        "carpa": _safe_get(inp, "carpa", 0.0),
        "c3d_positive": _coerce_bool(_safe_get(inp, "c3d_positive", False)),
        "alt_fold": _safe_get(inp, "alt_fold", 0.0),
        "ast_fold": _safe_get(inp, "ast_fold", 0.0),
        "bun_fold": _safe_get(inp, "bun_fold", 0.0),
        "cr_fold": _safe_get(inp, "cr_fold", 0.0),
        "il6": _safe_get(inp, "il6", 0.0),
        "tnfa": _safe_get(inp, "tnfa", 0.0),
        "ifna": _safe_get(inp, "ifna", 0.0),
        "pw_nfkb": _safe_get(inp, "pw_nfkb", 0.0),
        "pw_ifn": _safe_get(inp, "pw_ifn", 0.0),
        "pw_apop": _safe_get(inp, "pw_apop", 0.0),
        "pw_comp": _safe_get(inp, "pw_comp", 0.0),
        "pw_plt": _safe_get(inp, "pw_plt", 0.0),
        "pw_kyn": _safe_get(inp, "pw_kyn", 0.0),
        "is_lncrna_hub": _coerce_bool(_safe_get(inp, "is_lncrna_hub", False)),
        "micropeptide_risk": _safe_get(inp, "micropeptide_risk", 0.0),
        "human_specific": _coerce_bool(_safe_get(inp, "human_specific", False)),
    }


def _compute_mechanistic_proxies(inp: Any) -> Dict[str, Any]:
    seq = str(_safe_get(inp, "sequence", "") or "").upper().replace("T", "U")
    chemistry = str(_safe_get(inp, "chemistry", "") or "")
    conjugate = str(_safe_get(inp, "conjugate", "") or "")
    n = max(len(seq), 1)
    gc = (seq.count("G") + seq.count("C")) / n
    u_frac = seq.count("U") / n
    cpg = seq.count("CG")
    is_ps = "PS" in chemistry.upper()
    is_galnac = "GALNAC" in conjugate.upper()
    is_lna = "LNA" in chemistry.upper()

    return {
        "gc_fraction": round(gc, 4),
        "u_fraction": round(u_frac, 4),
        "cpg_count": cpg,
        "hybridization_risk_proxy": round(min(1.0, gc * 0.35 + cpg * 0.12), 4),
        "tlr78_risk_proxy": round(min(1.0, u_frac * 1.43), 4),
        "platelet_ps_risk_proxy": round(0.6857 if is_ps else 0.1, 4),
        "hepatic_gapmer_risk_proxy": round(0.45 if is_lna else 0.21, 4),
        "renal_accumulation_proxy": round(0.65 if is_galnac else 0.35, 4),
        "length_nt": len(seq),
    }


def _fallback_assays(mech: Dict[str, Any]) -> list[str]:
    assays = []
    if mech.get("platelet_ps_risk_proxy", 0) > 0.4:
        assays.append("Platelet activation and PF4-associated thrombocytopenia follow-up")
    if mech.get("renal_accumulation_proxy", 0) > 0.4:
        assays.append("Renal proximal tubule accumulation / injury follow-up")
    if mech.get("hybridization_risk_proxy", 0) > 0.15 or mech.get("hepatic_gapmer_risk_proxy", 0) > 0.3:
        assays.append("RNase H1 / hybridization off-target assessment")
    if mech.get("tlr78_risk_proxy", 0) > 0.25:
        assays.append("TLR7/8 innate immune activation assay (PBMC cytokine panel)")
    if mech.get("gc_fraction", 0) > 0.5 or mech.get("hepatic_gapmer_risk_proxy", 0) > 0.2:
        assays.append("ALT / AST and hepatocyte stress biomarker follow-up")
    return assays


def run_oligotox(inp: Any, omics_payload: dict | None = None) -> Dict[str, Any]:
    feature_dict = build_feature_dict_from_input(inp)

    base = score_candidate(
        candidate_features=feature_dict,
        processed_dir=PROCESSED_DATA_DIR,
        model_dir=MODEL_ARTIFACTS_DIR,
        config_dir=CONFIG_ARTIFACTS_DIR,
    )

    per_endpoint = base.get("per_endpoint", {}) if isinstance(base, dict) else {}
    if not isinstance(per_endpoint, dict):
        raise ValueError("score_candidate returned invalid per_endpoint structure")

    base_probs = {k: v["prob"] for k, v in per_endpoint.items() if isinstance(v, dict) and "prob" in v}

    omics_readout = parse_omics_payload(omics_payload) if omics_payload else None
    adjusted_probs, omics_adjustment = apply_adjustments(base_probs, omics_readout)

    for endpoint, prob in adjusted_probs.items():
        if endpoint in per_endpoint:
            per_endpoint[endpoint]["base_prob"] = per_endpoint[endpoint].get("prob")
            per_endpoint[endpoint]["prob"] = round(float(prob), 4)
            thr = per_endpoint[endpoint].get("threshold")
            per_endpoint[endpoint]["predicted_positive"] = bool(prob >= thr) if thr is not None else False

    mech = _compute_mechanistic_proxies(inp)
    proxy_boosts = {
        "thrombocytopenia": mech.get("platelet_ps_risk_proxy", 0.0),
        "renal_accumulation": mech.get("renal_accumulation_proxy", 0.0),
        "hepatic_vacuolation": mech.get("hepatic_gapmer_risk_proxy", 0.0),
        "innate_immune_activation": mech.get("tlr78_risk_proxy", 0.0),
        "complement_activation": mech.get("hybridization_risk_proxy", 0.0),
    }

    for ep, proxy_val in proxy_boosts.items():
        if ep in per_endpoint and isinstance(per_endpoint[ep], dict):
            current = float(per_endpoint[ep].get("prob", 0.0))
            blended = round(max(current, proxy_val * 0.6), 4)
            per_endpoint[ep]["prob"] = blended
            thr = per_endpoint[ep].get("threshold")
            per_endpoint[ep]["predicted_positive"] = bool(blended >= thr) if thr is not None else False
            per_endpoint[ep]["proxy_support"] = round(proxy_val, 4)

    final_probs = {
        k: v["prob"]
        for k, v in per_endpoint.items()
        if isinstance(v, dict) and "prob" in v
    }
    dominant_endpoint = max(final_probs, key=final_probs.get) if final_probs else None

    omics_assays = recommend_assays(omics_readout, final_probs) or []
    merged_assays = omics_assays if omics_assays else _fallback_assays(mech)

    base["dominant_endpoint"] = dominant_endpoint
    base["mechanistic_proxies"] = mech
    base["omics"] = {
        "adjustment_summary": omics_adjustment,
        "explanation": build_omics_explanation(omics_readout),
        "recommended_assays": merged_assays,
    }
    return base


def run_oligotox_unified(inp: Any) -> Dict[str, Any]:
    feature_dict = build_feature_dict_from_input(inp)

    sequence = _safe_get(inp, "sequence", "")
    if sequence:
        feature_dict["sequence"] = sequence

    v03 = score_candidate(
        candidate_features=feature_dict,
        processed_dir=PROCESSED_DATA_DIR,
        model_dir=MODEL_ARTIFACTS_DIR,
        config_dir=CONFIG_ARTIFACTS_DIR,
    )

    trans_cfg = _safe_get(inp, "translational_config", None)
    if not trans_cfg:
        trans_cfg = {
            "features": feature_dict,
            "species": _safe_get(inp, "species", "rat+NHP"),
            "disease_context": _safe_get(inp, "disease_context", "metabolic_liver"),
            "has_surrogate": _safe_get(inp, "has_surrogate", True),
            "dart_relevant": _safe_get(inp, "dart_relevant", False),
        }

    n1_cfg = _safe_get(inp, "n_of_1_config", None)
    if not n1_cfg:
        n1_cfg = {
            "target_gene": _safe_get(inp, "target_gene", "UNKNOWN"),
            "disease_context": _safe_get(inp, "disease_context", "rare genetic disease"),
            "timeline_weeks": _safe_get(inp, "timeline_weeks", 12),
        }

    v04 = {}

    return {
        "v03_model": v03,
        "v04_decision_support": v04,
    }


def score_candidate_payload(payload):
    import pandas as pd

    result = score_candidate(
        candidate_features=build_feature_dict_from_input(payload),
        processed_dir=PROCESSED_DATA_DIR,
        model_dir=MODEL_ARTIFACTS_DIR,
        config_dir=CONFIG_ARTIFACTS_DIR,
    )

    if not isinstance(result, dict):
        raise TypeError(f"Unsupported score_candidate output type: {type(result)}")

    per_endpoint = result.get("per_endpoint", {})
    if not isinstance(per_endpoint, dict) or not per_endpoint:
        raise ValueError("score_candidate_payload could not find per_endpoint scores")

    d1_map = {
        "hepatocellular_injury": 0.70,
        "hepatic_vacuolation": 0.55,
        "offtarget_gene_silencing": 0.60,
        "cns_toxicity": 0.35,
    }
    d2_map = {
        "innate_immune_activation": 0.75,
        "cytokine_release": 0.80,
        "complement_activation": 0.85,
        "kupffer_cell_activation": 0.60,
        "injection_site_reaction": 0.45,
    }
    d3_map = {
        "renal_tubular_injury": 0.80,
        "renal_accumulation": 0.85,
        "thrombocytopenia": 0.90,
        "coagulation_effect": 0.80,
        "cns_toxicity": 0.25,
    }

    confidence_map = {
        "low": "Low",
        "elevated": "Medium",
        "high": "High",
    }

    rows = []
    for endpoint, val in per_endpoint.items():
        if not isinstance(val, dict):
            continue

        prob = float(val.get("prob", 0.0))
        risk_band = str(val.get("risk_band", "low")).lower()
        confidence = confidence_map.get(risk_band, "Medium")

        rows.append(
            {
                "endpoint": str(endpoint),
                "risk_score": prob,
                "confidence": confidence,
                "predicted_positive": bool(val.get("predicted_positive", False)),
                "threshold": float(val.get("threshold", 0.0) or 0.0),
                "risk_band": str(val.get("risk_band", "low")),
                "D1": round(prob * d1_map.get(endpoint, 0.10), 3),
                "D2": round(prob * d2_map.get(endpoint, 0.10), 3),
                "D3": round(prob * d3_map.get(endpoint, 0.10), 3),
            }
        )

    out = pd.DataFrame(rows)
    if out.empty:
        raise ValueError("score_candidate_payload produced no endpoint rows")

    out["endpoint"] = out["endpoint"].astype(str)
    out["confidence"] = out["confidence"].astype(str)
    out["risk_band"] = out["risk_band"].astype(str)
    out["predicted_positive"] = out["predicted_positive"].astype(bool)

    numeric_cols = ["risk_score", "threshold", "D1", "D2", "D3"]
    for col in numeric_cols:
        out[col] = pd.to_numeric(out[col], errors="coerce").fillna(0.0)

    keep = [
        "endpoint",
        "risk_score",
        "confidence",
        "D1",
        "D2",
        "D3",
        "predicted_positive",
        "threshold",
        "risk_band",
    ]
    out = out[keep].copy()
    return out.sort_values("risk_score", ascending=False).reset_index(drop=True)