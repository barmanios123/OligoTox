from __future__ import annotations

import numpy as np
import pandas as pd


def _num(s: pd.Series, default: float = 0.0) -> pd.Series:
    return pd.to_numeric(s, errors="coerce").fillna(default)


def _flag_from_text(s: pd.Series, keywords: list[str]) -> pd.Series:
    x = s.fillna("").astype(str).str.lower()
    return x.apply(lambda v: 1.0 if any(k in v for k in keywords) else 0.0)


def build_aso_rich_features(df: pd.DataFrame) -> pd.DataFrame:
    out = pd.DataFrame(index=df.index)

    length_nt = _num(df["length_nt"]) if "length_nt" in df.columns else pd.Series(0.0, index=df.index)
    gc_content = _num(df["gc_content"]) if "gc_content" in df.columns else pd.Series(0.0, index=df.index)
    ps_fraction = _num(df["ps_fraction"]) if "ps_fraction" in df.columns else pd.Series(0.0, index=df.index)
    lna_fraction = _num(df["lna_fraction"]) if "lna_fraction" in df.columns else pd.Series(0.0, index=df.index)
    moe_fraction = _num(df["moe_fraction"]) if "moe_fraction" in df.columns else pd.Series(0.0, index=df.index)

    modality = df["modality"].fillna("").astype(str) if "modality" in df.columns else pd.Series("", index=df.index)
    chemistry = df["chemistry"].fillna("").astype(str) if "chemistry" in df.columns else pd.Series("", index=df.index)
    delivery = df["delivery"].fillna("").astype(str) if "delivery" in df.columns else pd.Series("", index=df.index)
    route = df["route"].fillna("").astype(str) if "route" in df.columns else pd.Series("", index=df.index)

    has_galnac = _num(df["has_galnac"]) if "has_galnac" in df.columns else pd.Series(0.0, index=df.index)
    lnp_formulated = _num(df["lnp_formulated"]) if "lnp_formulated" in df.columns else pd.Series(0.0, index=df.index)

    out["aso_length_nt"] = length_nt
    out["aso_gc_content"] = gc_content
    out["aso_ps_fraction"] = ps_fraction
    out["aso_lna_fraction"] = lna_fraction
    out["aso_moe_fraction"] = moe_fraction

    out["aso_ps_load"] = ps_fraction * length_nt
    out["aso_lna_load"] = lna_fraction * length_nt
    out["aso_moe_load"] = moe_fraction * length_nt
    out["aso_mod_total_fraction"] = (ps_fraction + lna_fraction + moe_fraction).clip(lower=0.0)
    out["aso_mod_total_load"] = out["aso_mod_total_fraction"] * length_nt

    out["aso_is_long_oligo"] = (length_nt >= 18).astype(float)
    out["aso_is_short_oligo"] = (length_nt <= 12).astype(float)
    out["aso_gc_high"] = (gc_content >= 0.60).astype(float)
    out["aso_gc_low"] = (gc_content <= 0.35).astype(float)
    out["aso_ps_high"] = (ps_fraction >= 0.75).astype(float)
    out["aso_lna_high"] = (lna_fraction >= 0.30).astype(float)
    out["aso_moe_high"] = (moe_fraction >= 0.30).astype(float)

    modality_l = modality.str.lower()
    chemistry_l = chemistry.str.lower()
    delivery_l = delivery.str.lower()
    route_l = route.str.lower()

    out["aso_modality_gapmer"] = modality_l.str.contains("gapmer", na=False).astype(float)
    out["aso_modality_sirna"] = modality_l.str.contains("sirna", na=False).astype(float)
    out["aso_modality_pmo"] = modality_l.str.contains("pmo", na=False).astype(float)
    out["aso_modality_sso"] = modality_l.str.contains("sso", na=False).astype(float)
    out["aso_modality_mrna"] = modality_l.str.contains("mrna", na=False).astype(float)

    out["aso_chem_lna"] = chemistry_l.str.contains("lna", na=False).astype(float)
    out["aso_chem_moe"] = chemistry_l.str.contains("2moe|moe", na=False).astype(float)
    out["aso_chem_2ome"] = chemistry_l.str.contains("2ome", na=False).astype(float)
    out["aso_chem_ps"] = chemistry_l.str.contains("ps", na=False).astype(float)
    out["aso_chem_pmo"] = chemistry_l.str.contains("pmo", na=False).astype(float)
    out["aso_chem_galnac"] = chemistry_l.str.contains("galnac", na=False).astype(float)

    out["aso_delivery_galnac"] = ((has_galnac > 0) | delivery_l.str.contains("galnac", na=False)).astype(float)
    out["aso_delivery_lnp"] = ((lnp_formulated > 0) | delivery_l.str.contains("lnp", na=False)).astype(float)
    out["aso_delivery_naked"] = delivery_l.str.contains("naked", na=False).astype(float)

    out["aso_route_iv"] = route_l.eq("iv").astype(float)
    out["aso_route_sc"] = route_l.eq("sc").astype(float)
    out["aso_route_im"] = route_l.eq("im").astype(float)
    out["aso_route_it"] = route_l.eq("it").astype(float)
    out["aso_route_in_vitro"] = route_l.eq("in_vitro").astype(float)

    out["aso_gapmer_like"] = (
        (out["aso_modality_gapmer"] == 1.0)
        | ((out["aso_chem_lna"] == 1.0) & (out["aso_chem_ps"] == 1.0) & (length_nt >= 14))
    ).astype(float)

    out["aso_galnac_ps_combo"] = ((out["aso_delivery_galnac"] == 1.0) & (out["aso_ps_high"] == 1.0)).astype(float)
    out["aso_lnp_gc_combo"] = ((out["aso_delivery_lnp"] == 1.0) & (out["aso_gc_high"] == 1.0)).astype(float)
    out["aso_lna_ps_combo"] = ((out["aso_chem_lna"] == 1.0) & (out["aso_chem_ps"] == 1.0)).astype(float)
    out["aso_moe_ps_combo"] = ((out["aso_chem_moe"] == 1.0) & (out["aso_chem_ps"] == 1.0)).astype(float)

    out["aso_risk_backbone_density"] = out["aso_ps_load"] / (length_nt.replace(0, np.nan)).fillna(1.0)
    out["aso_risk_mod_density"] = out["aso_mod_total_load"] / (length_nt.replace(0, np.nan)).fillna(1.0)
    out["aso_risk_lna_gc_interaction"] = out["aso_lna_fraction"] * out["aso_gc_content"]
    out["aso_risk_ps_length_interaction"] = out["aso_ps_fraction"] * out["aso_length_nt"]

    for col in out.columns:
        out[col] = pd.to_numeric(out[col], errors="coerce").fillna(0.0)

    return out
