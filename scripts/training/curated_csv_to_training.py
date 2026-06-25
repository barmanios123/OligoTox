from __future__ import annotations

import pandas as pd

from oligotox.schemas.data_models import (
    OligoInput,
    DeliveryInput,
    OligoToxInput,
    PlateletComplement,
    OmicsReadout,
    Modality,
    DeliveryVehicle,
    TargetClass,
)
from oligotox.features.engineer import extract_features
import re

RAW_PATH = "data/curated/master/curated_oligo_tox_template_v5.csv"
OUT_PATH = "data/processed/training_dataset.csv"

LABEL_COLS = [
    "hepatocellular_injury",
    "hepatic_vacuolation",
    "renal_tubular_injury",
    "renal_accumulation",
    "thrombocytopenia",
    "coagulation_effect",
    "complement_activation",
    "innate_immune_activation",
    "cytokine_release",
    "kupffer_cell_activation",
    "cns_toxicity",
    "offtarget_gene_silencing",
    "injection_site_reaction",
]

def to_bool(x):
    if pd.isna(x):
        return False
    if isinstance(x, str):
        return x.strip().lower() in {"1", "true", "yes", "y"}
    return bool(x)


def _float_or_none(r, col):
    try:
        val = r[col]
    except KeyError:
        return None
    if pd.isna(val):
        return None
    try:
        return float(val)
    except Exception:
        return None


    if pd.isna(x):
        return False
    if isinstance(x, str):
        return x.strip().lower() in {"1", "true", "yes", "y"}
    return bool(x)


def _norm_str(x):
    if x is None:
        return ""
    s = str(x).strip()
    if s.lower() in ("", "nan", "none", "null"):
        return ""
    return s

def _slug(x):
    s = _norm_str(x).lower()
    s = s.replace("-", "_").replace(" ", "_").replace("/", "_")
    s = re.sub(r"_+", "_", s)
    return s.strip("_")

def normalize_vehicle(x):
    s = _slug(x)
    mapping = {
        "": "none",
        "none": "none",
        "naked": "none",
        "n_a": "none",
        "na": "none",
        "n_a_": "none",
        "n/a": "none",
        "null": "none",
        "nan": "none",
        "lnp": "LNP",
        "lipid_nanoparticle": "LNP",
        "lipid_nanoparticles": "LNP",
        "galnac": "GalNAc",
        "triantennary_galnac": "GalNAc",
        "saline": "saline",
    }
    return mapping.get(s, _norm_str(x))

def normalize_route(x):
    s = _slug(x)
    mapping = {
        "iv": "iv",
        "intravenous": "iv",
        "sc": "sc",
        "subcutaneous": "sc",
        "sub_q": "sc",
        "sq": "sc",
        "it": "it",
        "intrathecal": "it",
        "im": "im",
        "intramuscular": "im",
        "ip": "ip",
        "intraperitoneal": "ip",
        "po": "po",
        "oral": "po",
    }
    return mapping.get(s, _norm_str(x))

def normalize_modality(modality_raw, vehicle_raw=""):
    m = _slug(modality_raw)
    v = _slug(vehicle_raw)

    if m in {"aso", "gapmer", "aso_gapmer", "rnase_h_gapmer"}:
        return "ASO_gapmer"

    if m in {"aso_splice", "splice", "steric_block", "steric_blocking", "ssao"}:
        return "ASO_splice"

    if m in {"pmo", "phosphorodiamidate_morpholino"}:
        return "PMO"

    if m in {"sirna_lnp"}:
        return "siRNA_LNP"

    if m in {"sirna_galnac", "galnac_sirna"}:
        return "siRNA_GalNAc"

    if m in {"sirna", "si_rna"}:
        if v == "galnac":
            return "siRNA_GalNAc"
        if v == "lnp":
            return "siRNA_LNP"
        return "siRNA"

    if not m and v == "galnac":
        return "siRNA_GalNAc"

    if not m and v == "lnp":
        return "siRNA_LNP"

    return _norm_str(modality_raw)

def resolve_enum_inputs(row):
    vehicle = normalize_vehicle(row.get("vehicle", ""))
    route = normalize_route(row.get("route", ""))
    modality = normalize_modality(row.get("modality", ""), vehicle)
    return modality, vehicle, route

def main():
    df = pd.read_csv(RAW_PATH)
    rows = []

    for _, r in df.iterrows():
        _modality, _vehicle, _route = resolve_enum_inputs(r)

        try:
            modality = Modality(_modality)
            vehicle = DeliveryVehicle(_vehicle)
        except ValueError as e:
            print(
                f"Skipping row due to enum mismatch: "
                f"candidate_id={r.get('candidate_id', '')} "
                f"modality={_modality} vehicle={_vehicle} route={_route} error={e}"
            )
            continue

        oligo = OligoInput(
            sequence=str(r["sequence"]),
            modality=modality,
            target_gene=None if pd.isna(r["target_gene"]) else str(r["target_gene"]),
            target_class=TargetClass.CODING_RNA,
            is_human_specific_sequence=to_bool(r["is_human_specific_sequence"]),
            target_is_lncrna_hub=to_bool(r["target_is_lncrna_hub"]),
            micropeptide_co_depletion_risk=to_bool(r["micropeptide_co_depletion_risk"]),
        )

        delivery = DeliveryInput(
            vehicle=vehicle,
            route=_route,
            dose_mg_kg=float(r["dose_mg_kg"]),
            repeat_dose=to_bool(r["repeat_dose"]),
            peg_lipid_mol_pct=float(r["peg_lipid_mol_pct"]) if not pd.isna(r["peg_lipid_mol_pct"]) else None,
            particle_size_nm=float(r["particle_size_nm"]) if not pd.isna(r["particle_size_nm"]) else None,
            pdi=float(r["pdi"]) if not pd.isna(r["pdi"]) else None,
            galnac_valency=int(r["galnac_valency"]) if not pd.isna(r["galnac_valency"]) else None,
        )

        platelet = PlateletComplement(
            platelet_count_baseline=_float_or_none(r, "platelet_count_baseline"),
            platelet_count_treated=_float_or_none(r, "platelet_count_treated"),
            pf4_ps_elisa_od=_float_or_none(r, "pf4_ps_elisa_od"),
            c3a_pg_ml=_float_or_none(r, "c3a_pg_ml"),
            c5a_pg_ml=_float_or_none(r, "c5a_pg_ml"),
            sc5b9_ng_ml=_float_or_none(r, "sc5b9_ng_ml"),
            complement_activation_ex_vivo=to_bool(r["complement_activation_ex_vivo"]),
            platelet_c3d_c4d_positive=to_bool(r["platelet_c3d_c4d_positive"]),
        )

        omics = OmicsReadout(
            alt_iu_l=_float_or_none(r, "alt_iu_l"),
            ast_iu_l=_float_or_none(r, "ast_iu_l"),
            bun_mg_dl=_float_or_none(r, "bun_mg_dl"),
            creatinine_mg_dl=_float_or_none(r, "creatinine_mg_dl"),
            cytokine_il6=_float_or_none(r, "cytokine_il6"),
            cytokine_tnfa=_float_or_none(r, "cytokine_tnfa"),
            cytokine_ifna=_float_or_none(r, "cytokine_ifna"),
            pathway_scores={
                "pw_nfkb": _float_or_none(r, "pw_nfkb"),
                "pw_ifn": _float_or_none(r, "pw_ifn"),
                "pw_apop": _float_or_none(r, "pw_apop"),
                "pw_comp": _float_or_none(r, "pw_comp"),
                "pw_plt": _float_or_none(r, "pw_plt"),
                "pw_kyn": _float_or_none(r, "pw_kyn"),
            },
        )

        inp = OligoToxInput(
            oligo=oligo,
            delivery=delivery,
            omics=omics,
        )

        feats = extract_features(inp)

        out = {
            "candidate_id": r["candidate_id"],
            "source": r["source"],
            "notes": r["notes"],
            "series_id": r.get("series_id", ""),
        }
        out.update(feats)

        for c in LABEL_COLS:
            val = r.get(c, 0)
            try:
                out[c] = int(val)
            except (TypeError, ValueError):
                out[c] = 0

        rows.append(out)

    out_df = pd.DataFrame(rows)
    out_df.to_csv(OUT_PATH, index=False)
    print(out_df.head())
    print(out_df.columns.tolist())
    print(out_df.shape)
    print(f"Wrote {OUT_PATH}")

if __name__ == "__main__":
    main()