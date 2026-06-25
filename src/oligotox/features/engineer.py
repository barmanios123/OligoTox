"""Feature Engineering Layer"""
import math
import re
from typing import Dict


from oligotox.schemas.data_models import OligoToxInput, Modality, DeliveryVehicle


def _safe_fold(val, baseline):
    if val is None:
        return 0.0
    try:
        return max(val / baseline - 1.0, 0.0)
    except Exception:
        return 0.0


def _safe_val(val):
    if val is None:
        return 0.0
    try:
        return float(val)
    except Exception:
        return 0.0



def extract_features(inp: OligoToxInput) -> Dict:
    o = inp.oligo
    d = inp.delivery
    pc = inp.platelet_complement
    om = inp.omics

    seq = (o.sequence or "").upper().replace("T", "U")
    L = max(len(seq), 1)

    # Sequence features
    gc = (seq.count("G") + seq.count("C")) / L
    cpg = seq.replace("U", "T").count("CG") / max(L - 1, 1)
    poly_u = len(re.findall(r"U{4,}", seq)) / L
    g4 = len(re.findall(r"G{3,}", seq)) / L
    palindromes = sum(
        1
        for i in range(L - 5)
        if seq[i:i + 6] == seq[i:i + 6].translate(str.maketrans("AUGC", "UACG"))[::-1]
    )
    self_comp = min(palindromes / max(L - 5, 1) * 4, 1.0)

    # Chemistry features
    ps_frac = 0.0
    lna_frac = 0.0
    if o.nucleotides and len(o.nucleotides) > 0:
        ps_frac = sum(1 for n in o.nucleotides if n.backbone_mod == "PS") / len(o.nucleotides)
        lna_frac = sum(1 for n in o.nucleotides if n.sugar_mod == "LNA") / len(o.nucleotides)
    else:
        defaults = {
            Modality.ASO_GAPMER: (0.90, 0.25),
            Modality.ASO_SSO: (0.80, 0.20),
            Modality.SIRNA_LNP: (0.20, 0.00),
            Modality.SIRNA_GALNAC: (0.20, 0.00),
            Modality.SIRNA_GENERIC: (0.20, 0.00),
            Modality.LNCRNA_TARGETING: (0.75, 0.15),
            Modality.MIRNA_MIMIC: (0.20, 0.00),
            Modality.MIRNA_INHIBITOR: (0.70, 0.10),
            Modality.CIRCRNA: (0.10, 0.00),
            Modality.PMO: (0.10, 0.00),
            Modality.APTAMER: (0.05, 0.00),
            Modality.ANTI_MIR: (0.70, 0.10),
            Modality.OTHER: (0.50, 0.05),
        }
        ps_frac, lna_frac = defaults.get(o.modality, (0.50, 0.05))

    lna_5flank = lna_frac * 0.6
    ps_entropy = -(
        ps_frac * math.log(ps_frac + 1e-9)
        + (1 - ps_frac) * math.log(1 - ps_frac + 1e-9)
    )

    # Delivery features
    route = (d.route or "").lower()
    is_lnp = d.vehicle == DeliveryVehicle.LNP
    is_galnac = d.vehicle == DeliveryVehicle.GALNAC
    is_it = ("intrathecal" in route) or (route == "it")
    is_iv = ("intravenous" in route) or (route == "iv")
    peg_pct = d.peg_lipid_mol_pct or 0.0
    particle_size = d.particle_size_nm or 0.0
    pdi_val = d.pdi or 0.0
    galnac_val = d.galnac_valency or 0
    dose = d.dose_mg_kg or 0.0
    repeat = float(d.repeat_dose)

    # Platelet/complement biomarkers
    plt_drop = max(
        0.0,
        (pc.platelet_count_baseline - pc.platelet_count_treated)
        / max(pc.platelet_count_baseline, 1),
    )
    pf4_od = pc.pf4_ps_elisa_od
    c3a = pc.c3a_pg_ml
    c5a = pc.c5a_pg_ml
    sc5b9 = pc.sc5b9_ng_ml
    carpa = float(pc.complement_activation_ex_vivo)
    c3d_pos = float(pc.platelet_c3d_c4d_positive)

    # Omics features
    alt_fold = _safe_fold(om.alt_iu_l, 40.0)
    ast_fold = _safe_fold(om.ast_iu_l, 40.0)
    bun_fold = _safe_fold(om.bun_mg_dl, 25.0)
    cr_fold = _safe_fold(om.creatinine_mg_dl, 1.2)
    il6 = _safe_val(om.cytokine_il6)
    tnfa = _safe_val(om.cytokine_tnfa)
    ifna = _safe_val(om.cytokine_ifna)

    pw_nfkb = om.pathway_scores.get("pw_nfkb", om.pathway_scores.get("nfkb", 0.0))
    pw_ifn = om.pathway_scores.get("pw_interferon", om.pathway_scores.get("ifn", 0.0))
    pw_apop = om.pathway_scores.get("pw_apoptosis", om.pathway_scores.get("apoptosis", 0.0))
    pw_comp = om.pathway_scores.get("pw_complement", om.pathway_scores.get("complement", 0.0))
    pw_plt = om.pathway_scores.get("pw_platelet", om.pathway_scores.get("platelet", 0.0))
    pw_kyn = om.pathway_scores.get("pw_kynurenine", om.pathway_scores.get("kyn", 0.0))

    # Target/biology flags
    is_lncrna_hub = float(o.target_is_lncrna_hub)
    micropeptide_risk = float(o.micropeptide_co_depletion_risk)
    human_specific = float(o.is_human_specific_sequence)

    return dict(
        seq_length=L,
        gc_content=gc,
        cpg_density=cpg,
        poly_u_density=poly_u,
        g4_motif_density=g4,
        self_complementarity=self_comp,
        ps_fraction=ps_frac,
        lna_fraction=lna_frac,
        lna_5flank=lna_5flank,
        ps_entropy=ps_entropy,
        is_lnp=float(is_lnp),
        is_galnac=float(is_galnac),
        is_it=float(is_it),
        is_iv=float(is_iv),
        peg_pct=peg_pct,
        particle_size=particle_size,
        pdi=pdi_val,
        galnac_valency=float(galnac_val),
        dose_mg_kg=dose,
        repeat_dose=repeat,
        platelet_drop_fraction=plt_drop,
        pf4_od=pf4_od,
        c3a=c3a,
        c5a=c5a,
        sc5b9=sc5b9,
        carpa=carpa,
        c3d_positive=c3d_pos,
        alt_fold=alt_fold,
        ast_fold=ast_fold,
        bun_fold=bun_fold,
        cr_fold=cr_fold,
        il6=il6,
        tnfa=tnfa,
        ifna=ifna,
        pw_nfkb=pw_nfkb,
        pw_ifn=pw_ifn,
        pw_apop=pw_apop,
        pw_comp=pw_comp,
        pw_plt=pw_plt,
        pw_kyn=pw_kyn,
        is_lncrna_hub=is_lncrna_hub,
        micropeptide_risk=micropeptide_risk,
        human_specific=human_specific,
    )