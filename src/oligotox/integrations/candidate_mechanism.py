from __future__ import annotations

from typing import Any, Dict, List

import pandas as pd

from oligotox.integrations.seed_match import score_seed_matches
from oligotox.integrations.tlr_immune import score_tlr_innate
from oligotox.integrations.gsea_pathway import score_pathways


def normalize_sequence(seq: str) -> str:
    return "".join((seq or "").strip().upper().split()).replace("T", "U")


def infer_deg_list_from_endpoints(top_endpoints: List[str]) -> List[str]:
    degs: List[str] = []

    endpoint_to_degs = {
        "innate_immune_activation": ["IL6", "CXCL8", "NFKBIA", "ISG15", "OAS1"],
        "cytokine_release":         ["IL6", "CXCL8", "NFKBIA", "ISG15", "OAS1"],
        "complement_activation":    ["C3", "C5", "CFB"],
        "thrombocytopenia":         ["PF4", "VWF", "SERPINE1"],
        "coagulation_effect":       ["PF4", "VWF", "SERPINE1"],
        "hepatocellular_injury":    ["BAX", "DDIT3", "ATF4"],
        "renal_tubular_injury":     ["KIM1", "LCN2"],
    }

    for ep in top_endpoints:
        if ep in endpoint_to_degs:
            degs.extend(endpoint_to_degs[ep])

    if not degs:
        degs = ["IL6", "CXCL8", "NFKBIA"]

    return list(dict.fromkeys(degs))


def build_candidate_mechanistic_evidence(
    payload: Dict[str, Any],
    risk_df: pd.DataFrame,
) -> Dict[str, Any]:
    seq = normalize_sequence(payload.get("sequence", ""))
    ps_fraction = float(payload.get("ps_fraction", 0.0) or 0.0)

    # Seed-match scoring
    try:
        seed_res = score_seed_matches(seq)
    except Exception as e:
        raise RuntimeError(f"seed_match failed: {e}") from e

    # Innate immune scoring
    try:
        tlr_res = score_tlr_innate(
            seq,
            ps_fraction=ps_fraction,
            is_dsrna=bool(seed_res.dsrna_probability >= 0.5),
        )
    except Exception as e:
        raise RuntimeError(f"tlr_immune failed: {e}") from e

    # Pathway scoring based on inferred DEG list
    top_endpoints = risk_df.head(5)["endpoint"].tolist()
    degs = infer_deg_list_from_endpoints(top_endpoints)

    try:
        pathway_res = score_pathways(degs)
    except Exception as e:
        raise RuntimeError(f"gsea_pathway failed: {e}") from e

    # Build plain-English bullets
    bullets = [
        f"D2 seed risk {seed_res.d2_risk_contribution:.3f} with "
        f"{seed_res.n_perfect_seed_hits} predicted 7-mer hits and "
        f"{seed_res.n_partial_seed_hits} 6-mer hits.",
        f"D3 innate score {tlr_res.composite_innate_score:.3f}; "
        f"dominant receptor {tlr_res.dominant_receptor}.",
    ]

    if pathway_res.top_pathways:
        top_pw = pathway_res.top_pathways[0]
        bullets.append(
            f"Top pathway context: {top_pw.pathway} ({top_pw.enrichment_score:.3f}) "
            f"with overlap {', '.join(top_pw.overlap_genes[:3])}."
        )

    pathway_rows = [
        {
            "pathway": p.pathway,
            "score": round(float(p.enrichment_score), 3),
            "overlap_genes": ", ".join(p.overlap_genes[:5]),
            "n_overlap": int(p.n_overlap),
        }
        for p in pathway_res.top_pathways[:5]
    ]

    return {
        "seed": {
            "seed_region": seed_res.seed_region,
            "seed_rc": seed_res.seed_rc,
            "n_perfect_seed_hits": int(seed_res.n_perfect_seed_hits),
            "n_partial_seed_hits": int(seed_res.n_partial_seed_hits),
            "seed_match_score": float(seed_res.seed_match_score),
            "dsrna_probability": float(seed_res.dsrna_probability),
            "d2_risk_contribution": float(seed_res.d2_risk_contribution),
        },
        "innate": {
            "tlr7_score": float(tlr_res.tlr7_score),
            "tlr8_score": float(tlr_res.tlr8_score),
            "tlr9_score": float(tlr_res.tlr9_score),
            "rigi_score": float(tlr_res.rigi_score),
            "mda5_score": float(tlr_res.mda5_score),
            "composite_innate_score": float(tlr_res.composite_innate_score),
            "dominant_receptor": tlr_res.dominant_receptor,
            "mechanistic_note": tlr_res.mechanistic_note,
        },
        "pathways": pd.DataFrame(pathway_rows),
        "bullets": bullets,
        "deg_list": degs,
    }