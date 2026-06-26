import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from oligotox.inference.api import run_oligotox


def test_omics_decision_layer():
    candidate = {
        "sequence": "GCUAUGCUAGCUAGCUAGCUAGCUAGC",
        "is_lnp": True,
        "is_galnac": False,
        "is_iv": True,
        "is_it": False,
        "dose_mg_kg": 1.5,
        "repeat_dose": False,
        "ps_fraction": 0.7,
        "lna_fraction": 0.2,
        "particle_size": 85,
        "pdi": 0.18,
        "peg_pct": 1.5,
        "c3a": 0.2,
        "c5a": 0.1,
        "sc5b9": 0.25,
        "alt_fold": 1.1,
        "ast_fold": 1.0,
        "bun_fold": 1.0,
        "cr_fold": 1.0,
        "il6": 0.1,
        "tnfa": 0.05,
        "ifna": 0.05,
        "pw_nfkb": 0.2,
        "pw_ifn": 0.1,
        "pw_apop": 0.1,
        "pw_comp": 0.25,
        "pw_plt": 0.15,
        "pw_kyn": 0.05,
        "micropeptide_risk": 0.1,
        "human_specific": True,
    }

    omics_payload = {
        "candidate_id": "demo_001",
        "modality": "transcriptomics+proteomics",
        "tissue": "liver",
        "species": "human",
        "system": "primary_hepatocyte",
        "platform": "rna_seq",
        "timepoint_hr": 24,
        "dose_level": "mid",
        "dose_value": 3.0,
        "dose_unit": "uM",
        "replicate_count": 3,
        "coverage_score": 0.9,
        "quality_score": 0.85,
        "concordance_score": 0.75,
        "translational_relevance": 0.8,
        "pw_interferon": 0.72,
        "pw_nfkb": 0.66,
        "pw_tlr_innate": 0.61,
        "pw_complement": 0.22,
        "pw_platelet": 0.18,
        "pw_apoptosis": 0.71,
        "pw_hepatocyte_injury": 0.74,
        "pw_kupffer_activation": 0.59,
        "pw_rnaseh1_long_transcript_risk": 0.64,
        "pw_hybridization_offtarget": 0.55,
        "sentinel_ISG15": 0.88,
        "sentinel_CXCL10": 0.81,
        "sentinel_IL6": 0.62,
        "sentinel_PF4": 0.21,
        "sentinel_C3": 0.18,
        "sentinel_RNASEH1_DEP_LONGTX": 0.67,
        "sentinel_OFFTARGET_HYB_SCORE": 0.58,
    }

    result = run_oligotox(candidate, omics_payload=omics_payload)

    assert "omics" in result
    assert "adjustment_summary" in result["omics"]
    assert "recommended_assays" in result["omics"]
    assert "explanation" in result["omics"]

    print("=" * 80)
    print("OMICS DECISION LAYER DEMO")
    print("=" * 80)
    print("Dominant endpoint:", result["dominant_endpoint"])
    print()

    for endpoint, info in sorted(result["per_endpoint"].items()):
        base_prob = info.get("base_prob", info["prob"])
        print(
            f"{endpoint:28s} "
            f"base={base_prob:.3f} "
            f"final={info['prob']:.3f} "
            f"positive={info['predicted_positive']}"
        )

    print()
    print("Omics summary:")
    print(result["omics"]["explanation"]["summary"])
    print()
    print("Assay recommendations:")
    for rec in result["omics"]["recommended_assays"]:
        print(f"- [{rec['priority']}] {rec['assay']}: {rec['rationale']}")


if __name__ == "__main__":
    test_omics_decision_layer()
