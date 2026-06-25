import os
import sys
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from oligotox.omics.deg_to_readout import deg_table_to_readout
from api import run_oligotox


def test_deg_to_omics_readout_and_decision_layer():
    deg = pd.DataFrame(
        {
            "gene_symbol": [
                "ISG15",
                "OAS1",
                "MX1",
                "IFI44",
                "IL6",
                "CXCL10",
                "RNASEH1",
                "FOXP1",
                "DMD",
                "TTN",
            ],
            "log2fc": [2.1, 1.8, 1.7, 1.6, 1.2, 2.0, 1.3, 1.1, 1.0, 1.4],
            "pval_adj": [1e-6, 2e-5, 3e-5, 4e-5, 1e-3, 2e-4, 1e-2, 5e-3, 0.02, 0.01],
        }
    )

    readout = deg_table_to_readout(
        deg,
        gene_col="gene_symbol",
        effect_col="log2fc",
        pval_col="pval_adj",
        metadata={
            "candidate_id": "deg_demo_001",
            "tissue": "liver",
            "species": "human",
            "system": "primary_hepatocyte",
            "platform": "rna_seq",
            "timepoint_hr": 24,
            "dose_level": "high",
            "dose_value": 10.0,
            "dose_unit": "uM",
            "replicate_count": 3,
        },
    )

    candidate = {
        "sequence": "GCUAUGCUAGCUAGCUAGCUAGCUAGC",
        "is_lnp": True,
        "is_iv": True,
        "dose_mg_kg": 1.5,
        "repeat_dose": False,
    }

    omics_payload = {
        "candidate_id": readout.metadata.candidate_id,
        "tissue": readout.metadata.tissue,
        "species": readout.metadata.species,
        "system": readout.metadata.system,
        "platform": readout.metadata.platform,
        "timepoint_hr": readout.metadata.timepoint_hr,
        "dose_level": readout.metadata.dose_level,
        "dose_value": readout.metadata.dose_value,
        "dose_unit": readout.metadata.dose_unit,
        "replicate_count": readout.metadata.replicate_count,
        "coverage_score": readout.confidence.coverage_score,
        "quality_score": readout.confidence.quality_score,
        "concordance_score": readout.confidence.concordance_score,
        "translational_relevance": readout.confidence.translational_relevance,
    }

    for k, v in readout.pathway_scores.items():
        omics_payload[f"pw_{k}"] = v
    for k, v in readout.sentinel_markers.items():
        omics_payload[f"sentinel_{k}"] = v

    result = run_oligotox(candidate, omics_payload=omics_payload)

    assert "omics" in result
    assert result["omics"]["adjustment_summary"]["used_omics"] is True
    assert result["per_endpoint"]["innate_immune_activation"]["prob"] >= result["per_endpoint"]["innate_immune_activation"]["base_prob"]
    assert result["per_endpoint"]["offtarget_gene_silencing"]["prob"] >= result["per_endpoint"]["offtarget_gene_silencing"]["base_prob"]

    print("=" * 80)
    print("DEG → OmicsReadout → OligoTox demo")
    print("=" * 80)
    print("Omics summary:", result["omics"]["explanation"]["summary"])
    print()
    for endpoint, info in sorted(result["per_endpoint"].items()):
        print(
            f"{endpoint:28s} "
            f"base={info.get('base_prob', info['prob']):.3f} "
            f"final={info['prob']:.3f} "
            f"positive={info['predicted_positive']}"
        )


if __name__ == "__main__":
    test_deg_to_omics_readout_and_decision_layer()
