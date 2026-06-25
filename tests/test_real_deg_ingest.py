import os
import sys
import json

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from oligotox.omics.real_data_ingest import ingest_deg_file_to_readout, readout_to_payload
from api import run_oligotox


def test_real_style_deg_ingestion():
    path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "omics",
        "examples",
        "example_real_deg.csv",
    )

    readout, qc = ingest_deg_file_to_readout(
        path,
        metadata={
            "candidate_id": "real_deg_demo_001",
            "study_id": "demo_study",
            "comparison_label": "treated_vs_control",
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
        qc_profile="relaxed",
        fail_on_qc=True,
    )

    assert qc["passes_basic_qc"] is True
    payload = readout_to_payload(readout)

    candidate = {
        "sequence": "GCUAUGCUAGCUAGCUAGCUAGCUAGC",
        "is_lnp": True,
        "is_iv": True,
        "dose_mg_kg": 1.5,
        "repeat_dose": False,
    }

    result = run_oligotox(candidate, omics_payload=payload)

    assert "omics" in result
    assert result["omics"]["adjustment_summary"]["used_omics"] is True

    print("=" * 80)
    print("REAL-STYLE DEG INGEST DEMO")
    print("=" * 80)
    print("QC:", json.dumps(qc, indent=2))
    print("Omics summary:", result["omics"]["explanation"]["summary"])


if __name__ == "__main__":
    test_real_style_deg_ingestion()
