import os
import sys
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from oligotox.omics.batch_ingest import batch_ingest_manifest
from oligotox.omics.training_join import build_candidate_level_omics_table


def test_batch_manifest_ingest():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    manifest = os.path.join(root, "manifests", "omics_study_manifest_demo.csv")
    out_table = os.path.join(root, "data", "processed", "omics", "omics_training_table.csv")
    out_qc = os.path.join(root, "data", "processed", "omics", "omics_ingest_qc.csv")
    payload_dir = os.path.join(root, "output", "omics_payloads")

    df, qc = batch_ingest_manifest(
        manifest_path=manifest,
        output_table_path=out_table,
        output_qc_path=out_qc,
        payload_dir=payload_dir,
        fail_fast=True,
    )

    assert len(df) == 2
    assert len(qc) == 2
    assert "candidate_id" in df.columns
    assert "omics_pw_interferon" in df.columns or "omics_pw_renal_tubular_stress" in df.columns

    cand = build_candidate_level_omics_table(
        omics_training_table_path=out_table,
        output_path=os.path.join(root, "data", "processed", "omics", "candidate_level_omics_features.csv"),
    )
    assert len(cand) == 2

    print("=" * 80)
    print("BATCH MANIFEST INGEST TEST")
    print("=" * 80)
    print(df.head().to_string(index=False))
    print()
    print(qc.to_string(index=False))


if __name__ == "__main__":
    test_batch_manifest_ingest()
