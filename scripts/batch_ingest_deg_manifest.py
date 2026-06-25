from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from oligotox.omics.batch_ingest import batch_ingest_manifest


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--output-table", default="data/processed/omics/omics_training_table.csv")
    ap.add_argument("--output-qc", default="data/processed/omics/omics_ingest_qc.csv")
    ap.add_argument("--payload-dir", default="output/omics_payloads")
    ap.add_argument("--fail-fast", action="store_true")
    args = ap.parse_args()

    out_df, qc_df = batch_ingest_manifest(
        manifest_path=args.manifest,
        output_table_path=args.output_table,
        output_qc_path=args.output_qc,
        payload_dir=args.payload_dir,
        fail_fast=args.fail_fast,
    )

    print("=" * 80)
    print("BATCH OMICS INGEST COMPLETE")
    print("=" * 80)
    print(f"Successful rows: {len(out_df)}")
    print(f"QC rows total  : {len(qc_df)}")
    if not qc_df.empty:
        print()
        print(qc_df.to_string(index=False))


if __name__ == "__main__":
    main()
