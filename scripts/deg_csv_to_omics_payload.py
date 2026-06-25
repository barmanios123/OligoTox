from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from omics import ingest_deg_file_to_readout, readout_to_payload


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="Path to DEG CSV/TSV")
    ap.add_argument("--output", required=True, help="Output JSON payload path")
    ap.add_argument("--candidate-id", required=True)
    ap.add_argument("--study-id", default="unspecified_study")
    ap.add_argument("--comparison-label", default="treatment_vs_control")
    ap.add_argument("--tissue", default="liver")
    ap.add_argument("--species", default="human")
    ap.add_argument("--system", default="primary_hepatocyte")
    ap.add_argument("--platform", default="rna_seq")
    ap.add_argument("--timepoint-hr", type=float, default=24.0)
    ap.add_argument("--dose-level", default="mid")
    ap.add_argument("--dose-value", type=float, default=0.0)
    ap.add_argument("--dose-unit", default="uM")
    ap.add_argument("--replicate-count", type=int, default=3)
    ap.add_argument("--qc-profile", choices=["strict", "standard", "relaxed"], default="standard")
    ap.add_argument("--allow-qc-fail", action="store_true")
    args = ap.parse_args()

    metadata = {
        "candidate_id": args.candidate_id,
        "study_id": args.study_id,
        "comparison_label": args.comparison_label,
        "tissue": args.tissue,
        "species": args.species,
        "system": args.system,
        "platform": args.platform,
        "timepoint_hr": args.timepoint_hr,
        "dose_level": args.dose_level,
        "dose_value": args.dose_value,
        "dose_unit": args.dose_unit,
        "replicate_count": args.replicate_count,
    }

    readout, qc = ingest_deg_file_to_readout(
        args.input,
        metadata=metadata,
        qc_profile=args.qc_profile,
        fail_on_qc=not args.allow_qc_fail,
    )
    payload = readout_to_payload(readout)

    out = {
        "metadata": metadata,
        "qc": qc,
        "omics_payload": payload,
        "omics_summary": {
            "top_pathways": sorted(readout.pathway_scores.items(), key=lambda x: x[1], reverse=True)[:8],
            "top_sentinels": sorted(readout.sentinel_markers.items(), key=lambda x: x[1], reverse=True)[:8],
            "confidence_overall": readout.confidence.overall(),
        },
    }

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(out, indent=2))
    print(f"Wrote: {output_path}")
    print(json.dumps(qc, indent=2))


if __name__ == "__main__":
    main()
