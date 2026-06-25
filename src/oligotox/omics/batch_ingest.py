from __future__ import annotations

from pathlib import Path
import json
import pandas as pd

from oligotox.omics.batch_manifest import load_manifest, manifest_to_rows
from oligotox.omics.real_data_ingest import ingest_deg_file_to_readout, readout_to_payload


def _flatten_readout(candidate_id: str, study_id: str, comparison_label: str, readout, qc: dict) -> dict:
    row = {
        "candidate_id": candidate_id,
        "study_id": study_id,
        "comparison_label": comparison_label,
        "modality": readout.metadata.modality,
        "tissue": readout.metadata.tissue,
        "species": readout.metadata.species,
        "system": readout.metadata.system,
        "platform": readout.metadata.platform,
        "timepoint_hr": readout.metadata.timepoint_hr,
        "dose_level": readout.metadata.dose_level,
        "dose_value": readout.metadata.dose_value,
        "dose_unit": readout.metadata.dose_unit,
        "replicate_count": readout.metadata.replicate_count,
        "omics_confidence": readout.confidence.overall(),
        "omics_coverage_score": readout.confidence.coverage_score,
        "omics_quality_score": readout.confidence.quality_score,
        "omics_concordance_score": readout.confidence.concordance_score,
        "omics_translational_relevance": readout.confidence.translational_relevance,
        "qc_profile": qc.get("qc_profile"),
        "qc_passes_basic": qc.get("passes_basic_qc"),
        "qc_n_rows": qc.get("n_rows"),
        "qc_n_unique_genes": qc.get("n_unique_genes"),
        "qc_significant_gene_count": qc.get("significant_gene_count"),
    }

    for k, v in readout.pathway_scores.items():
        row[f"omics_pw_{k}"] = v

    for k, v in readout.sentinel_markers.items():
        row[f"omics_sentinel_{k}"] = v

    return row


def batch_ingest_manifest(
    manifest_path: str | Path,
    output_table_path: str | Path,
    output_qc_path: str | Path,
    payload_dir: str | Path,
    fail_fast: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    manifest_df = load_manifest(manifest_path)
    rows = manifest_to_rows(manifest_df)

    payload_dir = Path(payload_dir)
    payload_dir.mkdir(parents=True, exist_ok=True)

    good_rows = []
    qc_rows = []

    for row in rows:
        metadata = {
            "candidate_id": row.candidate_id,
            "study_id": row.study_id,
            "comparison_label": row.comparison_label,
            "tissue": row.tissue,
            "species": row.species,
            "system": row.system,
            "platform": row.platform,
            "timepoint_hr": row.timepoint_hr,
            "dose_level": row.dose_level,
            "dose_value": row.dose_value,
            "dose_unit": row.dose_unit,
            "replicate_count": row.replicate_count,
        }

        try:
            readout, qc = ingest_deg_file_to_readout(
                row.file_path,
                metadata=metadata,
                qc_profile=row.qc_profile,
                fail_on_qc=True,
            )

            payload = readout_to_payload(readout)
            payload_path = payload_dir / f"{row.candidate_id}__{row.study_id}.json"
            payload_path.write_text(json.dumps(payload, indent=2))

            flat = _flatten_readout(
                candidate_id=row.candidate_id,
                study_id=row.study_id,
                comparison_label=row.comparison_label,
                readout=readout,
                qc=qc,
            )
            flat["payload_json"] = str(payload_path)

            qc_rows.append({
                "candidate_id": row.candidate_id,
                "study_id": row.study_id,
                "file_path": row.file_path,
                "status": "ok",
                **qc,
            })
            good_rows.append(flat)

        except Exception as e:
            qc_rows.append({
                "candidate_id": row.candidate_id,
                "study_id": row.study_id,
                "file_path": row.file_path,
                "status": "failed",
                "error": str(e),
            })
            if fail_fast:
                raise

    out_df = pd.DataFrame(good_rows)
    qc_df = pd.DataFrame(qc_rows)

    output_table_path = Path(output_table_path)
    output_qc_path = Path(output_qc_path)
    output_table_path.parent.mkdir(parents=True, exist_ok=True)
    output_qc_path.parent.mkdir(parents=True, exist_ok=True)

    out_df.to_csv(output_table_path, index=False)
    qc_df.to_csv(output_qc_path, index=False)

    return out_df, qc_df
