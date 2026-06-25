from __future__ import annotations

REQUIRED_METADATA_FIELDS = [
    "candidate_id",
    "study_id",
    "comparison_label",
    "tissue",
    "species",
    "system",
    "platform",
    "timepoint_hr",
    "dose_level",
    "dose_value",
    "dose_unit",
]

OPTIONAL_METADATA_FIELDS = [
    "replicate_count",
    "lab",
    "source_file",
    "normalization_method",
    "dge_method",
    "reference_group",
    "treatment_group",
]
