from __future__ import annotations

from oligotox.omics.mechanistic_schema import OmicsMetadata, OmicsConfidence, OmicsReadout
from oligotox.omics.feature_dictionary import PATHWAY_FEATURES, SENTINEL_MARKERS


def _num(v, default=0.0) -> float:
    try:
        if v is None:
            return float(default)
        return float(v)
    except Exception:
        return float(default)


def parse_omics_payload(payload: dict) -> OmicsReadout:
    metadata = OmicsMetadata(
        candidate_id=payload.get("candidate_id"),
        modality=payload.get("modality"),
        tissue=payload.get("tissue"),
        species=payload.get("species"),
        system=payload.get("system"),
        platform=payload.get("platform"),
        timepoint_hr=_num(payload.get("timepoint_hr"), 0.0),
        dose_level=payload.get("dose_level"),
        dose_value=_num(payload.get("dose_value"), 0.0),
        dose_unit=payload.get("dose_unit"),
        replicate_count=int(_num(payload.get("replicate_count"), 0)),
    )

    confidence = OmicsConfidence(
        coverage_score=_num(payload.get("coverage_score"), 0.0),
        quality_score=_num(payload.get("quality_score"), 0.0),
        concordance_score=_num(payload.get("concordance_score"), 0.0),
        translational_relevance=_num(payload.get("translational_relevance"), 0.0),
    )

    pathway_scores = {}
    sentinel_markers = {}

    for key in PATHWAY_FEATURES:
        pathway_scores[key] = _num(payload.get(f"pw_{key}", payload.get(key, 0.0)), 0.0)

    for key in SENTINEL_MARKERS:
        sentinel_markers[key] = _num(payload.get(f"sentinel_{key}", payload.get(key, 0.0)), 0.0)

    return OmicsReadout(
        metadata=metadata,
        pathway_scores=pathway_scores,
        sentinel_markers=sentinel_markers,
        confidence=confidence,
        raw_summary={k: v for k, v in payload.items() if k not in pathway_scores and k not in sentinel_markers},
    )
