from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


@dataclass
class OmicsMetadata:
    candidate_id: str | None = None
    modality: str | None = None
    tissue: str | None = None
    species: str | None = None
    system: str | None = None
    platform: str | None = None
    timepoint_hr: float | None = None
    dose_level: str | None = None
    dose_value: float | None = None
    dose_unit: str | None = None
    replicate_count: int | None = None


@dataclass
class OmicsConfidence:
    coverage_score: float = 0.0
    quality_score: float = 0.0
    concordance_score: float = 0.0
    translational_relevance: float = 0.0

    def overall(self) -> float:
        vals = [
            float(self.coverage_score),
            float(self.quality_score),
            float(self.concordance_score),
            float(self.translational_relevance),
        ]
        return sum(vals) / len(vals)


@dataclass
class OmicsReadout:
    metadata: OmicsMetadata
    pathway_scores: dict[str, float] = field(default_factory=dict)
    sentinel_markers: dict[str, float] = field(default_factory=dict)
    confidence: OmicsConfidence = field(default_factory=OmicsConfidence)
    raw_summary: dict[str, Any] = field(default_factory=dict)

    def to_feature_dict(self) -> dict[str, float]:
        out = {}
        for k, v in self.pathway_scores.items():
            out[f"omics_pw_{k}"] = float(v)
        for k, v in self.sentinel_markers.items():
            out[f"omics_sentinel_{k}"] = float(v)
        out["omics_confidence"] = float(self.confidence.overall())
        out["omics_coverage_score"] = float(self.confidence.coverage_score)
        out["omics_quality_score"] = float(self.confidence.quality_score)
        out["omics_concordance_score"] = float(self.confidence.concordance_score)
        out["omics_translational_relevance"] = float(self.confidence.translational_relevance)
        out["omics_timepoint_hr"] = float(self.metadata.timepoint_hr or 0.0)
        out["omics_replicate_count"] = float(self.metadata.replicate_count or 0.0)
        return out
