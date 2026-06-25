from __future__ import annotations

from oligotox.omics.mechanistic_schema import OmicsReadout


def top_items(d: dict[str, float], k: int = 5) -> list[tuple[str, float]]:
    return sorted(d.items(), key=lambda x: x[1], reverse=True)[:k]


def build_omics_explanation(readout: OmicsReadout | None) -> dict:
    if readout is None:
        return {
            "summary": "No omics evidence provided.",
            "top_pathways": [],
            "top_sentinels": [],
        }

    top_pathways = top_items(readout.pathway_scores, k=6)
    top_sentinels = top_items(readout.sentinel_markers, k=6)

    summary_parts = []
    if top_pathways:
        summary_parts.append(
            "Top pathways: " + ", ".join([f"{k}={v:.2f}" for k, v in top_pathways[:4]])
        )
    if top_sentinels:
        summary_parts.append(
            "Top sentinels: " + ", ".join([f"{k}={v:.2f}" for k, v in top_sentinels[:4]])
        )

    return {
        "summary": " | ".join(summary_parts) if summary_parts else "Omics evidence present.",
        "top_pathways": top_pathways,
        "top_sentinels": top_sentinels,
        "metadata": {
            "tissue": readout.metadata.tissue,
            "species": readout.metadata.species,
            "system": readout.metadata.system,
            "timepoint_hr": readout.metadata.timepoint_hr,
            "platform": readout.metadata.platform,
        },
        "confidence": {
            "coverage_score": readout.confidence.coverage_score,
            "quality_score": readout.confidence.quality_score,
            "concordance_score": readout.confidence.concordance_score,
            "translational_relevance": readout.confidence.translational_relevance,
            "overall": readout.confidence.overall(),
        },
    }
