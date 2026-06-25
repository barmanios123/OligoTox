from __future__ import annotations

from typing import Any

from oligotox.omics.mechanistic_schema import OmicsReadout


def recommend_assays(readout: OmicsReadout | None, endpoint_probs: dict[str, float]) -> list[dict[str, Any]]:
    recs: list[dict[str, Any]] = []

    def add(name: str, rationale: str, priority: str = "medium") -> None:
        recs.append({"assay": name, "priority": priority, "rationale": rationale})

    if endpoint_probs.get("innate_immune_activation", 0.0) > 0.35 or endpoint_probs.get("cytokine_release", 0.0) > 0.35:
        add("PBMC cytokine release panel", "Innate / cytokine liability appears elevated", "high")
        add("TLR / IFN pathway reporter assay", "Clarify innate immune mechanism", "high")

    if endpoint_probs.get("complement_activation", 0.0) > 0.30:
        add("Human complement activation panel (C3a/C5a/sC5b-9)", "Complement liability signal elevated", "high")

    if endpoint_probs.get("thrombocytopenia", 0.0) > 0.25 or endpoint_probs.get("coagulation_effect", 0.0) > 0.25:
        add("Platelet activation / aggregation assay", "Platelet or coagulation risk is non-trivial", "high")
        add("PF4 / coagulation biomarker panel", "Resolve thrombocytopenia vs coagulation liability", "medium")

    if endpoint_probs.get("hepatocellular_injury", 0.0) > 0.35:
        add("Primary hepatocyte tox time-course", "Hepatocellular injury risk elevated", "high")
        add("Liver transcriptomics pathway panel", "Mechanistically resolve liver injury program", "high")

    if endpoint_probs.get("renal_tubular_injury", 0.0) > 0.30:
        add("Renal proximal tubule assay", "Renal tubular liability signal elevated", "high")

    if endpoint_probs.get("offtarget_gene_silencing", 0.0) > 0.30:
        add("Transcriptome-wide off-target profiling", "Potential hybridization-dependent off-target risk", "high")

    if readout is not None and readout.confidence.overall() < 0.4:
        add("Repeat omics with better coverage / replicates", "Current omics evidence confidence is limited", "medium")

    return recs
