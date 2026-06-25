from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Any

from oligotox.omics.mechanistic_schema import OmicsReadout


@dataclass
class EndpointAdjustment:
    endpoint: str
    delta: float
    reasons: list[str]


def _pw(readout: OmicsReadout, key: str) -> float:
    return float(readout.pathway_scores.get(key, 0.0))


def _sent(readout: OmicsReadout, key: str) -> float:
    return float(readout.sentinel_markers.get(key, 0.0))


def compute_endpoint_adjustments(readout: OmicsReadout) -> dict[str, EndpointAdjustment]:
    confidence = readout.confidence.overall()

    out: dict[str, EndpointAdjustment] = {}

    def add(endpoint: str, delta: float, reason: str) -> None:
        scaled = float(delta) * confidence
        if endpoint not in out:
            out[endpoint] = EndpointAdjustment(endpoint=endpoint, delta=0.0, reasons=[])
        out[endpoint].delta += scaled
        out[endpoint].reasons.append(reason)

    if _pw(readout, "interferon") > 0.5 or _pw(readout, "tlr_innate") > 0.5:
        add("innate_immune_activation", 0.18, "High interferon / innate pathway activity")
        add("cytokine_release", 0.10, "Innate signaling suggests cytokine liability")

    if _pw(readout, "nfkb") > 0.5 or _pw(readout, "cytokine_signaling") > 0.5:
        add("cytokine_release", 0.14, "NF-kB / cytokine pathways elevated")
        add("kupffer_cell_activation", 0.08, "Inflammatory signaling may activate liver macrophage programs")

    if _pw(readout, "complement") > 0.5 or _sent(readout, "C3") > 0.5 or _sent(readout, "C5") > 0.5:
        add("complement_activation", 0.20, "Complement program elevated")
        add("thrombocytopenia", 0.06, "Complement-linked platelet liability plausible")
        add("coagulation_effect", 0.06, "Complement/coagulation cross-talk signal")

    if _pw(readout, "platelet") > 0.5 or _sent(readout, "PF4") > 0.5:
        add("thrombocytopenia", 0.18, "Platelet biology signature elevated")
        add("coagulation_effect", 0.10, "Platelet / coagulation biology perturbed")

    if _pw(readout, "hepatocyte_injury") > 0.5 or _pw(readout, "apoptosis") > 0.5:
        add("hepatocellular_injury", 0.18, "Hepatocyte injury / apoptosis program elevated")

    if _pw(readout, "cholestasis_bile_acid") > 0.5:
        add("hepatic_vacuolation", 0.08, "Bile acid / hepatobiliary stress pattern")

    if _pw(readout, "renal_tubular_stress") > 0.5:
        add("renal_tubular_injury", 0.18, "Renal tubular stress program elevated")
        add("renal_accumulation", 0.10, "Renal liability / accumulation pattern")

    if _pw(readout, "rnaseh1_long_transcript_risk") > 0.5 or _sent(readout, "RNASEH1_DEP_LONGTX") > 0.5:
        add("hepatocellular_injury", 0.12, "RNase H1 / long-transcript liability signal")
        add("offtarget_gene_silencing", 0.18, "Long-transcript depletion risk elevated")

    if _pw(readout, "hybridization_offtarget") > 0.5 or _sent(readout, "OFFTARGET_HYB_SCORE") > 0.5:
        add("offtarget_gene_silencing", 0.22, "Hybridization-dependent off-target signature")

    if _pw(readout, "kupffer_activation") > 0.5:
        add("kupffer_cell_activation", 0.18, "Kupffer / macrophage activation signature")

    if _pw(readout, "oxidative_stress") > 0.5 or _pw(readout, "mitochondrial_stress") > 0.5:
        add("hepatocellular_injury", 0.10, "Oxidative / mitochondrial stress elevation")
        add("renal_tubular_injury", 0.05, "General cellular stress may affect kidney")

    return out


def apply_adjustments(
    base_probs: Dict[str, float],
    readout: OmicsReadout | None,
) -> tuple[Dict[str, float], Dict[str, Any]]:
    if readout is None:
        return base_probs, {"used_omics": False, "adjustments": {}}

    adjustments = compute_endpoint_adjustments(readout)
    adjusted = dict(base_probs)

    for endpoint, adj in adjustments.items():
        base = adjusted.get(endpoint, 0.0)
        adjusted[endpoint] = max(0.0, min(1.0, base + adj.delta))

    explanation = {
        "used_omics": True,
        "omics_confidence": readout.confidence.overall(),
        "adjustments": {
            endpoint: {
                "delta": adj.delta,
                "reasons": adj.reasons,
            }
            for endpoint, adj in adjustments.items()
        },
    }
    return adjusted, explanation
