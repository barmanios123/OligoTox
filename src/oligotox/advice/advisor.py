"""In-silico single-substitution scan over the gated models.

For a candidate sequence, mutate every position to every other base, score
each variant through `score_candidate` and the per-assay potency models, and
report which edits move which numbers.

The output is deliberately two-sided. An edit that lowers predicted CNS risk
often lowers predicted knockdown too, and a recommendation that hides the
trade-off is not advice, it is marketing. Every suggested edit carries both
effects.

Cost: a 20-mer gives 60 variants, each a full scoring pass. That is a second
or two, which is why the scan is exhaustive rather than greedy.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from oligotox.efficacy.potency import score_efficacy
from oligotox.features.candidate import candidate_features
from oligotox.pipelines.score_candidate import score_candidate

BASES = ("A", "C", "G", "T")
# Below this, a change is noise rather than a finding worth acting on.
MIN_MEANINGFUL_DELTA = 0.01


@dataclass
class Edit:
    position: int           # 1-based, as a bench scientist counts
    from_base: str
    to_base: str
    sequence: str
    toxicity_change: dict[str, float] = field(default_factory=dict)
    potency_change: dict[str, float] = field(default_factory=dict)

    @property
    def total_toxicity_reduction(self) -> float:
        """Positive means the edit lowers risk."""
        return -sum(self.toxicity_change.values())

    @property
    def mean_potency_change(self) -> float:
        if not self.potency_change:
            return 0.0
        return sum(self.potency_change.values()) / len(self.potency_change)

    def describe(self) -> str:
        tox = ", ".join(
            f"{k.replace('_', ' ')} {v:+.3f}" for k, v in sorted(self.toxicity_change.items())
        ) or "no validated endpoint moved"
        line = f"{self.from_base}{self.position}{self.to_base}: {tox}"
        if self.potency_change:
            line += f"; predicted knockdown {self.mean_potency_change:+.1f}%"
        return line


@dataclass
class SequenceAdvice:
    sequence: str
    modality: str
    validated_endpoints: list[str]
    baseline_toxicity: dict[str, float]
    baseline_potency: dict[str, float]
    edits: list[Edit]
    notes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "sequence": self.sequence,
            "modality": self.modality,
            "validated_endpoints": self.validated_endpoints,
            "baseline_toxicity": self.baseline_toxicity,
            "baseline_potency": self.baseline_potency,
            "notes": self.notes,
            "edits": [
                {
                    "edit": f"{e.from_base}{e.position}{e.to_base}",
                    "sequence": e.sequence,
                    "toxicity_change": e.toxicity_change,
                    "potency_change": e.potency_change,
                    "toxicity_reduction": round(e.total_toxicity_reduction, 4),
                    "mean_potency_change": round(e.mean_potency_change, 2),
                }
                for e in self.edits
            ],
        }


def _score(sequence: str, modality: str, chemistry: str, ps_fraction: float,
           model_dir: Path, config_dir: Path) -> tuple[dict, dict]:
    feats = candidate_features(sequence, modality, chemistry, ps_fraction=ps_fraction)
    tox = score_candidate(
        candidate_features=feats,
        processed_dir=model_dir.parent.parent / "data" / "processed",
        model_dir=model_dir,
        config_dir=config_dir,
    )
    validated = {
        name: info["prob"]
        for name, info in tox["per_endpoint"].items()
        if info.get("evidence_tier") == "validated"
    }
    efficacy = score_efficacy(sequence, modality, model_dir, config_dir)
    potency = {
        name: pred["score"]
        for name, pred in efficacy.get("predictions", {}).items()
        if pred.get("tier") == "context"
    }
    return validated, potency


def explain(
    sequence: str,
    modality: str = "ASO",
    chemistry: str = "PS",
    ps_fraction: float = 1.0,
    model_dir: str | Path = "artifacts/models",
    config_dir: str | Path = "artifacts/config",
) -> dict[str, Any]:
    """Current predictions, restricted to what the evidence supports."""
    model_dir, config_dir = Path(model_dir), Path(config_dir)
    tox, potency = _score(sequence, modality, chemistry, ps_fraction, model_dir, config_dir)
    return {
        "sequence": sequence.strip().upper().replace("U", "T"),
        "validated_toxicity": tox,
        "per_assay_potency": potency,
        "note": "Only endpoints tiered validated are shown. The other ten toxicity "
                "models run but are not presented as predictions.",
    }


def suggest(
    sequence: str,
    modality: str = "ASO",
    chemistry: str = "PS",
    ps_fraction: float = 1.0,
    model_dir: str | Path = "artifacts/models",
    config_dir: str | Path = "artifacts/config",
    top_n: int = 10,
) -> SequenceAdvice:
    """Every single-base substitution, scored and ranked by risk reduction."""
    model_dir, config_dir = Path(model_dir), Path(config_dir)
    clean = sequence.strip().upper().replace("U", "T")
    if not clean or set(clean) - set(BASES):
        raise ValueError(f"sequence must be A/C/G/T(/U), got {sequence!r}")

    base_tox, base_potency = _score(clean, modality, chemistry, ps_fraction,
                                    model_dir, config_dir)
    notes = []
    if not base_tox:
        notes.append(
            "No validated toxicity endpoint applies here, so no edit can be "
            "recommended on toxicity grounds. Potency effects are still shown."
        )

    edits: list[Edit] = []
    for index, original in enumerate(clean):
        for replacement in BASES:
            if replacement == original:
                continue
            variant = clean[:index] + replacement + clean[index + 1:]
            tox, potency = _score(variant, modality, chemistry, ps_fraction,
                                  model_dir, config_dir)
            edit = Edit(
                position=index + 1, from_base=original, to_base=replacement,
                sequence=variant,
                toxicity_change={k: round(tox[k] - base_tox[k], 5)
                                 for k in base_tox if abs(tox[k] - base_tox[k]) >= MIN_MEANINGFUL_DELTA},
                potency_change={k: round(potency[k] - base_potency[k], 3)
                                for k in base_potency if abs(potency[k] - base_potency[k]) >= 0.5},
            )
            if edit.toxicity_change or edit.potency_change:
                edits.append(edit)

    edits.sort(key=lambda e: (-e.total_toxicity_reduction, -e.mean_potency_change))

    if base_tox and not any(e.total_toxicity_reduction > 0 for e in edits):
        notes.append(
            "No single substitution lowers a validated toxicity score by more than "
            f"{MIN_MEANINGFUL_DELTA}. The sequence may already be near the model's floor."
        )
    notes.append(
        "Each number is the difference between two real predictions, not a "
        "coefficient. Predictions carry the limits of the model that made them: "
        "see the evidence string for each endpoint."
    )

    return SequenceAdvice(
        sequence=clean, modality=modality,
        validated_endpoints=sorted(base_tox),
        baseline_toxicity=base_tox, baseline_potency=base_potency,
        edits=edits[:top_n], notes=notes,
    )
