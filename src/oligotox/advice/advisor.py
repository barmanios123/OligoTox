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


# --- Full written summary -------------------------------------------------
#
# `explain` returns the numbers. This turns them into something a person can
# read without knowing the codebase, and it reports the sequence's measured
# parameters next to the predictions so the two are never confused. Every
# parameter here is computed from the sequence by the same functions training
# used. Nothing is estimated and nothing is rounded into a claim.

# What each parameter is, and why it is in the model at all. Shown with the
# number so a reader is not left guessing what "dist_3p_to_g" means.
PARAMETER_NOTES: dict[str, tuple[str, str]] = {
    "seq_length":     ("length", "nucleotides"),
    "gc_content":     ("GC content", "fraction of G or C"),
    "seed_gc":        ("seed GC", "positions 2-8, which drive off-target load"),
    "cpg_density":    ("CpG density", "CpG dinucleotides per position"),
    "tlr7_score":     ("TLR7 score", "innate immune motif content"),
    "g_count":        ("G count", "Hagedorn 2022 links G content to acute neurotoxicity"),
    "g_fraction":     ("G fraction", "G nucleotides per position"),
    "dist_3p_to_g":   ("3' distance to nearest G", "nucleotides from the 3' end"),
    "longest_g_run":  ("longest G run", "consecutive Gs"),
    "struct_mfe_proxy": ("folding energy", "kcal/mol, real ViennaRNA MFE"),
    "struct_hairpin_score": ("hairpin score", "self-complementarity"),
    "struct_palindrome_fraction": ("palindrome fraction",
                                   "4-mers whose reverse complement is present"),
    "struct_longest_homopolymer": ("longest homopolymer",
                                   "consecutive identical bases"),
    "struct_dinucleotide_repeats": ("dinucleotide repeats", "repeated 2-mers"),
    "ps_fraction":    ("phosphorothioate fraction", "backbone linkages modified"),
    "lna_fraction":   ("LNA fraction", "sugar-modified positions"),
}


def _format_probability(value: float) -> str:
    """Small probabilities must not be rounded into zero.

    An offtarget score of 8.5e-06 printed as "0.000" reads as impossible
    rather than as very unlikely, which is a different claim.
    """
    return f"{value:.3f}" if value >= 0.001 else f"{value:.2e}"


def _describe(features: dict[str, float]) -> tuple[list[dict[str, Any]], list[str]]:
    """The parameters worth showing, plus any this function asked for and did
    not get.

    The missing list matters more than it looks. Three parameters were named
    here with the wrong keys and simply never appeared in the output, which is
    indistinguishable from a sequence that has no hairpin. A name that stops
    matching the feature pipeline has to be visible, not silently dropped.
    """
    out, missing = [], []
    for key, (label, meaning) in PARAMETER_NOTES.items():
        if key in features:
            out.append({
                "key": key, "label": label, "meaning": meaning,
                "value": round(float(features[key]), 4),
            })
        else:
            missing.append(key)
    return out, missing


def summarize(
    sequence: str,
    modality: str = "ASO",
    chemistry: str = "PS",
    ps_fraction: float = 1.0,
    model_dir: str | Path = "artifacts/models",
    config_dir: str | Path = "artifacts/config",
) -> dict[str, Any]:
    """Everything known about one sequence, with a readable write-up.

    Three parts, kept apart on purpose:

      parameters   measured properties of the sequence itself. Descriptive,
                   not predictive, and true regardless of any model.
      predictions  only endpoints tiered `validated`. The other ten models
                   run but are not shown, because a number the evidence does
                   not support is worse than no number.
      text         the same thing in prose, including what is missing.
    """
    model_dir, config_dir = Path(model_dir), Path(config_dir)
    seq = sequence.strip().upper().replace("U", "T")
    features = candidate_features(seq, modality, chemistry, ps_fraction=ps_fraction)
    tox, potency = _score(seq, modality, chemistry, ps_fraction, model_dir, config_dir)
    parameters, missing = _describe(features)

    lines = [
        f"Sequence: {seq}",
        f"{len(seq)} nt, {modality}, {chemistry} backbone at "
        f"{ps_fraction:.0%} phosphorothioate.",
        "",
        "Measured parameters of the sequence",
    ]
    for item in parameters:
        lines.append(f"  {item['label']}: {item['value']} ({item['meaning']})")

    lines += ["", "Predictions the evidence supports"]
    if tox:
        for name, prob in sorted(tox.items()):
            lines.append(f"  {name}: {_format_probability(prob)}")
        lines.append("  These are the only endpoints tiered validated. Ten other "
                     "toxicity models run and are deliberately not shown.")
    else:
        lines.append("  None. No toxicity endpoint reached the validated tier for "
                     "this input, so no toxicity prediction is presented.")

    lines += ["", "Per-assay potency"]
    if potency:
        for name, score in sorted(potency.items()):
            lines.append(f"  {name}: {score:.3f}")
        lines.append("  Scored per assay rather than pooled, because pooling the "
                     "assays made the prediction worse.")
    else:
        lines.append("  None available for this modality.")

    if missing:
        lines += ["", "Parameters this summary expected and did not find",
                  f"  {', '.join(missing)}",
                  "  These names no longer match the feature pipeline, so the "
                  "summary is incomplete rather than wrong."]

    lines += [
        "",
        "What this does not tell you",
        "  The parameters above are properties of the oligo alone. Target-site "
        "accessibility and off-target burden both need the target transcript, "
        "which is not supplied here.",
        "  A probability is only as good as its calibration. Endpoints that rank "
        "but do not calibrate are labelled as such rather than quoted as risk.",
    ]

    return {
        "sequence": seq,
        "modality": modality,
        "chemistry": chemistry,
        "ps_fraction": ps_fraction,
        "parameters": parameters,
        "validated_toxicity": tox,
        "per_assay_potency": potency,
        "parameters_unavailable": missing,
        "text": "\n".join(lines),
    }
