"""Sequence-based potency scoring for the modality in hand.

The Efficacy tab was removed once because it held untrained model slots whose
outputs were not predictions. These models are trained and validated: see
reports/results/efficacy.json and the evidence string carried in
artifacts/config/efficacy_models.json, which the app displays next to every
number it shows.

Only a model whose modality matches the candidate is used. Scoring an ASO
with the siRNA potency model would be meaningless, so it returns nothing
rather than a number.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict

import joblib
import pandas as pd

from oligotox.features.candidate import g_motif_features, sequence_features


def guide_features(sequence: str) -> Dict[str, float]:
    """Design features of the guide strand.

    Must stay identical to scripts/evaluate_efficacy.guide_features; the
    tests assert the two agree.
    """
    seq = sequence.strip().upper().replace("U", "T")
    feats = sequence_features(seq)
    feats.pop("ps_x_len", None)
    feats.pop("lna_x_gc", None)
    feats.update(g_motif_features(seq))

    n = len(seq)
    five, three = seq[:5], seq[-5:]
    feats["five_prime_au"] = sum(c in "AT" for c in five) / max(len(five), 1)
    feats["three_prime_au"] = sum(c in "AT" for c in three) / max(len(three), 1)
    feats["asymmetry"] = feats["five_prime_au"] - feats["three_prime_au"]
    for base in "ACGT":
        feats[f"frac_{base}"] = seq.count(base) / n if n else 0.0
    # Position-specific identity. The literature's strongest siRNA potency
    # signals are positional, not compositional: the 5' end of the guide
    # (position 1), the cleavage-site region near position 10, and the 3' end.
    for pos, label in ((0, "p1"), (9, "p10"), (n - 1 if n else 0, "p_last")):
        base_at = seq[pos] if 0 <= pos < n else ""
        for base in "ACGT":
            feats[f"{label}_{base}"] = float(base_at == base)
    # Seed region AU content (positions 2-8), which drives off-target load and
    # duplex-end stability.
    seed = seq[1:8]
    feats["seed_au"] = sum(c in "AT" for c in seed) / max(len(seed), 1)
    # Duplex-end stability difference, the classic asymmetry rule in nucleotides.
    feats["end_stability_diff"] = (
        sum(c in "GC" for c in seq[-2:]) - sum(c in "GC" for c in seq[:2])
    )
    return feats


def load_efficacy_config(config_dir: str | Path) -> dict:
    path = Path(config_dir) / "efficacy_models.json"
    return json.loads(path.read_text()) if path.exists() else {}


def score_efficacy(
    sequence: str,
    modality: str,
    model_dir: str | Path,
    config_dir: str | Path,
) -> Dict[str, Any]:
    """Potency for the matching modality, or an explanation of why not.

    The returned score is the model's raw readout, in the units the training
    data used, and those units are returned alongside it. They are not a
    probability and not comparable across modalities.
    """
    config = load_efficacy_config(config_dir)
    model_dir = Path(model_dir)

    matched = {k: v for k, v in config.items()
               if v.get("modality", "").lower() == str(modality).lower()}
    if not matched:
        return {
            "available": False,
            "reason": f"No efficacy model is validated for modality {modality!r}. "
                      f"Models exist for: {sorted({v['modality'] for v in config.values()})}.",
            "predictions": {},
        }

    feats = guide_features(sequence)
    predictions = {}
    for name, spec in matched.items():
        path = model_dir / spec["model_file"]
        if not path.exists():
            continue
        cols = spec["feature_columns"]
        missing = [c for c in cols if c not in feats]
        if missing:
            raise KeyError(f"{name} needs features not computed here: {missing}")
        row = pd.DataFrame([[float(feats[c]) for c in cols]], columns=cols)
        predictions[name] = {
            "score": float(joblib.load(path).predict(row)[0]),
            "units": spec["units"],
            "readout": spec["readout"],
            "tier": spec["tier"],
            "evidence": spec["evidence"],
            "caveat": spec.get("caveat", ""),
            "training_sequences": spec["training_sequences"],
        }

    return {"available": bool(predictions), "reason": "", "predictions": predictions}
