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


def load_context_models(config_dir: str | Path) -> dict:
    """Per-assay potency models, from scripts/train_context_potency.py.

    ASO Atlas pools percent inhibition across 81 cell lines, and pooled it
    predicts worse than the far smaller single-study model. Each assay that
    clears its own baselines gets a model instead, so a prediction always
    names the assay it belongs to.
    """
    path = Path(config_dir) / "potency_contexts.json"
    return json.loads(path.read_text()) if path.exists() else {}


def _accessibility_features(sequence: str, transcript: str, cols: list[str]):
    """Accessibility for a located site, or None when it does not locate.

    Returning None rather than zeros is the whole point. A zero here would
    mean "not computable" while the model reads it as "fully paired", which
    is the source-indicator failure in NEGATIVE_RESULTS entry 13.
    """
    from oligotox.features.accessibility import site_accessibility

    site = site_accessibility(sequence, transcript)
    if not site.site_found:
        return None
    feats = site.as_features()
    return {c: feats[c] for c in cols if c in feats}


def score_efficacy(
    sequence: str,
    modality: str,
    model_dir: str | Path,
    config_dir: str | Path,
    target_transcript: str | None = None,
) -> Dict[str, Any]:
    """Potency for the matching modality, or an explanation of why not.

    The returned score is the model's raw readout, in the units the training
    data used, and those units are returned alongside it. They are not a
    probability and not comparable across modalities.

    `target_transcript` is optional. When a context ships an accessibility
    model and the oligo's binding site locates in the transcript, that model
    is used instead and the prediction says so. Without a transcript, or when
    the site does not locate, the design-only model is used unchanged. The
    accessibility columns are never filled in with zeros to make up the
    difference.
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

    # Per-assay models, each predicting knockdown in its own context only.
    contexts = load_context_models(config_dir).get("contexts", {})
    for name, spec in contexts.items():
        if spec.get("modality", "ASO").lower() != str(modality).lower():
            continue
        path = model_dir / spec["model_file"]
        if not path.exists():
            continue
        cols = spec["feature_columns"]
        missing = [c for c in cols if c not in feats]
        if missing:
            continue
        # Prefer the accessibility model when a transcript was supplied and
        # the site locates. Everything needed is resolved first, so the
        # fallback path stays exactly what it was before this existed.
        use_path, use_cols, use_feats = path, cols, feats
        used_accessibility = False
        evidence = spec["evidence"]

        access_spec = spec.get("accessibility_model")
        if target_transcript and access_spec:
            acc_path = model_dir / access_spec["model_file"]
            acc_cols = access_spec["feature_columns"]
            site_feats = (
                _accessibility_features(sequence, target_transcript, acc_cols)
                if acc_path.exists() else None
            )
            if site_feats is not None:
                merged = {**feats, **site_feats}
                if all(c in merged for c in acc_cols):
                    use_path, use_cols, use_feats = acc_path, acc_cols, merged
                    used_accessibility = True
                    evidence = access_spec["evidence"]

        row = pd.DataFrame([[float(use_feats[c]) for c in use_cols]], columns=use_cols)
        predictions[f"potency_{name}"] = {
            "score": float(joblib.load(use_path).predict(row)[0]),
            "used_accessibility": used_accessibility,
            "units": spec["units"],
            "readout": spec["readout"],
            "tier": "context",
            "context": spec["context"],
            "evidence": evidence,
            "caveat": spec.get("caveat", ""),
            "training_sequences": spec["sequences"],
        }

    return {"available": bool(predictions), "reason": "", "predictions": predictions}
