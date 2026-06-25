"""oligotox.ml_predict
~~~~~~~~~~~~~~~~~~~~
Loads v4 RF models from artifacts/models_v4/ and artifacts/models_neurotox/
and generates per-endpoint toxicity probability scores for a given nucleotide sequence.

Usage:
    from oligotox.ml_predict import predict_sequence
    results = predict_sequence("TGGCATGATGCCTCGTAGA")
    # -> [{'endpoint': 'hepatocellular_injury', 'prob': 0.34, ...}, ...]
"""
from __future__ import annotations
from itertools import product
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

try:
    import joblib as _joblib
    _HAS_JOBLIB = True
except ImportError:
    import pickle as _pickle  # type: ignore
    _HAS_JOBLIB = False

# Resolved relative to repo root
REPO_ROOT = Path(__file__).resolve().parents[2]
MODEL_DIRS = [
    REPO_ROOT / "artifacts" / "models_v4",
    REPO_ROOT / "artifacts" / "models_neurotox",
]

# ----------- feature engineering -----------
def _kmer(seq: str, k: int) -> dict[str, float]:
    nts = "ACGTU"
    kmers = ["".join(p) for p in product(nts, repeat=k)]
    n = len(seq)
    return {f"k{k}_{km}": seq.count(km) / max(n - k + 1, 1) for km in kmers}


def featurise(sequence: str) -> dict[str, float]:
    """Return a flat dict of k-mer + physicochemical features for one sequence."""
    seq = sequence.upper().replace(" ", "").replace("-", "")
    n = len(seq)
    gc = (seq.count("G") + seq.count("C")) / max(n, 1)
    t = max(n // 3, 1)
    feats: dict[str, float] = {
        "length_nt": float(n),
        "gc_content": gc,
        "gc_5p": (seq[:t].count("G") + seq[:t].count("C")) / max(t, 1),
        "gc_mid": (seq[t:2*t].count("G") + seq[t:2*t].count("C")) / max(t, 1),
        "gc_3p": (seq[2*t:].count("G") + seq[2*t:].count("C")) / max(len(seq[2*t:]), 1),
    }
    feats.update(_kmer(seq, 2))
    feats.update(_kmer(seq, 3))
    return feats

_featurise = featurise  # backwards-compat alias

# ------------ model loading (cached) ------------
_MODEL_CACHE: dict[str, Any] = {}


def _load_models() -> dict[str, Any]:
    """Return dict of {endpoint_label: model} lazily loaded from disk."""
    global _MODEL_CACHE
    if _MODEL_CACHE:
        return _MODEL_CACHE
    for model_dir in MODEL_DIRS:
        if not model_dir.exists():
            continue
        # Support both .joblib (fast) and .pkl (legacy)
        for ext in ("*.joblib", "*.pkl"):
            for model_file in sorted(model_dir.glob(ext)):
                label = model_file.stem
                if label in _MODEL_CACHE:
                    continue  # first found wins (joblib preferred)
                try:
                    if _HAS_JOBLIB:
                        model = _joblib.load(model_file)
                    else:
                        with open(model_file, "rb") as f:
                            model = _pickle.load(f)
                        # Unwrap dict-wrapped models (e.g. neurotox {model, feature_names, threshold})
                    if isinstance(model, dict) and "model" in model:
                        import numpy as _np
                        feat_arr = model.get("feature_names")
                        inner = model["model"]
                        if feat_arr is not None and not hasattr(inner, "feature_names_in_"):
                            inner.feature_names_in_ = _np.array(list(feat_arr))
                        model = inner
                    _MODEL_CACHE[label] = model
                except Exception as e:
                    print(f"  [warn] could not load {model_file.name}: {e}")
    return _MODEL_CACHE


# ------------ public API ------------
def predict_sequence(sequence: str) -> list[dict]:
    """
    Predict toxicity probabilities for a single nucleotide sequence.

    Returns a list of dicts, one per loaded model endpoint:
        [
          {'endpoint': 'hepatocellular_injury', 'prob': 0.42,
           'risk_level': 'Medium', 'model_dir': 'models_v4'},
          ...
        ]
    """
    feat = featurise(sequence)
    models = _load_models()
    if not models:
        return []

    results = []
    for label, model in models.items():
        try:
            feat_df = pd.DataFrame([feat])
            # Align to model's expected features
            if hasattr(model, "feature_names_in_"):
                cols = list(model.feature_names_in_)
                missing = {c: 0.0 for c in cols if c not in feat_df.columns}
                if missing:
                    feat_df = pd.concat([feat_df, pd.DataFrame([missing])], axis=1)
                feat_df = feat_df[cols]
            else:
                feat_df = feat_df.reindex(
                    sorted(feat_df.columns), axis=1, fill_value=0.0
                )

            if hasattr(model, "predict_proba"):
                prob = float(model.predict_proba(feat_df)[0, 1])
            else:
                prob = float(model.predict(feat_df)[0])

            if prob >= 0.70:
                risk = "High"
            elif prob >= 0.40:
                risk = "Medium"
            else:
                risk = "Low"

            results.append(
                {
                    "endpoint": label,
                    "prob": round(prob, 4),
                    "prob_%": round(prob * 100, 1),
                    "risk_level": risk,
                }
            )
        except Exception as e:
            results.append(
                {"endpoint": label, "prob": None, "prob_%": None,
                 "risk_level": "Error", "error": str(e)}
            )

    results.sort(key=lambda x: (x["prob"] or 0), reverse=True)
    return results


def predict_batch(sequences: list[str]) -> pd.DataFrame:
    """Run predict_sequence on a list; return tidy DataFrame."""
    rows = []
    for seq in sequences:
        for r in predict_sequence(seq):
            r["sequence"] = seq
            rows.append(r)
    return pd.DataFrame(rows)


if __name__ == "__main__":
    import sys
    seq = sys.argv[1] if len(sys.argv) > 1 else "GCCTGTTGTGAGCCTGATAGTGACAG"
    print(f"Sequence: {seq}")
    for r in predict_sequence(seq):
        prob_str = f"{r['prob']:.3f}" if r['prob'] is not None else 'N/A'; print(f"  {r['endpoint']:<35s}  prob={prob_str}  risk={r['risk_level']}")


def get_top_features(sequence: str, endpoint: str, n: int = 10) -> list[dict]:
    """Return the top-n most important features for a given endpoint model.

    Args:
        sequence: Nucleotide sequence string.
        endpoint: Model label (e.g. 'hepatocellular_injury').
        n: Number of top features to return.

    Returns:
        List of dicts: [{"feature": str, "importance": float, "value": float}]
    """
    import numpy as _np
    models = _load_models()
    model = models.get(endpoint)
    if model is None:
        return []
    try:
        feat = featurise(sequence)
        feat_df = pd.DataFrame([feat]).fillna(0.0)
        if hasattr(model, "feature_names_in_"):
            cols = list(model.feature_names_in_)
            missing = {c: 0.0 for c in cols if c not in feat_df.columns}
            if missing:
                feat_df = pd.concat([feat_df, pd.DataFrame([missing])], axis=1)
            feat_df = feat_df[cols]
        # Support both tree-based (feature_importances_) and linear (coef_) models
        importances = getattr(model, "feature_importances_", None)
        if importances is None and hasattr(model, "coef_"):
            importances = np.abs(model.coef_).flatten()
        if importances is None:
            return []
        feature_names = list(feat_df.columns)
        top_idx = np.argsort(importances)[::-1][:n]
        return [
            {
                "feature": feature_names[i],
                "importance": float(importances[i]),
                "value": float(feat_df.iloc[0, i]),
            }
            for i in top_idx
        ]
    except Exception:
        return []
