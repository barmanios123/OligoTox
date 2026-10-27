"""Design-time features for a single candidate, matching the trained columns.

Reuses the same builders as the training pipeline so the app and the models
see features computed the same way.
"""
from __future__ import annotations

from typing import Any, Dict

import pandas as pd

from oligotox.features.aso_rich_features import build_aso_rich_features
from oligotox.features.embedding_features import build_embedding_features
from oligotox.features.structure_features import build_structure_features
from oligotox.integrations.tlr_immune import score_tlr_innate


def _gc(seq: str) -> float:
    return (seq.count("G") + seq.count("C")) / len(seq) if seq else 0.0


def candidate_features(
    sequence: str,
    modality: str = "siRNA",
    chemistry: str = "PS",
    ps_fraction: float = 0.0,
    lna_fraction: float = 0.0,
    dose_mg_kg: float = 0.0,
) -> Dict[str, Any]:
    seq = sequence.strip().upper().replace("U", "T")
    length = len(seq)
    gc = _gc(seq)
    row = pd.DataFrame([{
        "sequence": seq,
        "modality": modality,
        "chemistry": chemistry,
        "length_nt": length,
        "gc_content": gc,
        "ps_fraction": ps_fraction,
        "lna_fraction": lna_fraction,
    }])

    feats: Dict[str, Any] = {
        "seq_length": float(length),
        "gc_content": gc,
        "ps_fraction": ps_fraction,
        "lna_fraction": lna_fraction,
        # seed region = positions 2-8
        "seed_gc": _gc(seq[1:8]),
        "cpg_density": sum(seq[i:i + 2] == "CG" for i in range(length - 1)) / max(length - 1, 1),
        "tlr7_score": score_tlr_innate(seq, ps_fraction=ps_fraction).tlr7_score,
        "ps_x_len": ps_fraction * length,
        "lna_x_gc": lna_fraction * gc,
        "dose_mg_kg": dose_mg_kg,
    }
    for block in (
        build_aso_rich_features(row),
        build_structure_features(row, cache_path=None),
        build_embedding_features(row, cache_path=None),
    ):
        feats.update(block.iloc[0].to_dict())
    return feats
