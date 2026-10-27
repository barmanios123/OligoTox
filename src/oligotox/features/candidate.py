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


def sequence_features(sequence: str, ps_fraction: float = 0.0, lna_fraction: float = 0.0) -> Dict[str, float]:
    """Features derived from sequence and backbone fractions alone.

    Training recomputes these for every row with this function, because the
    upstream sources filled them in inconsistently.
    """
    seq = sequence.strip().upper().replace("U", "T")
    length = len(seq)
    gc = _gc(seq)
    return {
        "seq_length": float(length),
        "gc_content": gc,
        # seed region = positions 2-8
        "seed_gc": _gc(seq[1:8]),
        "cpg_density": sum(seq[i:i + 2] == "CG" for i in range(length - 1)) / max(length - 1, 1),
        "tlr7_score": score_tlr_innate(seq, ps_fraction=ps_fraction).tlr7_score,
        "ps_x_len": ps_fraction * length,
        "lna_x_gc": lna_fraction * gc,
    }


def recompute_sequence_features(df: pd.DataFrame) -> pd.DataFrame:
    """Overwrite sequence-derived columns in a training table, row by row."""
    df = df.copy()
    has_seq = df["sequence"].fillna("").astype(str).str.strip().str.len() > 0
    ps = pd.to_numeric(df.get("ps_fraction", 0.0), errors="coerce").fillna(0.0)
    lna = pd.to_numeric(df.get("lna_fraction", 0.0), errors="coerce").fillna(0.0)
    rows = [
        sequence_features(s, p, l)
        for s, p, l in zip(df.loc[has_seq, "sequence"], ps[has_seq], lna[has_seq])
    ]
    if rows:
        recomputed = pd.DataFrame(rows, index=df.index[has_seq])
        for col in recomputed.columns:
            df.loc[has_seq, col] = recomputed[col]
    return df


def candidate_features(
    sequence: str,
    modality: str = "siRNA",
    chemistry: str = "PS",
    ps_fraction: float = 0.0,
    lna_fraction: float = 0.0,
    dose_mg_kg: float = 0.0,
) -> Dict[str, Any]:
    seq = sequence.strip().upper().replace("U", "T")
    row = pd.DataFrame([{
        "sequence": seq,
        "modality": modality,
        "chemistry": chemistry,
        "length_nt": len(seq),
        "gc_content": _gc(seq),
        "ps_fraction": ps_fraction,
        "lna_fraction": lna_fraction,
    }])

    feats: Dict[str, Any] = {
        "ps_fraction": ps_fraction,
        "lna_fraction": lna_fraction,
        "dose_mg_kg": dose_mg_kg,
        **sequence_features(seq, ps_fraction, lna_fraction),
        # Used by cns_toxicity, which has its own sequence-only feature list.
        **g_motif_features(seq),
    }
    for block in (
        build_aso_rich_features(row),
        build_structure_features(row, cache_path=None),
        build_embedding_features(row, cache_path=None),
    ):
        feats.update(block.iloc[0].to_dict())
    return feats

def g_motif_features(sequence: str) -> Dict[str, float]:
    """G-content features for CNS neurotoxicity.

    Hagedorn et al. (Nucleic Acid Ther. 2022) report that the number of G
    nucleotides and the distance from the 3' end to the nearest G both
    influence acute neurotoxicity. Training and scoring both call this, so the
    two stay in step.
    """
    seq = sequence.strip().upper().replace("U", "T")
    n = len(seq)
    g_count = seq.count("G")
    nearest = next((i for i, ch in enumerate(reversed(seq)) if ch == "G"), n)
    runs = "".join(c if c == "G" else " " for c in seq).split()
    return {
        "g_count": float(g_count),
        "g_fraction": g_count / n if n else 0.0,
        "dist_3p_to_g": float(nearest),
        "longest_g_run": float(max((len(r) for r in runs), default=0)),
    }


def cns_sequence_features(sequence: str) -> Dict[str, float]:
    """Every feature the externally trained cns_toxicity model uses.

    Sequence only: the external CNS assays carry no chemistry annotation, so
    this model never sees a chemistry or modality feature.
    """
    feats = sequence_features(sequence)
    feats.pop("ps_x_len", None)
    feats.pop("lna_x_gc", None)
    feats.update(g_motif_features(sequence))

    seq = sequence.strip().upper().replace("U", "T")
    struct = build_structure_features(pd.DataFrame({"sequence": [seq]}))
    for col in struct.columns:
        feats[col] = float(struct.iloc[0][col])
    return feats
