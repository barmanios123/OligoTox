from __future__ import annotations

from pathlib import Path
import hashlib
import numpy as np
import pandas as pd


def _find_sequence_column(df: pd.DataFrame) -> str | None:
    candidates = [
        "sequence", "oligo_sequence", "sense_sequence", "antisense_sequence",
        "guide_sequence", "seq", "raw_sequence"
    ]
    for c in candidates:
        if c in df.columns:
            return c
    return None


def _sanitize_seq(x: object) -> str:
    s = "" if x is None else str(x).upper().strip()
    allowed = set("ACGUTN")
    return "".join(ch for ch in s if ch in allowed)


def _hash_to_unit_vector(seq: str, dim: int = 16) -> np.ndarray:
    if not seq:
        return np.zeros(dim, dtype=float)
    vals = []
    for i in range(dim):
        h = hashlib.sha256(f"{seq}|{i}".encode()).hexdigest()
        vals.append(int(h[:8], 16) / 0xFFFFFFFF)
    arr = np.array(vals, dtype=float)
    norm = np.linalg.norm(arr)
    if norm == 0:
        return arr
    return arr / norm


def build_embedding_features(
    df: pd.DataFrame,
    cache_path: str | Path | None = "data/processed/embedding_features_cache.csv",
    model_name: str = "stub_hash_rna",
    dim: int = 16,
) -> pd.DataFrame:
    seq_col = _find_sequence_column(df)
    out = pd.DataFrame(index=df.index)

    if seq_col is None:
        for i in range(dim):
            out[f"embed_{model_name}_{i:02d}"] = 0.0
        out["embed_model_available"] = 0.0
        out["embed_seq_present"] = 0.0
        return out

    seqs = df[seq_col].map(_sanitize_seq)
    mat = np.vstack([_hash_to_unit_vector(seq, dim=dim) for seq in seqs])

    for i in range(dim):
        out[f"embed_{model_name}_{i:02d}"] = mat[:, i]

    out["embed_model_available"] = 0.0
    out["embed_seq_present"] = seqs.map(lambda s: 1.0 if len(s) > 0 else 0.0).astype(float)

    if cache_path is not None:
        cache_path = Path(cache_path)
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(cache_path, index=False)

    return out
