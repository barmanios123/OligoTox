from __future__ import annotations

from pathlib import Path
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


def _gc_fraction(seq: str) -> float:
    if not seq:
        return 0.0
    n = len(seq)
    gc = sum(1 for ch in seq if ch in {"G", "C"})
    return gc / n if n else 0.0


def _longest_homopolymer(seq: str) -> int:
    if not seq:
        return 0
    best = cur = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            cur += 1
            best = max(best, cur)
        else:
            cur = 1
    return best


def _reverse_complement_rna(seq: str) -> str:
    table = str.maketrans({"A": "U", "U": "A", "T": "A", "C": "G", "G": "C", "N": "N"})
    return seq.translate(table)[::-1]


def _palindrome_fraction(seq: str, k: int = 4) -> float:
    if len(seq) < k:
        return 0.0
    rc = _reverse_complement_rna(seq)
    hits = 0
    total = max(1, len(seq) - k + 1)
    windows = {seq[i:i+k] for i in range(len(seq) - k + 1)}
    for i in range(len(rc) - k + 1):
        if rc[i:i+k] in windows:
            hits += 1
    return hits / total


def _dinucleotide_repeats(seq: str) -> int:
    if len(seq) < 4:
        return 0
    count = 0
    for i in range(len(seq) - 3):
        if seq[i:i+2] == seq[i+2:i+4]:
            count += 1
    return count


def _simple_hairpin_score(seq: str) -> float:
    if len(seq) < 8:
        return 0.0
    left = seq[: len(seq) // 2]
    right = seq[len(seq) // 2 :]
    rc_right = _reverse_complement_rna(right)
    overlap = sum(1 for a, b in zip(left, rc_right) if a == b)
    return overlap / max(1, min(len(left), len(rc_right)))


def _placeholder_mfe(seq: str) -> float:
    if not seq:
        return 0.0
    gc = _gc_fraction(seq)
    hp = _simple_hairpin_score(seq)
    return -1.0 * len(seq) * (0.15 + 0.45 * gc + 0.40 * hp)


def build_structure_features(
    df: pd.DataFrame,
    cache_path: str | Path | None = "data/processed/structure_features_cache.csv",
    use_viennarna: bool = True,
) -> pd.DataFrame:
    seq_col = _find_sequence_column(df)
    out = pd.DataFrame(index=df.index)

    if seq_col is None:
        cols = [
            "struct_seq_len", "struct_gc_fraction", "struct_uracil_fraction",
            "struct_longest_homopolymer", "struct_palindrome_fraction",
            "struct_dinucleotide_repeats", "struct_hairpin_score",
            "struct_mfe_proxy", "struct_viennarna_used",
        ]
        for c in cols:
            out[c] = 0.0
        return out

    seqs = df[seq_col].map(_sanitize_seq)

    out["struct_seq_len"] = seqs.map(len).astype(float)
    out["struct_gc_fraction"] = seqs.map(_gc_fraction).astype(float)
    out["struct_uracil_fraction"] = seqs.map(lambda s: (s.count("U") + s.count("T")) / len(s) if len(s) else 0.0).astype(float)
    out["struct_longest_homopolymer"] = seqs.map(_longest_homopolymer).astype(float)
    out["struct_palindrome_fraction"] = seqs.map(_palindrome_fraction).astype(float)
    out["struct_dinucleotide_repeats"] = seqs.map(_dinucleotide_repeats).astype(float)
    out["struct_hairpin_score"] = seqs.map(_simple_hairpin_score).astype(float)
    out["struct_mfe_proxy"] = seqs.map(_placeholder_mfe).astype(float)
    out["struct_viennarna_used"] = 0.0

    if use_viennarna:
        # No silent fallback. A placeholder MFE where real thermodynamics were
        # expected is a feature drift nothing would catch: training and scoring
        # would disagree while both looked fine.
        try:
            import RNA  # type: ignore
        except ImportError as exc:  # pragma: no cover - requirements pin it
            raise ImportError(
                "ViennaRNA is required for structure features and is pinned in "
                "requirements.txt. Without it struct_mfe_proxy silently becomes a "
                "GC-derived placeholder, so training and scoring would drift apart. "
                "Install it, or pass use_viennarna=False to accept the placeholder "
                "deliberately."
            ) from exc
        try:
            out["struct_mfe_proxy"] = seqs.map(lambda s: float(RNA.fold(s.replace("T", "U"))[1]) if s else 0.0)
            out["struct_viennarna_used"] = 1.0
        except Exception:
            pass

    if cache_path is not None:
        cache_path = Path(cache_path)
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(cache_path, index=False)

    return out
