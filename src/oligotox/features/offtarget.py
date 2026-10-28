"""Hybridisation-dependent off-target burden.

The published mechanism for gapmer hepatotoxicity is not chemistry. Burel et
al. and Kasuya et al. showed it is RNase H1 cleaving unintended transcripts
the oligo is complementary to, particularly long pre-mRNAs. That makes
complementarity burden a design-time feature, computable from sequence
alone, and it is the one mechanistic feature this project has been missing.

What was here before did not compute it. `integrations/seed_match.py` takes
an optional transcriptome, nothing ever supplied one, and it fell through to
`est = 0.01 + gc * 0.02` — GC content rescaled and returned as a hit count.
Nothing trained on it, so nothing was wrong, but the capability was absent.

This counts real matches against a real background set:

    full_hits      exact matches of the full reverse complement
    gapmer_hits    matches of the central DNA gap region, which is what
                   RNase H1 actually requires to cleave
    seed_hits      matches of positions 2-8, the miRNA-like seed

Each is reported as a raw count and as a rate per megabase, so a background
of a different size stays comparable.

**This is a relative measure, not a genome-wide one.** The background is
whatever transcripts are supplied. Counts are comparable between oligos
scored against the same background and are not absolute off-target risk.
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from functools import lru_cache
from pathlib import Path
from typing import Iterable, Sequence

COMPLEMENT = str.maketrans("ACGTU", "TGCAA")
# The central DNA gap of a typical 5-10-5 gapmer, which RNase H1 needs.
DEFAULT_GAP = (5, 15)
SEED = (1, 8)


def reverse_complement(sequence: str) -> str:
    return sequence.upper().replace("U", "T").translate(COMPLEMENT)[::-1]


@dataclass
class OffTargetBurden:
    full_hits: int
    gapmer_hits: int
    seed_hits: int
    background_megabases: float
    full_per_mb: float
    gapmer_per_mb: float
    seed_per_mb: float
    background_transcripts: int

    def as_features(self) -> dict[str, float]:
        return {
            "offtarget_full_hits": float(self.full_hits),
            "offtarget_gapmer_hits": float(self.gapmer_hits),
            "offtarget_seed_hits": float(self.seed_hits),
            "offtarget_full_per_mb": self.full_per_mb,
            "offtarget_gapmer_per_mb": self.gapmer_per_mb,
            "offtarget_seed_per_mb": self.seed_per_mb,
        }

    def to_dict(self) -> dict:
        return asdict(self)


@lru_cache(maxsize=4)
def load_background(directory: str | Path) -> tuple[str, int, float]:
    """Concatenate cached transcripts into one searchable background.

    Joined with a separator so a match cannot straddle two transcripts.
    """
    directory = Path(directory)
    if not directory.exists():
        raise FileNotFoundError(
            f"no transcript cache at {directory}. Run scripts/fetch_transcripts.py first."
        )
    chunks, count = [], 0
    for path in sorted(directory.glob("*.json")):
        if path.name.startswith("_"):
            continue
        try:
            record = json.loads(path.read_text())
        except json.JSONDecodeError:
            continue
        sequence = (record.get("sequence") or "").upper().replace("U", "T")
        if sequence:
            chunks.append(sequence)
            count += 1
    if not chunks:
        raise ValueError(f"transcript cache at {directory} is empty")
    joined = "N".join(chunks)
    return joined, count, len(joined) / 1_000_000


def _count(haystack: str, needle: str) -> int:
    if not needle or len(needle) > len(haystack):
        return 0
    total = start = 0
    while True:
        found = haystack.find(needle, start)
        if found == -1:
            return total
        total += 1
        start = found + 1  # overlapping matches count


def burden(
    sequence: str,
    background: str,
    n_transcripts: int,
    megabases: float,
    gap: tuple[int, int] = DEFAULT_GAP,
) -> OffTargetBurden:
    """Count complementarity hits for one oligo against one background."""
    oligo = sequence.strip().upper().replace("U", "T")
    full_rc = reverse_complement(oligo)
    gap_start, gap_end = gap
    gap_rc = reverse_complement(oligo[gap_start:gap_end]) if len(oligo) > gap_end else ""
    seed_rc = reverse_complement(oligo[SEED[0]:SEED[1]])

    full = _count(background, full_rc)
    gapmer = _count(background, gap_rc) if gap_rc else 0
    seed = _count(background, seed_rc)
    scale = megabases if megabases > 0 else 1.0

    return OffTargetBurden(
        full_hits=full, gapmer_hits=gapmer, seed_hits=seed,
        background_megabases=round(scale, 3),
        full_per_mb=round(full / scale, 4),
        gapmer_per_mb=round(gapmer / scale, 4),
        seed_per_mb=round(seed / scale, 4),
        background_transcripts=n_transcripts,
    )


def offtarget_features(
    sequence: str,
    transcript_dir: str | Path = "data/external/transcripts",
    gap: tuple[int, int] = DEFAULT_GAP,
) -> dict[str, float]:
    """Design-time off-target features, or zeros when no background is cached.

    Returning zeros rather than raising keeps scoring working without the
    cache, but the caller should check `offtarget_available` before using
    them: all-zero features are absence of evidence, not absence of risk.
    """
    try:
        background, count, megabases = load_background(transcript_dir)
    except (FileNotFoundError, ValueError):
        return {
            **OffTargetBurden(0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0).as_features(),
            "offtarget_available": 0.0,
        }
    result = burden(sequence, background, count, megabases, gap)
    return {**result.as_features(), "offtarget_available": 1.0}


def burden_for_many(
    sequences: Iterable[str],
    transcript_dir: str | Path = "data/external/transcripts",
    gap: tuple[int, int] = DEFAULT_GAP,
) -> list[OffTargetBurden]:
    background, count, megabases = load_background(transcript_dir)
    return [burden(s, background, count, megabases, gap) for s in sequences]


__all__ = [
    "OffTargetBurden", "burden", "burden_for_many", "load_background",
    "offtarget_features", "reverse_complement",
]
