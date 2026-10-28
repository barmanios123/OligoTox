"""Target-site accessibility.

For an RNase H gapmer the oligo itself is largely unstructured. What
plausibly gates cleavage is whether the *target* site is open or base-paired,
which is secondary structure of the transcript rather than any property of
the oligo. That is the proposal's "accessible binding regions", and unlike
3D folding it is computable today and falsifiable against knockdown data.

Accessibility is the probability that the binding window is unpaired, from
ViennaRNA's local partition function (`pfl_fold_up`, the RNAplfold
algorithm). Local folding is used rather than global because a transcript of
several kilobases does not fold as one domain, and a global MFE over it
would be a number without a meaning.

Returned per site:

    accessibility       P(the whole binding window is unpaired)
    log_accessibility   log10, because the raw value spans orders of
                        magnitude and models should see that linearly
    flank_accessibility the same for a window just upstream, which
                        distinguishes a locally open region from a lucky site
    relative_position   where in the transcript the site falls

A site that cannot be located in the transcript returns `site_found = 0`.
That is absence of evidence and the caller must not read the zeros as
"inaccessible".
"""
from __future__ import annotations

import math
from dataclasses import asdict, dataclass
from functools import lru_cache

# RNAplfold parameters. Window 80 with a 40-nt max span is the usual choice
# for mRNA accessibility and keeps the computation tractable on kilobase
# transcripts.
WINDOW = 80
MAX_BP_SPAN = 40
MAX_UNPAIRED = 25
FLANK_OFFSET = 30


@dataclass
class SiteAccessibility:
    site_found: bool
    start: int
    length: int
    accessibility: float
    log_accessibility: float
    flank_accessibility: float
    relative_position: float
    transcript_length: int

    def as_features(self) -> dict[str, float]:
        return {
            "site_found": float(self.site_found),
            "site_accessibility": self.accessibility,
            "site_log_accessibility": self.log_accessibility,
            "site_flank_accessibility": self.flank_accessibility,
            "site_relative_position": self.relative_position,
        }

    def to_dict(self) -> dict:
        return asdict(self)


def _rc(sequence: str) -> str:
    return sequence.upper().replace("U", "T").translate(
        str.maketrans("ACGT", "TGCA")
    )[::-1]


def find_site(oligo: str, transcript: str) -> int:
    """Where the oligo binds, or -1.

    An ASO binds the reverse complement of its own sequence. Tries the full
    length, then progressively trims from the 5' end, because a corpus
    transcript may be a different isoform from the one screened.
    """
    target = transcript.upper().replace("U", "T")
    probe = _rc(oligo)
    for trim in range(0, max(len(probe) - 12, 0) + 1):
        candidate = probe[trim:]
        position = target.find(candidate)
        if position != -1:
            return position
    return -1


@lru_cache(maxsize=64)
def _unpaired_profile(transcript: str) -> tuple:
    """Local unpaired probabilities for a transcript, cached per sequence."""
    import RNA

    rna = transcript.upper().replace("T", "U")
    profile = RNA.pfl_fold_up(rna, MAX_UNPAIRED, WINDOW, MAX_BP_SPAN)
    return tuple(tuple(row) for row in profile)


def _window_probability(profile: tuple, end_position: int, length: int) -> float:
    """P(unpaired) for a window of `length` ending at `end_position`, 1-based."""
    if length > MAX_UNPAIRED:
        length = MAX_UNPAIRED
    if end_position < 1 or end_position >= len(profile):
        return 0.0
    row = profile[end_position]
    if length >= len(row):
        return 0.0
    value = row[length]
    # ViennaRNA returns None where a window does not fit.
    return float(value) if value is not None and value == value else 0.0


def site_accessibility(oligo: str, transcript: str) -> SiteAccessibility:
    """Accessibility of the oligo's binding site in one transcript."""
    oligo = oligo.strip().upper().replace("U", "T")
    transcript = transcript.strip().upper().replace("U", "T")
    start = find_site(oligo, transcript)
    if start < 0 or not transcript:
        return SiteAccessibility(False, -1, len(oligo), 0.0, 0.0, 0.0, 0.0, len(transcript))

    length = min(len(oligo), MAX_UNPAIRED)
    profile = _unpaired_profile(transcript)

    end = start + length  # 1-based end index into the profile
    probability = _window_probability(profile, end, length)
    flank = _window_probability(profile, max(end - FLANK_OFFSET, length), length)

    return SiteAccessibility(
        site_found=True, start=start, length=length,
        accessibility=round(probability, 6),
        # Floor keeps the log finite for sites the model says are fully paired.
        log_accessibility=round(math.log10(max(probability, 1e-9)), 4),
        flank_accessibility=round(flank, 6),
        relative_position=round(start / max(len(transcript), 1), 4),
        transcript_length=len(transcript),
    )


def accessibility_features(oligo: str, transcript: str | None) -> dict[str, float]:
    """Features for a scoring pass, with zeros when no transcript is known."""
    if not transcript:
        return SiteAccessibility(False, -1, len(oligo), 0.0, 0.0, 0.0, 0.0, 0).as_features()
    return site_accessibility(oligo, transcript).as_features()


__all__ = [
    "SiteAccessibility", "accessibility_features", "find_site", "site_accessibility",
]
