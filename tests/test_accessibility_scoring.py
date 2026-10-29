"""Accessibility may change a prediction only when it was really computed.

The rule these enforce is the one from NEGATIVE_RESULTS entry 13. A sequence
whose binding site does not locate has *unknown* accessibility, not zero
accessibility, and a zero fed to the model would read as "fully paired". So a
missing or unlocatable target must fall back to the design-only model rather
than filling the columns in.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from oligotox.efficacy.potency import load_context_models, score_efficacy

ROOT = Path(__file__).resolve().parent.parent
import sys  # noqa: E402
sys.path.insert(0, str(ROOT / "scripts"))
MODELS = ROOT / "artifacts" / "models"
CONFIG = ROOT / "artifacts" / "config"
GENOMIC = ROOT / "data" / "external" / "genomic"


def _context_with_accessibility():
    contexts = load_context_models(CONFIG).get("contexts", {})
    for name, spec in contexts.items():
        if spec.get("accessibility_model"):
            return name, spec
    return None, None


def _a_located_pair():
    """A (sequence, transcript) whose site actually locates, from the corpus."""
    from oligotox.features.accessibility import find_site

    for path in sorted(GENOMIC.glob("*.json"))[:40]:
        if path.name.startswith("_"):
            continue
        record = json.loads(path.read_text())
        seq = (record.get("sequence") or "").upper().replace("U", "T")
        if len(seq) < 2000:
            continue
        # Take a window out of the transcript and reverse complement it, which
        # is by construction an oligo that binds there.
        window = seq[1000:1018]
        oligo = window.translate(str.maketrans("ACGT", "TGCA"))[::-1]
        if find_site(oligo, seq) >= 0:
            return oligo, seq
    return None, None


@pytest.fixture(scope="module")
def shipped():
    name, spec = _context_with_accessibility()
    if not name:
        pytest.skip("no context ships an accessibility model yet")
    return name, spec


def test_without_a_transcript_nothing_changes():
    """The default path must be exactly what it was before this feature."""
    out = score_efficacy("ACGTACGTACGTACGTACGT", "ASO", MODELS, CONFIG)
    for pred in out.get("predictions", {}).values():
        assert pred.get("used_accessibility") in (False, None)


def test_an_unlocatable_site_falls_back_rather_than_zero_filling(shipped):
    """The entry 13 rule, enforced at the scoring boundary."""
    # Poly-A against a transcript it cannot bind.
    oligo = "A" * 20
    transcript = "GC" * 5000
    out = score_efficacy(oligo, "ASO", MODELS, CONFIG, target_transcript=transcript)
    for pred in out.get("predictions", {}).values():
        assert not pred.get("used_accessibility"), (
            "a site that does not locate must not be scored as if it had"
        )


def test_a_located_site_uses_the_accessibility_model(shipped):
    oligo, transcript = _a_located_pair()
    if oligo is None:
        pytest.skip("no locatable oligo could be constructed from the cache")
    out = score_efficacy(oligo, "ASO", MODELS, CONFIG, target_transcript=transcript)
    used = [n for n, p in out.get("predictions", {}).items()
            if p.get("used_accessibility")]
    assert used, "a located site should reach at least one accessibility model"


def test_the_accessibility_model_beat_the_plain_one_on_the_same_rows(shipped):
    """It ships only if it won a like-for-like comparison at training time."""
    _, spec = shipped
    acc = spec["accessibility_model"]
    assert acc["gain"] > 0, acc
    assert acc["spearman"] > acc["design_only_on_same_rows"]
    assert acc["rows"] > 0 and acc["genes"] >= 5


def test_only_the_contexts_that_survived_a_reseed_have_one():
    """Four contexts gained once and lost it on a different sample seed."""
    from train_context_potency import ACCESS_CANDIDATES

    contexts = load_context_models(CONFIG).get("contexts", {})
    with_access = {n for n, s in contexts.items() if s.get("accessibility_model")}
    assert with_access <= ACCESS_CANDIDATES, (
        f"{with_access - ACCESS_CANDIDATES} ship accessibility without having "
        "held across seeds"
    )
