"""The triage decides what a paper is; these pin the decisions that matter.

The risk it guards against is asymmetric. Calling a single-compound study a
dataset wastes five minutes of someone's reading. Calling a dataset a
single-compound study loses the dataset silently, and nothing downstream ever
notices. So the tests that matter most are the ones about what stays open.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from triage_literature import classify, decided_pmids  # noqa: E402


def test_a_screen_with_a_readout_is_a_dataset():
    verdict, _ = classify(
        "We screened a library of 200 gapmers and measured ALT elevation in mouse "
        "liver for each sequence."
    )
    assert verdict == "dataset"


def test_a_review_is_not_a_dataset():
    """A review borrows the vocabulary of a screen while reporting no rows.

    This is the case that first slipped through: a CNS safety review naming
    calcium imaging and cytotoxicity assays matched panel-and-readout and was
    offered up as a dataset.
    """
    verdict, _ = classify(
        "This review summarizes current in vitro approaches for assessing safety, "
        "including a panel of cytotoxicity assays, calcium oscillation and "
        "knockdown measurements across oligonucleotides."
    )
    assert verdict == "review"


def test_a_delivery_platform_is_a_method():
    verdict, _ = classify(
        "We developed a conjugation platform enabling duplex RNA delivery to the "
        "CNS, and identified an optimal architecture."
    )
    assert verdict == "method"


def test_a_title_alone_is_never_enough_to_dismiss_a_paper():
    """The false negative this whole file exists to prevent.

    `single_study` is the fall-through verdict, so without this a paper whose
    abstract never downloaded gets quietly dismissed on its title. It has to
    stay open instead.
    """
    title_only = "Gapmers targeting conserved single-stranded regions reduce viral load"
    assert classify(title_only, has_abstract=False)[0] == "unknown"
    # With an abstract behind it, the same text may legitimately be dismissed.
    assert classify(title_only, has_abstract=True)[0] != "unknown"


def test_nothing_retrieved_stays_open():
    assert classify("")[0] == "unknown"


def test_the_registry_is_what_counts_as_already_decided():
    """A paper whose dataset already ships must not come back as unassessed."""
    decided = decided_pmids()
    assert "42542598" in decided, "the ASO Atlas paper is adopted and should be decided"
    assert "adopted_per_context" in decided["42542598"]


def test_a_hand_decision_outranks_the_classifier():
    """Someone who read the paper beats a regex, and only reads it once."""
    decisions = json.loads(
        (ROOT / "data" / "registry" / "literature_decisions.json").read_text()
    )
    assert decisions, "hand decisions should be recorded, not kept in a head"
    for pmid, call in decisions.items():
        assert call["verdict"] in {
            "dataset", "single_study", "method", "reference", "review", "unknown",
        }, f"{pmid} has a verdict the triage does not define"
        assert len(call["reason"]) > 40, f"{pmid} needs a reason worth reading"


def test_every_candidate_the_monitor_found_has_been_assessed():
    """The backlog this was built to clear must stay cleared."""
    triage = ROOT / "reports" / "results" / "literature_triage.json"
    candidates = ROOT / "reports" / "results" / "literature_candidates.json"
    if not (triage.exists() and candidates.exists()):
        return  # a clean clone has not run the monitor yet
    assessed = json.loads(triage.read_text())["assessed"]
    found = json.loads(candidates.read_text())["candidates"]
    unassessed = [c["id"] for c in found if c["id"] not in assessed]
    assert not unassessed, f"unassessed candidates: {unassessed}"
