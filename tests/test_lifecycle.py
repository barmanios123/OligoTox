"""An autonomous loop is only safe if promotion and demotion are both real.

The ratchet: evidence can move an endpoint either way, but the gates never
loosen and nothing can pin a tier.
"""
import json

import pytest

from oligotox.lifecycle import (
    EndpointState,
    apply_transitions,
    describe_change,
    load_history,
    record_history,
)
from oligotox.lifecycle.transitions import previous_from_history


def _state(name, tier, **kw):
    return EndpointState(endpoint=name, tier=tier, **kw)


def test_more_evidence_promotes():
    before = {"ep": _state("ep", "insufficient_data", positive_sequences=4)}
    after = {"ep": _state("ep", "validated", positive_sequences=40)}
    (t,) = apply_transitions(before, after)
    assert t.direction == "promoted"
    assert t.now_presented and not t.was_presented
    assert "PROMOTED" in describe_change(t)


def test_losing_evidence_demotes_without_a_grace_period():
    """An endpoint that stops clearing its gates loses its tier immediately."""
    before = {"ep": _state("ep", "validated", positive_sequences=40)}
    after = {"ep": _state("ep", "source_confounded", positive_sequences=40)}
    (t,) = apply_transitions(before, after)
    assert t.direction == "demoted"
    assert t.was_presented and not t.now_presented
    assert "no longer presented" in describe_change(t)


def test_an_endpoint_that_vanishes_is_a_demotion():
    before = {"ep": _state("ep", "validated")}
    (t,) = apply_transitions(before, {})
    assert t.direction == "demoted" and t.after == "not_evaluated"
    assert t.was_presented and not t.now_presented


def test_a_new_endpoint_is_new_not_promoted():
    (t,) = apply_transitions({}, {"ep": _state("ep", "weak")})
    assert t.direction == "new"


def test_unchanged_endpoints_are_reported_but_not_logged(tmp_path):
    states = {"a": _state("a", "validated"), "b": _state("b", "weak")}
    transitions = apply_transitions(states, states)
    assert all(t.direction == "unchanged" for t in transitions)

    log = tmp_path / "history.jsonl"
    entry = record_history(log, transitions, trigger="test")
    assert entry["transitions"] == [], "only changes belong in the log"
    assert entry["presented"] == ["a"]


def test_history_is_append_only_and_replayable(tmp_path):
    log = tmp_path / "history.jsonl"
    record_history(log, apply_transitions({}, {"ep": _state("ep", "insufficient_data")}), "run 1")
    record_history(
        log,
        apply_transitions({"ep": _state("ep", "insufficient_data")},
                          {"ep": _state("ep", "validated", positive_sequences=40)}),
        "run 2",
    )
    entries = load_history(log)
    assert len(entries) == 2, "each run appends; nothing is rewritten"
    assert entries[1]["promoted"] == ["ep"]

    # The latest tier per endpoint is recoverable from the log alone.
    assert previous_from_history(log)["ep"].tier == "validated"


def test_the_shipped_history_matches_the_shipped_config(root):
    """The log must describe the config the app actually reads."""
    log = root / "reports/results/endpoint_history.jsonl"
    if not log.exists():
        pytest.skip("pipeline has not been run in this checkout")
    entries = load_history(log)
    assert entries, "history exists but is empty"

    evidence = json.loads(
        (root / "artifacts/config/endpoint_evidence.json").read_text()
    )["endpoints"]
    presented = sorted(e for e, v in evidence.items() if v["tier"] == "validated")
    assert entries[-1]["presented"] == presented


def test_there_is_no_way_to_force_a_tier():
    """Guard against someone adding an override later."""
    from oligotox.lifecycle import transitions as module

    source = module.__file__
    text = open(source).read().lower()
    for forbidden in ("def force_", "def pin_", "def override_", "manual_tier"):
        assert forbidden not in text, f"{forbidden} would break the ratchet"
