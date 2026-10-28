"""Every health check must fail when the thing it watches is broken.

A monitoring script is only worth its runtime if it goes red. The failure
these guard against is the one the script itself was written to catch: a
check that reads the wrong field, finds nothing, and reports success. The
first draft of `evidence_consistent` did exactly that — it looked for a
`presented` key that does not exist, found no endpoints, and said "all
validated: []".

So each check is driven here against a deliberately broken state, and the
assertion is that it returns False.
"""
from __future__ import annotations

import json
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

import loop_health  # noqa: E402


@pytest.fixture
def sandbox(tmp_path, monkeypatch):
    """Point the module at a throwaway tree so nothing real is touched."""
    monkeypatch.setattr(loop_health, "ROOT", tmp_path)
    monkeypatch.setattr(loop_health, "RESULTS", tmp_path / "reports" / "results")
    monkeypatch.setattr(loop_health, "HISTORY", tmp_path / "history.jsonl")
    monkeypatch.setattr(loop_health, "REGISTRY", tmp_path / "sources.json")
    monkeypatch.setattr(loop_health, "INVENTORY", tmp_path / "inventory.json")
    monkeypatch.setattr(loop_health, "TRIAGE", tmp_path / "triage.json")
    monkeypatch.setattr(loop_health, "CANDIDATES", tmp_path / "candidates.json")
    monkeypatch.setattr(loop_health, "EVIDENCE", tmp_path / "evidence.json")
    monkeypatch.setattr(loop_health, "HEALTH", tmp_path / "health.json")
    return tmp_path


def _history(path: Path, *, days_ago: float, presented: list[str]) -> None:
    stamp = datetime.now(timezone.utc) - timedelta(days=days_ago)
    path.write_text(json.dumps({
        "timestamp": stamp.isoformat(), "trigger": "test",
        "presented": presented, "promoted": [], "demoted": [], "transitions": [],
    }) + "\n")


def test_a_stalled_loop_is_caught(sandbox):
    """The whole point: a schedule that stopped running must go red."""
    _history(loop_health.HISTORY, days_ago=40, presented=["a"])
    ok, detail = loop_health.check_history_fresh(14)
    assert not ok and "40 days ago" in detail


def test_a_fresh_run_passes(sandbox):
    _history(loop_health.HISTORY, days_ago=2, presented=["a"])
    assert loop_health.check_history_fresh(14)[0]


def test_a_loop_that_never_ran_is_caught(sandbox):
    assert not loop_health.check_history_fresh(14)[0]


def test_lost_data_is_caught(sandbox):
    """A source declared present and missing is how disk loss shows up."""
    (sandbox / "sources.json").write_text(json.dumps({"sources": [
        {"name": "gone", "status": "present", "local_path": "data/gone.csv"},
    ]}))
    ok, detail = loop_health.check_sources_on_disk()
    assert not ok and "gone" in detail


def test_a_source_waiting_on_a_human_is_not_a_failure(sandbox):
    """`fetchable` means a script cannot get it unattended. That is known, not broken."""
    (sandbox / "sources.json").write_text(json.dumps({"sources": [
        {"name": "manual", "status": "fetchable", "local_path": "data/manual.csv"},
    ]}))
    ok, detail = loop_health.check_sources_on_disk()
    assert ok and "manual" in detail


def test_having_nothing_to_train_on_is_caught(sandbox):
    (sandbox / "inventory.json").write_text(
        json.dumps({"usable_sources": 0, "usable_sequences": 0}))
    assert not loop_health.check_training_ran()[0]


def test_an_unassessed_backlog_is_caught(sandbox):
    (sandbox / "candidates.json").write_text(json.dumps(
        {"candidates": [{"id": "1"}, {"id": "2"}]}))
    (sandbox / "triage.json").write_text(json.dumps({"assessed": {"1": {}}}))
    ok, detail = loop_health.check_triage_clear()
    assert not ok and "1 of 2" in detail


def test_a_schema_change_fails_rather_than_passing_vacuously(sandbox):
    """The bug this file exists for.

    If the config stops looking the way the check expects, the check must go
    red. Reading no endpoints and concluding the gate holds is the one outcome
    that must be impossible.
    """
    (sandbox / "evidence.json").write_text(json.dumps({"something_else": {}}))
    ok, detail = loop_health.check_evidence_consistent()
    assert not ok and "schema changed" in detail


def test_an_endpoint_with_no_tier_is_caught(sandbox):
    (sandbox / "evidence.json").write_text(json.dumps(
        {"endpoints": {"a": {"tier": "validated"}, "b": {}}}))
    ok, detail = loop_health.check_evidence_consistent()
    assert not ok and "b" in detail


def test_config_and_log_disagreeing_is_caught(sandbox):
    """Two artifacts written by one pipeline must not drift apart."""
    (sandbox / "evidence.json").write_text(json.dumps(
        {"endpoints": {"a": {"tier": "validated"}, "b": {"tier": "weak"}}}))
    _history(loop_health.HISTORY, days_ago=1, presented=["a", "b"])
    ok, detail = loop_health.check_evidence_consistent()
    assert not ok and "disagree" in detail


def test_agreement_passes(sandbox):
    (sandbox / "evidence.json").write_text(json.dumps(
        {"endpoints": {"a": {"tier": "validated"}, "b": {"tier": "weak"}}}))
    _history(loop_health.HISTORY, days_ago=1, presented=["a"])
    assert loop_health.check_evidence_consistent()[0]


def test_nothing_validated_needs_a_human(sandbox):
    (sandbox / "evidence.json").write_text(json.dumps(
        {"endpoints": {"a": {"tier": "weak"}}}))
    ok, detail = loop_health.check_evidence_consistent()
    assert not ok and "presents nothing" in detail


def test_the_exit_status_is_what_a_scheduled_job_branches_on(sandbox):
    """Unhealthy must exit 1, or the workflow cannot tell the difference."""
    _history(loop_health.HISTORY, days_ago=99, presented=["a"])
    assert loop_health.main([]) == 1
    written = json.loads((sandbox / "health.json").read_text())
    assert written["healthy"] is False
    assert "history_fresh" in written["failed"]
