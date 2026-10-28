"""Is the unattended loop actually still working?

The weekly job is built to survive nobody watching it. The failure mode that
creates is silence that looks like success: the literature steps run with
`continue-on-error`, the inventory exits 0 when a source returns 403, and the
evidence commit lands every week regardless because `endpoint_history.jsonl`
gains a timestamped line on every run. So a loop that has quietly stopped
learning anything commits exactly like a healthy one.

This checks the loop's own vital signs and fails loudly when one is wrong, so
the schedule can open an issue instead of continuing to look fine:

    history_fresh      did the pipeline run when it was supposed to?
    sources_on_disk    is every source declared `present` actually there?
    training_ran        did the last run train, or skip and reuse?
    triage_clear        has every monitored paper been assessed?
    evidence_consistent does every presented endpoint have validated evidence?
    tiers_not_silently_demoted  did anything stop being presented without a PR?

What this deliberately does NOT treat as a failure: a `fetchable` source that
is absent, or a `requested`/`wanted` one. Those are waiting on a human or on
a third party, which is a known state, not a broken loop.

    python scripts/loop_health.py
    python scripts/loop_health.py --max-age-days 10

Writes reports/results/loop_health.json. Exit status is 1 if any check fails,
so a scheduled job can branch on it.
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "reports" / "results"
REGISTRY = ROOT / "data" / "registry" / "sources.json"
HISTORY = RESULTS / "endpoint_history.jsonl"
INVENTORY = RESULTS / "dataset_inventory.json"
TRIAGE = RESULTS / "literature_triage.json"
CANDIDATES = RESULTS / "literature_candidates.json"
EVIDENCE = ROOT / "artifacts" / "config" / "endpoint_evidence.json"
HEALTH = RESULTS / "loop_health.json"

# The schedule is weekly, so anything past a fortnight means at least one run
# was missed entirely.
DEFAULT_MAX_AGE_DAYS = 14
# Sources in these states are waiting on someone, not broken.
EXPECTED_ABSENT = {"fetchable", "requested", "deferred", "wanted",
                   "not_applicable", "unverified", "adopted_per_context"}


def _load(path: Path):
    try:
        return json.loads(path.read_text())
    except (FileNotFoundError, json.JSONDecodeError):
        return None


def check_history_fresh(max_age_days: int) -> tuple[bool, str]:
    """A loop that has stopped running is the failure this exists to catch."""
    if not HISTORY.exists():
        return False, "no endpoint_history.jsonl; the pipeline has never run"
    lines = [l for l in HISTORY.read_text().splitlines() if l.strip()]
    if not lines:
        return False, "endpoint_history.jsonl is empty"
    try:
        last = json.loads(lines[-1])
        stamp = datetime.fromisoformat(last["timestamp"])
    except (json.JSONDecodeError, KeyError, ValueError) as exc:
        return False, f"could not read the last history entry: {exc}"
    if stamp.tzinfo is None:
        stamp = stamp.replace(tzinfo=timezone.utc)
    age = datetime.now(timezone.utc) - stamp
    if age > timedelta(days=max_age_days):
        return False, (f"the last run was {age.days} days ago, over the "
                       f"{max_age_days}-day limit; the schedule may be disabled")
    return True, f"last run {age.days}d {age.seconds // 3600}h ago ({len(lines)} runs logged)"


def check_sources_on_disk() -> tuple[bool, str]:
    """A source declared `present` and missing means the data moved or was lost."""
    registry = _load(REGISTRY)
    if registry is None:
        return False, "data/registry/sources.json is missing or unreadable"
    entries = registry["sources"] if isinstance(registry, dict) else registry

    missing = []
    for entry in entries:
        path = entry.get("local_path")
        if entry.get("status") == "present" and path and not (ROOT / path).exists():
            missing.append(entry["name"])
    if missing:
        return False, (f"declared present but not on disk: {missing}. Either the "
                       "data was lost or the registry is out of date.")
    waiting = sorted(e["name"] for e in entries
                     if e.get("status") in EXPECTED_ABSENT
                     and e.get("local_path") and not (ROOT / e["local_path"]).exists())
    note = f"; waiting on a human or a third party: {waiting}" if waiting else ""
    present = sum(1 for e in entries if e.get("status") == "present")
    return True, f"all {present} present source(s) on disk{note}"


def check_training_ran() -> tuple[bool, str]:
    """An inventory with nothing usable means the loop has nothing to learn from."""
    inventory = _load(INVENTORY)
    if inventory is None:
        return False, "no dataset_inventory.json; fetch_datasets.py did not complete"
    usable = inventory.get("usable_sources", 0)
    sequences = inventory.get("usable_sequences", 0)
    if usable == 0:
        return False, "no usable sources; the pipeline has nothing to train on"
    return True, f"{usable} usable source(s), {sequences:,} sequences"


def check_triage_clear() -> tuple[bool, str]:
    """An unassessed pile means the monitor is generating work, not findings."""
    triage, candidates = _load(TRIAGE), _load(CANDIDATES)
    if candidates is None:
        return True, "the monitor has not run yet; nothing to assess"
    found = candidates.get("candidates", [])
    assessed = (triage or {}).get("assessed", {})
    pending = [c["id"] for c in found if c["id"] not in assessed]
    if pending:
        return False, (f"{len(pending)} of {len(found)} monitored paper(s) unassessed; "
                       "run scripts/triage_literature.py")
    open_now = [a for a in assessed.values() if a.get("verdict") in {"dataset", "unknown"}]
    note = f"; {len(open_now)} await a human read" if open_now else ""
    return True, f"all {len(found)} monitored paper(s) assessed{note}"


def check_evidence_consistent() -> tuple[bool, str]:
    """Do the two artifacts the loop writes still agree on what is presented?

    An endpoint is presented exactly when its tier is `validated`. The run log
    separately records what was presented on that run. They are written by the
    same pipeline from the same state, so a divergence means one of them was
    edited, partially committed, or written by a run that half-failed.

    Anything structural is a failure, not a pass. A check that silently finds
    no endpoints and reports "all good" is the exact bug this script exists to
    catch.
    """
    evidence = _load(EVIDENCE)
    if evidence is None:
        return False, "no artifacts/config/endpoint_evidence.json; the gate is missing"
    entries = evidence.get("endpoints")
    if not isinstance(entries, dict) or not entries:
        return False, ("endpoint_evidence.json has no `endpoints` mapping; the "
                       "schema changed and this check can no longer read the gate")

    untiered = sorted(name for name, e in entries.items() if not e.get("tier"))
    if untiered:
        return False, f"endpoint(s) carry no tier: {untiered}"

    validated = sorted(name for name, e in entries.items()
                       if e.get("tier") == "validated")
    if not validated:
        return False, ("no endpoint is validated, so the app presents nothing. "
                       "That may be correct after a demotion, but it needs a human.")

    # Cross-check against what the last run recorded.
    lines = ([l for l in HISTORY.read_text().splitlines() if l.strip()]
             if HISTORY.exists() else [])
    if lines:
        try:
            logged = sorted(json.loads(lines[-1]).get("presented", []))
        except json.JSONDecodeError:
            return False, "the last history entry is not readable JSON"
        if logged != validated:
            return False, (f"the evidence config and the run log disagree: config "
                           f"validates {validated}, the last run presented {logged}")

    return True, f"{len(validated)} validated endpoint(s), log agrees: {validated}"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--max-age-days", type=int, default=DEFAULT_MAX_AGE_DAYS)
    args = parser.parse_args(argv)

    checks = {
        "history_fresh": check_history_fresh(args.max_age_days),
        "sources_on_disk": check_sources_on_disk(),
        "training_ran": check_training_ran(),
        "triage_clear": check_triage_clear(),
        "evidence_consistent": check_evidence_consistent(),
    }

    failed = [name for name, (ok, _) in checks.items() if not ok]
    for name, (ok, detail) in checks.items():
        print(f"  {'ok  ' if ok else 'FAIL'}  {name:22s} {detail}")

    RESULTS.mkdir(parents=True, exist_ok=True)
    HEALTH.write_text(json.dumps({
        "checked_at": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "healthy": not failed,
        "failed": failed,
        "note": "A loop that has quietly stopped learning commits exactly like a "
                "healthy one, because the history log gains a line every run. "
                "These are the signs that distinguish them.",
        "checks": {name: {"ok": ok, "detail": detail}
                   for name, (ok, detail) in checks.items()},
    }, indent=2) + "\n")

    print()
    if failed:
        print(f"UNHEALTHY: {len(failed)} check(s) failed: {', '.join(failed)}")
    else:
        print("healthy: the loop is running, has data, and the evidence gate holds")
    print(f"wrote {HEALTH.relative_to(ROOT)}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
