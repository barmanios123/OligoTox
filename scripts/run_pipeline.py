"""The continuous improvement loop, as one command.

Point it at the repository after new data lands and it re-does everything
that decides what the platform claims:

    python scripts/run_pipeline.py

    normalize    declared sources to labelled tables, directions enforced
    train        12 internal endpoints; externally trained ones are skipped
    baselines    source-only, GC-only, length-only, with bootstrap CIs
    calibration  Brier skill against the base rate, per endpoint
    tiers        recompute every endpoint's tier from the evidence above
    transitions  compare against the last run and record what changed
    requirements what each failing endpoint still needs

Nothing here can promote an endpoint on its own. Tiers come from
`derive_evidence_tiers.py`, which reads the measurements and applies the same
gates every time. This script's job is to make sure those measurements are
current and that any change is visible.

Options:

    --skip-slow      leave out the baseline step, which takes about 95 s
    --check          run everything and fail if any endpoint changed tier,
                     which is what you want in CI
    --trigger TEXT   what caused this run, recorded in the history log
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from oligotox.lifecycle import (  # noqa: E402
    EndpointState,
    apply_transitions,
    describe_change,
    record_history,
)
from oligotox.lifecycle.transitions import previous_from_history  # noqa: E402

CONFIG = ROOT / "artifacts" / "config" / "endpoint_evidence.json"
HISTORY = ROOT / "reports" / "results" / "endpoint_history.jsonl"

# (label, argv, slow). Order matters: each step consumes the previous output.
STEPS = [
    ("normalize external data", [sys.executable, "scripts/normalize_external_tox.py"], False),
    ("train internal endpoints", [sys.executable, "-m", "oligotox.training.train"], False),
    ("publish training summary", None, False),  # handled in-process
    ("baselines", [sys.executable, "scripts/evaluate_baselines.py"], True),
    ("calibration", [sys.executable, "scripts/evaluate_calibration.py"], True),
    # The gate reads this, so it has to be current. It used to be an
    # optional script run by hand, which is fine for a report and not fine
    # for an input to a tier.
    ("paired baselines", [sys.executable, "scripts/evaluate_paired_baselines.py"], True),
    ("evidence tiers", [sys.executable, "scripts/derive_evidence_tiers.py"], False),
    ("endpoint requirements", [sys.executable, "scripts/endpoint_requirements.py"], False),
]


def _current_states() -> dict[str, EndpointState]:
    endpoints = json.loads(CONFIG.read_text())["endpoints"]
    return {name: EndpointState.from_evidence(name, info) for name, info in endpoints.items()}


def _publish_training_summary() -> None:
    """Training writes to artifacts/metrics; the evaluation reads reports/results."""
    source = ROOT / "artifacts" / "metrics" / "training_summary.json"
    target = ROOT / "reports" / "results" / "training_summary_after.json"
    if source.exists():
        target.write_text(source.read_text())


def run_step(label: str, argv: list[str] | None) -> bool:
    print(f"\n=== {label}")
    start = time.monotonic()
    if argv is None:
        _publish_training_summary()
        print(f"    ok ({time.monotonic() - start:.1f}s)")
        return True
    result = subprocess.run(argv, cwd=ROOT, capture_output=True, text=True,
                            env={**__import__("os").environ,
                                 "PYTHONPATH": str(ROOT / "src")})
    elapsed = time.monotonic() - start
    if result.returncode != 0:
        print(f"    FAILED ({elapsed:.1f}s)")
        print((result.stderr or result.stdout)[-1500:])
        return False
    tail = [ln for ln in (result.stdout or "").splitlines() if ln.strip()][-3:]
    for line in tail:
        print(f"    {line[:160]}")
    print(f"    ok ({elapsed:.1f}s)")
    return True


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--skip-slow", action="store_true")
    parser.add_argument("--check", action="store_true",
                        help="fail if any endpoint changed tier")
    parser.add_argument("--trigger", default="manual run")
    args = parser.parse_args(argv)

    previous = previous_from_history(HISTORY) or _current_states()

    for label, step_argv, slow in STEPS:
        if slow and args.skip_slow:
            print(f"\n=== {label}\n    skipped")
            continue
        if not run_step(label, step_argv):
            print(f"\npipeline stopped at: {label}")
            return 1

    current = _current_states()
    transitions = apply_transitions(previous, current)
    changed = [t for t in transitions if t.direction != "unchanged"]

    print("\n=== endpoint lifecycle")
    for transition in transitions:
        marker = "  " if transition.direction == "unchanged" else "* "
        print(f"{marker}{describe_change(transition)}")

    entry = record_history(HISTORY, transitions, trigger=args.trigger)
    presented = entry["presented"]
    print(f"\npresented as predictions: {len(presented)} of {len(current)} "
          f"({', '.join(presented) or 'none'})")
    if entry["promoted"]:
        print(f"promoted: {', '.join(entry['promoted'])}")
    if entry["demoted"]:
        print(f"demoted:  {', '.join(entry['demoted'])}")
    print(f"recorded in {HISTORY.relative_to(ROOT)}")

    if args.check and changed:
        print(f"\n--check: {len(changed)} endpoint(s) changed tier since the last run")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
