"""Compare the tiers a run produced against the tiers that were shipped.

Every pipeline run recomputes each endpoint's tier from the current evidence.
This turns the before and after into explicit transitions, so a promotion or
a demotion is a recorded event with a reason rather than a silent diff in a
JSON file.

Two properties the rest of the system depends on:

  promotion is never implicit   an endpoint reaching `validated` is recorded
                                with the evidence that got it there
  demotion is automatic         an endpoint that stops clearing its gates
                                loses its tier on the next run, with no
                                grace period and no manual override

There is deliberately no mechanism to pin or force a tier. If an endpoint
should be presented as a prediction, the way to do that is to give it
evidence.
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

# Ordered best to worst. Index comparison decides promotion from demotion.
TIER_RANK = {
    "validated": 0,
    "weak": 1,
    "source_confounded": 2,
    "insufficient_data": 3,
    "not_evaluated": 4,
}

PRESENTED_AS_PREDICTION = "validated"


@dataclass
class EndpointState:
    """What one endpoint looked like at one point in time."""

    endpoint: str
    tier: str
    positive_sequences: int | None = None
    pr_auc: float | None = None
    brier_skill: float | None = None
    probabilities_calibrated: bool | None = None
    reason: str = ""

    @classmethod
    def from_evidence(cls, endpoint: str, info: dict) -> "EndpointState":
        return cls(
            endpoint=endpoint,
            tier=info.get("tier", "not_evaluated"),
            positive_sequences=info.get("positive_sequences"),
            pr_auc=info.get("pr_auc"),
            brier_skill=info.get("brier_skill"),
            probabilities_calibrated=info.get("probabilities_calibrated"),
            reason=info.get("reason", ""),
        )


@dataclass
class Transition:
    endpoint: str
    before: str
    after: str
    direction: str  # "promoted", "demoted", "unchanged", "new"
    now_presented: bool
    was_presented: bool
    reason: str = ""
    evidence: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict:
        return asdict(self)


def _direction(before: str, after: str) -> str:
    if before == after:
        return "unchanged"
    rank_before = TIER_RANK.get(before, TIER_RANK["not_evaluated"])
    rank_after = TIER_RANK.get(after, TIER_RANK["not_evaluated"])
    return "promoted" if rank_after < rank_before else "demoted"


def apply_transitions(
    previous: dict[str, EndpointState],
    current: dict[str, EndpointState],
) -> list[Transition]:
    """Every endpoint in `current`, with how it changed.

    An endpoint absent from `previous` is "new" rather than promoted: it has
    no history to have been promoted from.
    """
    transitions = []
    for endpoint in sorted(current):
        now = current[endpoint]
        before = previous.get(endpoint)
        if before is None:
            direction, before_tier = "new", "not_evaluated"
        else:
            before_tier = before.tier
            direction = _direction(before_tier, now.tier)

        transitions.append(
            Transition(
                endpoint=endpoint,
                before=before_tier,
                after=now.tier,
                direction=direction,
                now_presented=now.tier == PRESENTED_AS_PREDICTION,
                was_presented=before_tier == PRESENTED_AS_PREDICTION,
                reason=now.reason,
                evidence={
                    "positive_sequences": now.positive_sequences,
                    "pr_auc": now.pr_auc,
                    "brier_skill": now.brier_skill,
                    "probabilities_calibrated": now.probabilities_calibrated,
                },
            )
        )

    # An endpoint that disappeared entirely stops being presented, which is a
    # demotion and has to be visible.
    for endpoint in sorted(set(previous) - set(current)):
        before = previous[endpoint]
        transitions.append(
            Transition(
                endpoint=endpoint,
                before=before.tier,
                after="not_evaluated",
                direction="demoted",
                now_presented=False,
                was_presented=before.tier == PRESENTED_AS_PREDICTION,
                reason="no longer produced by the pipeline",
            )
        )
    return transitions


def describe_change(transition: Transition) -> str:
    """One line a human can read in a log."""
    arrow = f"{transition.before} -> {transition.after}"
    if transition.direction == "unchanged":
        return f"{transition.endpoint}: {transition.after} (unchanged)"
    if transition.direction == "new":
        return f"{transition.endpoint}: new, {transition.after}"
    verb = "PROMOTED" if transition.direction == "promoted" else "DEMOTED"
    note = ""
    if transition.now_presented and not transition.was_presented:
        note = " and is now presented as a prediction"
    elif transition.was_presented and not transition.now_presented:
        note = " and is no longer presented as a prediction"
    return f"{transition.endpoint}: {verb} {arrow}{note}"


def load_history(path: str | Path) -> list[dict]:
    """Every run recorded so far, oldest first."""
    path = Path(path)
    if not path.exists():
        return []
    entries = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if line:
            entries.append(json.loads(line))
    return entries


def record_history(
    path: str | Path,
    transitions: list[Transition],
    trigger: str = "pipeline run",
) -> dict:
    """Append this run to the log. Append-only, never rewritten."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    entry = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "trigger": trigger,
        "presented": sorted(t.endpoint for t in transitions if t.now_presented),
        "promoted": sorted(t.endpoint for t in transitions if t.direction == "promoted"),
        "demoted": sorted(t.endpoint for t in transitions if t.direction == "demoted"),
        "transitions": [t.to_dict() for t in transitions if t.direction != "unchanged"],
    }
    with path.open("a") as handle:
        handle.write(json.dumps(entry) + "\n")
    return entry


def previous_from_history(path: str | Path) -> dict[str, EndpointState]:
    """Reconstruct the last known tiers from the log.

    Only changes are recorded, so the current state is the latest tier seen
    for each endpoint across the whole log.
    """
    state: dict[str, EndpointState] = {}
    for entry in load_history(path):
        for raw in entry.get("transitions", []):
            state[raw["endpoint"]] = EndpointState(
                endpoint=raw["endpoint"],
                tier=raw["after"],
                reason=raw.get("reason", ""),
                **{k: v for k, v in (raw.get("evidence") or {}).items()
                   if k in {"positive_sequences", "pr_auc", "brier_skill",
                            "probabilities_calibrated"}},
            )
    return state
