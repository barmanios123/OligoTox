"""Endpoint lifecycle: what an endpoint is allowed to become, and when.

This exists so the platform can take in new data continuously without the
failure that started this project. An endpoint's tier is never set by hand
and never carried forward out of habit. It is recomputed from evidence on
every run, and every change is recorded with the reason.

The rule is a ratchet. Evidence can promote an endpoint and evidence can
demote it, but the bar itself never moves: promotion requires passing every
gate in `oligotox.audit`, and a gate is never waived because an endpoint used
to pass. An endpoint that stops clearing its baselines is demoted on the next
run whether or not anyone is watching.

    from oligotox.lifecycle import apply_transitions, load_history
"""
from .transitions import (
    EndpointState,
    Transition,
    apply_transitions,
    describe_change,
    load_history,
    record_history,
)

__all__ = [
    "EndpointState",
    "Transition",
    "apply_transitions",
    "describe_change",
    "load_history",
    "record_history",
]
