"""One endpoint's confidence interval must not depend on another's data.

`evaluate_baselines.py` created a single `default_rng(42)` before the loop
and passed it to every endpoint in turn, so each endpoint's bootstrap drew
from wherever the previous ones had left the stream. Adding 1,008
hepatocellular rows moved `offtarget_gene_silencing`'s lower bound from
0.209 to 0.188 without touching its model, its data or its score, and
demoted it, because that bound sits 0.004 from the number it is compared
against.

A tier is a claim about an endpoint. It must not be a claim about what else
happened to be in the corpus that day.
"""
from __future__ import annotations

import sys
import zlib
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))


def _seeded(name: str) -> np.random.Generator:
    return np.random.default_rng([42, zlib.crc32(name.encode())])


def test_each_endpoint_draws_its_own_stream():
    """Two endpoints must not share a sequence of draws."""
    a, b = _seeded("hepatocellular_injury"), _seeded("offtarget_gene_silencing")
    assert not np.array_equal(a.random(20), b.random(20))


def test_the_stream_is_the_same_whatever_ran_before_it():
    """The whole point: reproducible per endpoint, order-independent."""
    first = _seeded("offtarget_gene_silencing").random(20)
    noise = _seeded("hepatocellular_injury")
    noise.random(5000)
    second = _seeded("offtarget_gene_silencing").random(20)
    assert np.array_equal(first, second)


def test_the_evaluators_seed_inside_their_loops():
    """A generator built before the loop is the defect itself."""
    for name in ("evaluate_baselines.py", "evaluate_motif_replication.py",
                 "evaluate_hepatotox_model.py"):
        source = (ROOT / "scripts" / name).read_text()
        assert "zlib.crc32" in source, (
            f"{name} does not seed its bootstrap per unit of comparison"
        )


def test_a_failed_baseline_check_reports_the_paired_numbers():
    """The reason must carry the test the decision was actually made on.

    This check previously forbade the phrase "not distinguishable", because
    the gate reached it by comparing a marginal interval to a point
    estimate, which is not a test of a difference. The gate now runs the
    paired test, so the phrase is earned; what matters is that the numbers
    behind it are in the reason.
    """
    import json

    evidence = json.loads(
        (ROOT / "artifacts/config/endpoint_evidence.json").read_text())
    for name, spec in evidence["endpoints"].items():
        checks = spec.get("checks", {})
        if checks.get("beats_simple_baseline_paired") is not False:
            continue
        # Only where this check is the one that decided the tier. An
        # endpoint stopped earlier, by too few sequences or by source
        # confounding, never reaches the baseline comparison.
        if not (checks.get("has_enough_sequences")
                and checks.get("beats_source_baseline")):
            continue
        assert "paired difference" in spec["reason"], name
        assert "P(better)" in spec["reason"], name
