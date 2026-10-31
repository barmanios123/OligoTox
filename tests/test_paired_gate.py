"""The simple-baseline check is a test of a difference, and fails closed.

It used to compare the model's marginal 95% interval lower bound to the
baseline's point estimate. Both are scored on the same rows with the same
folds, so their errors are correlated, and that comparison is not a test of
whether one beats the other: it is much more conservative than one.

The replacement asks directly, by resampling whole sequences and scoring
both models on each resample. It is not a looser gate. On the run that
introduced it, it promoted hepatocellular_injury (+0.119 [+0.036, +0.194]
against GC content) and demoted offtarget_gene_silencing (+0.355, but
[-0.075, +0.688] on 12 positive sequences).

The part worth guarding is what happens when the paired result is missing.
A gate that passes an endpoint because its evidence was not computed is
worse than the gate it replaced.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

import derive_evidence_tiers as tiers  # noqa: E402

GOOD = {"pr_auc": 0.60, "ci_low": 0.40, "ci_high": 0.75}


def _classify(endpoint="hepatocellular_injury"):
    return tiers.classify(GOOD, source_pr_auc=0.20, best_simple=0.45,
                          base_rate=0.30, n_sequences=100,
                          within_beats_base=True, endpoint=endpoint)


def test_a_missing_paired_result_fails_the_check(monkeypatch):
    """Absent evidence is not a pass."""
    monkeypatch.setattr(tiers, "_paired", lambda endpoint: None)
    tier, reason, checks = _classify()
    assert checks["beats_simple_baseline_paired"] is False
    assert tier == "weak"
    assert "no paired comparison" in reason


def test_a_difference_that_spans_zero_fails(monkeypatch):
    """+0.355 [-0.075, +0.688] is what off-target gene silencing looks like."""
    monkeypatch.setattr(tiers, "_paired", lambda endpoint: {
        "difference": 0.355, "ci": [-0.075, 0.688],
        "p_better": 0.917, "clears_paired": False})
    tier, reason, checks = _classify()
    assert checks["beats_simple_baseline_paired"] is False
    assert tier == "weak"
    assert "does not exclude zero" in reason
    assert "P(better) 0.917" in reason


def test_a_difference_that_excludes_zero_passes(monkeypatch):
    """+0.119 [+0.036, +0.194] is hepatocellular injury against GC content."""
    monkeypatch.setattr(tiers, "_paired", lambda endpoint: {
        "difference": 0.119, "ci": [0.036, 0.194],
        "p_better": 0.996, "clears_paired": True})
    tier, reason, checks = _classify()
    assert checks["beats_simple_baseline_paired"] is True
    assert tier == "validated"
    assert "paired" in reason


def test_the_marginal_check_is_gone(monkeypatch):
    """The old key must not linger, or two gates disagree in one config."""
    monkeypatch.setattr(tiers, "_paired", lambda endpoint: {
        "difference": 0.119, "ci": [0.036, 0.194],
        "p_better": 0.996, "clears_paired": True})
    _, _, checks = _classify()
    assert "ci_clears_simple_baselines" not in checks


def test_an_earlier_failure_still_wins(monkeypatch):
    """Source confounding is decided before the baseline comparison."""
    monkeypatch.setattr(tiers, "_paired", lambda endpoint: {
        "difference": 0.5, "ci": [0.3, 0.7],
        "p_better": 1.0, "clears_paired": True})
    tier, _, _ = tiers.classify(GOOD, source_pr_auc=0.90, best_simple=0.45,
                                base_rate=0.30, n_sequences=100,
                                within_beats_base=True, endpoint="x")
    assert tier == "source_confounded"


def test_the_gate_runs_the_paired_step_itself():
    """It cannot depend on somebody remembering to run a script."""
    source = (ROOT / "scripts" / "run_pipeline.py").read_text()
    # The STEPS list, not the module docstring, which names both scripts.
    start = source.index("STEPS = [")
    steps = source[start:source.index("\n]", start)]
    assert "evaluate_paired_baselines.py" in steps, (
        "the gate's evidence is not computed by the pipeline that uses it"
    )
    assert steps.index("evaluate_paired_baselines.py") < \
        steps.index("derive_evidence_tiers.py"), "it must run before tiers"


def test_the_shipped_config_records_the_paired_check():
    config = json.loads(
        (ROOT / "artifacts/config/endpoint_evidence.json").read_text())
    # cns_toxicity is trained and gated externally, on its own checks.
    scored = [v for v in config["endpoints"].values()
              if "beats_source_baseline" in (v.get("checks") or {})]
    if not scored:
        pytest.skip("no endpoint goes through the internal gate in this build")
    for spec in scored:
        assert "beats_simple_baseline_paired" in spec["checks"]
        assert "ci_clears_simple_baselines" not in spec["checks"]
