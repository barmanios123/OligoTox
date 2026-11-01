"""A tier decided inside the bootstrap's own noise must say so.

hepatocellular_injury was demoted on a paired CI lower bound of -0.001. At
1,500 bootstrap draws that bound has a Monte Carlo standard deviation of
0.0032 and ran from -0.0023 to +0.0110 across twelve seeds: the gate passed
on eleven of them and the shipped run had drawn the twelfth. The criterion
was fine; the estimate of it was noisier than the margin being decided.

Raising the draws to 6,000 narrows that to +0.0028..+0.0077 and no other
endpoint's verdict moves. But precision is not certainty, so a bound still
within Monte Carlo error of zero is now labelled as unsettled rather than
reported as a result.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

import derive_evidence_tiers as tiers  # noqa: E402
import evaluate_paired_baselines as paired  # noqa: E402


def test_the_bootstrap_is_precise_enough_for_the_margins_it_decides():
    assert paired.N_BOOT >= 6000, (
        "1,500 draws put more noise on the lower bound than the margins "
        "this gate turns on"
    )


def test_a_bound_inside_the_noise_is_called_unsettled():
    spec = {"difference": 0.108, "ci": [0.003, 0.205], "p_better": 0.979,
            "knife_edge": True}
    sentence = tiers._paired_sentence(spec, 0.38)
    assert "not stable across resampling seeds" in sentence


def test_a_bound_clear_of_the_noise_is_not_hedged():
    spec = {"difference": 0.53, "ci": [0.229, 0.784], "p_better": 0.999,
            "knife_edge": False}
    sentence = tiers._paired_sentence(spec, 0.28)
    assert "not stable" not in sentence


def test_the_shipped_report_records_its_own_monte_carlo_error():
    path = ROOT / "reports/results/paired_baselines.json"
    if not path.exists():
        pytest.skip("paired baselines not computed in this build")
    for name, spec in json.loads(path.read_text())["endpoints"].items():
        assert "ci_low_monte_carlo_sd" in spec, name
        assert "knife_edge" in spec, name
        # The flag must agree with the numbers it is derived from.
        expected = abs(spec["ci"][0]) < paired.KNIFE_EDGE_SDS * spec["ci_low_monte_carlo_sd"]
        assert spec["knife_edge"] == expected, name


def test_a_knife_edge_endpoint_says_so_wherever_its_tier_is_explained():
    config = json.loads(
        (ROOT / "artifacts/config/endpoint_evidence.json").read_text())
    report = ROOT / "reports/results/paired_baselines.json"
    if not report.exists():
        pytest.skip("paired baselines not computed in this build")
    paired_by_endpoint = json.loads(report.read_text())["endpoints"]
    for name, spec in config["endpoints"].items():
        info = paired_by_endpoint.get(name)
        if not info or not info.get("knife_edge"):
            continue
        if "paired difference" not in spec.get("reason", ""):
            continue  # an earlier check decided this tier
        assert "not stable across resampling seeds" in spec["reason"], name
