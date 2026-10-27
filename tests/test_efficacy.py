"""The Efficacy tab was removed once for showing untrained stubs.

These tests keep what replaced it honest: a model only scores the modality it
was trained on, a `validated` tier has to be backed by cross-study transfer,
and the app's feature builder has to match the one used in training.
"""
import json
import sys
from pathlib import Path

import pytest

from oligotox.efficacy.potency import guide_features, load_efficacy_config, score_efficacy

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
from evaluate_efficacy import guide_features as training_guide_features  # noqa: E402

SIRNA = "AAAACAAAGAAGGCCACCATT"


@pytest.fixture
def config(root):
    return load_efficacy_config(root / "artifacts/config")


@pytest.fixture
def results(root):
    return json.loads((root / "reports/results/efficacy.json").read_text())["results"]


def test_app_and_training_features_agree():
    """A drift here would score candidates with the wrong columns, silently."""
    assert guide_features(SIRNA) == training_guide_features(SIRNA)


def test_a_model_only_scores_its_own_modality(root):
    out = score_efficacy(SIRNA, "PMO", root / "artifacts/models", root / "artifacts/config")
    assert out["available"] is False
    assert "PMO" in out["reason"]

    sirna = score_efficacy(SIRNA, "siRNA", root / "artifacts/models", root / "artifacts/config")
    assert set(sirna["predictions"]) == {"sirna_potency"}


def test_every_shipped_model_states_units_and_evidence(config):
    assert config, "no efficacy models configured"
    for name, spec in config.items():
        assert spec["units"] and spec["evidence"], name
        assert spec["tier"] in {"validated", "weak"}, name


def test_validated_tier_requires_cross_study_transfer(config, results):
    transfers = [r for r in results if r["evaluation"] == "cross-study transfer"
                 and r["model"] == "All guide features"]
    for name, spec in config.items():
        if spec["tier"] != "validated":
            continue
        assert transfers, f"{name} claims validated with no transfer result"
        for r in transfers:
            assert r["ci_low"] > 0, f"{name}: transfer CI reaches zero"


def test_sequence_model_beats_the_gene_mean_baseline(results):
    """Guards against predicting the target rather than the oligo."""
    for dataset in ("Ichihara_2007_1", "McQuisten_2007"):
        rows = {r["model"]: r for r in results
                if r["dataset"] == dataset and r["evaluation"] == "gene-grouped CV"}
        assert rows["All guide features"]["spearman"] > rows["Target gene mean (no sequence)"]["spearman"]
        assert rows["All guide features"]["ci_low"] > rows["GC only"]["spearman"], dataset


def test_scores_are_reported_in_training_units(root):
    out = score_efficacy(SIRNA, "siRNA", root / "artifacts/models", root / "artifacts/config")
    pred = out["predictions"]["sirna_potency"]
    # Percent knockdown, not a probability.
    assert "percent" in pred["units"]
    assert 0.0 <= pred["score"] <= 150.0
