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


def test_context_models_name_the_assay_they_predict(root):
    """A per-assay model must never read as general potency."""
    from oligotox.efficacy.potency import load_context_models

    config = load_context_models(root / "artifacts/config")
    contexts = config.get("contexts", {})
    if not contexts:
        pytest.skip("no context models trained in this checkout")

    for name, spec in contexts.items():
        assert spec["context"]["cell_line"], name
        # The units say where the number applies.
        assert spec["context"]["cell_line"] in spec["units"], name
        # And the caveat says what it is not.
        assert "general potency" in spec["caveat"].lower(), name


def test_every_shipped_context_beat_its_own_baselines(root):
    """The gate is per context, not an average over contexts."""
    import csv

    path = root / "reports/results/context_potency.csv"
    if not path.exists():
        pytest.skip("context potency not evaluated in this checkout")

    with path.open() as fh:
        rows = list(csv.DictReader(fh))
    assert rows, "evaluation produced no contexts"
    for row in rows:
        if row["ships"].lower() != "true":
            continue
        rho = float(row["spearman"])
        assert rho > float(row["gc_only"]), row["context"]
        assert rho > float(row["gene_mean"]), row["context"]
        assert rho >= 0.15, row["context"]


def test_a_context_that_failed_is_not_shipped(root):
    """Something must actually have been rejected, or the gate is decorative."""
    import csv
    import json

    path = root / "reports/results/context_potency.csv"
    if not path.exists():
        pytest.skip("context potency not evaluated in this checkout")
    with path.open() as fh:
        rows = list(csv.DictReader(fh))

    rejected = [r for r in rows if r["ships"].lower() != "true"]
    assert rejected, "no context was rejected; check the gate still bites"
    shipped = json.loads(
        (root / "artifacts/config/potency_contexts.json").read_text()
    )["contexts"]
    for row in rejected:
        assert row["context"] not in shipped, row["context"]
        assert row["blocked_by"], row["context"]
