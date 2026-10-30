"""A row nobody assayed is not a negative result.

The trainer took `Y[endpoint].fillna(0)`, so every row that was never
measured for an endpoint became a row asserting that endpoint did not occur.
Thrombocytopenia was trained on 4 real positives and 2,030 compounds nobody
had assayed for platelets, and because those rows came along with whatever
else was in the corpus, adding an unrelated source moved its score.

That is how a source of liver measurements promoted two endpoints that have
no liver data: kupffer_cell_activation went 0.478 -> 0.606 PR-AUC and
injection_site_reaction 0.209 -> 0.317 purely from 444 new rows that are
blank for both.

These pin the rule across the trainer and every evaluator, because a model
scored on real labels against a baseline scored on invented ones is worse
than either alone.
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parent.parent
PROCESSED = ROOT / "data" / "processed"
RESULTS = ROOT / "reports" / "results"


@pytest.fixture(scope="module")
def labels():
    return pd.read_csv(PROCESSED / "label_matrix.csv")


@pytest.fixture(scope="module")
def summary():
    return json.loads((RESULTS / "training_summary_after.json").read_text())


def test_each_endpoint_was_trained_only_on_its_measured_rows(labels, summary):
    """`samples` must be the labelled count, not the whole corpus."""
    for endpoint, spec in summary.items():
        if not isinstance(spec, dict) or spec.get("samples") is None:
            continue
        if endpoint not in labels.columns:
            continue  # trained externally, on its own matrix
        measured = int(pd.to_numeric(labels[endpoint], errors="coerce").notna().sum())
        assert spec["samples"] == measured, (
            f"{endpoint} was trained on {spec['samples']} rows but only "
            f"{measured} carry its label; the rest were invented negatives"
        )


def test_no_endpoint_claims_the_whole_corpus(labels, summary):
    """The signature of the bug: every endpoint reporting one sample count."""
    counts = {s["samples"] for s in summary.values()
              if isinstance(s, dict) and s.get("samples")}
    assert len(counts) > 1, (
        "every endpoint reports the same number of samples, which means they "
        f"are all being trained on the full corpus: {counts}"
    )


def test_the_zero_fill_is_gone_from_the_trainer_and_the_evaluators():
    """One of these scoring on invented negatives is enough to break the gate."""
    offenders = []
    for path in [ROOT / "src/oligotox/training/train.py",
                 ROOT / "scripts/evaluate_baselines.py",
                 ROOT / "scripts/evaluate_calibration.py",
                 ROOT / "scripts/endpoint_requirements.py"]:
        for line in path.read_text().splitlines():
            stripped = line.strip()
            if stripped.startswith("#"):
                continue
            if 'errors="coerce").fillna(0)' in stripped:
                offenders.append(f"{path.name}: {stripped}")
    assert not offenders, offenders


def test_a_sparse_endpoint_is_not_padded_to_look_well_evidenced(labels, summary):
    """Thrombocytopenia has 28 measured rows and must report 28."""
    if "thrombocytopenia" not in labels.columns:
        pytest.skip("endpoint not in this build")
    measured = pd.to_numeric(labels["thrombocytopenia"], errors="coerce")
    assert measured.notna().sum() < 100, "fixture assumption changed"
    spec = summary.get("thrombocytopenia", {})
    if spec.get("samples") is not None:
        assert spec["samples"] == int(measured.notna().sum())
