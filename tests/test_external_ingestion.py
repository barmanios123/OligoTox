"""Guards on the external ingestion, where a label direction was inverted."""
import json
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
from normalize_external_tox import HIGHER_IS_TOXIC, LOWER_IS_TOXIC, SOURCES  # noqa: E402


def test_every_source_declares_a_direction_or_is_efficacy():
    for name, spec in SOURCES.items():
        if spec.kind == "toxicity":
            assert spec.direction in {HIGHER_IS_TOXIC, LOWER_IS_TOXIC}, name
            assert spec.threshold is not None, name
        else:
            assert spec.direction is None and spec.threshold is None, name
        assert spec.citation, name


def test_hagedorn_is_lower_is_toxic():
    """Reduced calcium oscillation is the toxic direction, not raised."""
    assert SOURCES["Hagedorn_2022"].direction == LOWER_IS_TOXIC


def test_potency_datasets_are_not_toxicity():
    for name in ("McQuisten_2007", "Ichihara_2007_1", "Ichihara_2007_2"):
        assert SOURCES[name].kind == "efficacy", name


def test_normalized_outputs_agree_with_the_manifest(root):
    manifest_path = root / "data/processed/external_normalization_manifest.json"
    if not manifest_path.exists():
        pytest.skip("normalizer has not been run in this checkout")
    manifest = json.loads(manifest_path.read_text())
    tox = pd.read_csv(root / "data/processed/external_tox_normalized.csv")
    eff = pd.read_csv(root / "data/processed/external_efficacy_normalized.csv")

    assert len(tox) == manifest["toxicity_sequences"]
    assert len(eff) == manifest["efficacy_sequences"]
    # No efficacy readout may carry a toxicity label.
    assert "toxic" not in eff.columns
    # Every toxicity row records how its label was made.
    assert tox["direction"].notna().all() and tox["threshold"].notna().all()
    # Sequences are unique within a source, so no replicate spans a split.
    assert not tox.duplicated(["source", "sequence"]).any()
