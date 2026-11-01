"""The app must be able to say what chemistry a candidate actually is.

The Modality dropdown offered only siRNA and ASO, and `lna_fraction` was
never set, so every candidate reached the models with lna_fraction 0.0 and
none of the gapmer flags. The hepatocellular model consumes both, and its
training rows range from 0.125 to 0.625 LNA, so an LNA gapmer was scored as
a chemistry that does not occur in the corpus.

The external screen was scored the same way. Its 149 compounds are uniform
3-8-3 LNA gapmers, and scoring them as plain ASOs is most of why the
endpoint looked like chance on them: ROC-AUC 0.457 told that way, 0.539
told correctly.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import oligotox_app as app  # noqa: E402


def test_the_screen_s_own_design_round_trips():
    """3-8-3 on a 14-mer is 6 of 14, which is what Table S1 encodes."""
    assert app._lna_fraction("LNA_gapmer", "3-8-3", 14) == 6 / 14


def test_a_moe_gapmer_carries_no_lna():
    """MOE wings are not LNA. The patent rows are ASO_gapmer at zero."""
    assert app._lna_fraction("ASO_gapmer", "5-10-5", 20) == 0.0


def test_a_non_gapmer_carries_no_lna():
    assert app._lna_fraction("siRNA", "3-8-3", 20) == 0.0


def test_an_unparseable_motif_does_not_invent_a_fraction():
    for motif in ("", "abc", "3-8", "lots"):
        assert app._lna_fraction("LNA_gapmer", motif, 14) == 0.0


def test_the_fraction_reaches_the_features_the_model_sees():
    """The whole point: it must not stop at the sidebar."""
    seq = "GTATTAGCCCCTAG"
    plain = app._sequence_features(seq, "PS", 3.0, "ASO")
    gapmer = app._sequence_features(seq, "PS", 3.0, "LNA_gapmer", "3-8-3")
    assert plain["lna_fraction"] == 0.0
    assert gapmer["lna_fraction"] == 6 / 14
    assert gapmer["lna_fraction"] != plain["lna_fraction"]


def test_the_gapmer_modalities_are_offered():
    source = (ROOT / "oligotox_app.py").read_text()
    assert '"ASO_gapmer"' in source and '"LNA_gapmer"' in source, (
        "the app cannot express the chemistry of the compounds it is for"
    )


def test_the_external_check_is_reported_whichever_side_of_chance_it_lands():
    """It used to appear only when the result was below 0.5.

    The moment the number rose above chance, the only external check this
    endpoint has dropped out of its evidence entry. Silence reads as "it
    transfers", which 5 positives cannot support in either direction.
    """
    import json

    config = json.loads(
        (ROOT / "artifacts/config/endpoint_evidence.json").read_text())
    spec = config["endpoints"].get("hepatocellular_injury", {})
    external = spec.get("external_check")
    if not external:
        import pytest
        pytest.skip("no external check in this build")
    assert "independent screen" in spec["reason"]
    assert f"{external['roc_auc']:.2f}" in spec["reason"]
    assert "0.50 is chance" in spec["reason"]
    assert str(external["positives"]) in spec["reason"]
