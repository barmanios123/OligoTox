"""The written summary must not quietly drop or distort anything.

Two bugs in the first version of `summarize` are what these pin down. Three
parameters were listed under keys the feature pipeline does not use, so they
silently never appeared, which looks identical to a sequence that has none.
And an off-target probability of 8.5e-06 was printed as "0.000", which reads
as impossible rather than as very unlikely.
"""
from __future__ import annotations

import pytest

from oligotox.advice import summarize
from oligotox.advice.advisor import PARAMETER_NOTES, _format_probability

SEQ = "GCATTGGTACGTCATTGCAA"


@pytest.fixture(scope="module")
def result():
    return summarize(SEQ)


def test_every_named_parameter_is_actually_produced(result):
    """A key that stops matching the pipeline must not vanish in silence."""
    assert result["parameters_unavailable"] == [], (
        "these parameter names no longer match candidate_features: "
        f"{result['parameters_unavailable']}"
    )
    assert len(result["parameters"]) == len(PARAMETER_NOTES)


def test_each_parameter_carries_its_meaning(result):
    """A number called dist_3p_to_g means nothing on its own."""
    for item in result["parameters"]:
        assert item["meaning"], f"{item['key']} has no explanation"
        assert item["label"] != item["key"], f"{item['key']} was never given a label"


def test_a_small_probability_is_not_rounded_to_zero():
    assert _format_probability(8.47e-06) == "8.47e-06"
    assert _format_probability(0.000470) == "4.70e-04"
    # Ordinary probabilities stay readable as decimals.
    assert _format_probability(0.6225) == "0.623"
    assert _format_probability(0.001) == "0.001"


def test_the_text_never_shows_a_bare_zero_probability(result):
    """Guards the regression directly in the rendered output."""
    for endpoint, value in result["validated_toxicity"].items():
        if value < 0.001:
            assert f"{endpoint}: 0.000" not in result["text"]


def test_only_validated_endpoints_are_presented(result):
    """The gate the whole project rests on, restated at the output layer."""
    assert set(result["validated_toxicity"]) <= {
        "cns_toxicity", "hepatocellular_injury", "offtarget_gene_silencing",
    }
    assert "deliberately not shown" in result["text"]


def test_the_summary_says_what_it_cannot_tell_you(result):
    """Accessibility and off-target burden need a transcript we were not given."""
    assert "does not tell you" in result["text"]
    assert "accessibility" in result["text"].lower()


def test_parameters_are_separated_from_predictions(result):
    """Descriptive and predictive numbers must not read as the same kind."""
    text = result["text"]
    assert text.index("Measured parameters") < text.index("Predictions the evidence")


def test_sequence_is_normalised(result):
    assert result["sequence"] == SEQ
    assert summarize("gcauuggu acgucauugcaa".replace(" ", ""))["sequence"] == SEQ
