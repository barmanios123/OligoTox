"""Design advice must only come from endpoints the evidence supports.

An edit recommended to move a number that means nothing is worse than no
advice, so these pin the gate and the measurement.
"""
import json

import pytest

from oligotox.advice import explain, suggest

SEQ = "CCGGCCCCTAGCGCGCGACT"  # G-rich, high predicted CNS risk


@pytest.fixture
def dirs(root):
    return {"model_dir": root / "artifacts/models", "config_dir": root / "artifacts/config"}


def test_only_validated_endpoints_are_advised_on(root, dirs):
    evidence = json.loads((root / "artifacts/config/endpoint_evidence.json").read_text())
    validated = {e for e, v in evidence["endpoints"].items() if v["tier"] == "validated"}

    out = explain(SEQ, "ASO", **dirs)
    assert set(out["validated_toxicity"]) == validated

    advice = suggest(SEQ, "ASO", top_n=5, **dirs)
    assert set(advice.validated_endpoints) == validated
    for edit in advice.edits:
        assert set(edit.toxicity_change) <= validated


def test_a_suggested_edit_really_produces_the_change_it_claims(dirs):
    """The claimed delta must equal rescoring the edited sequence."""
    advice = suggest(SEQ, "ASO", top_n=3, **dirs)
    assert advice.edits, "expected at least one edit on a G-rich sequence"

    edit = advice.edits[0]
    # The stated sequence is the stated edit applied to the original.
    assert edit.sequence[edit.position - 1] == edit.to_base
    assert advice.sequence[edit.position - 1] == edit.from_base
    assert edit.sequence != advice.sequence

    rescored = explain(edit.sequence, "ASO", **dirs)["validated_toxicity"]
    for endpoint, claimed in edit.toxicity_change.items():
        actual = rescored[endpoint] - advice.baseline_toxicity[endpoint]
        assert abs(actual - claimed) < 1e-4, endpoint


def test_trade_offs_are_reported_not_hidden(dirs):
    """An edit that costs potency must say so."""
    advice = suggest(SEQ, "ASO", top_n=12, **dirs)
    assert advice.baseline_potency, "per-assay potency should apply to an ASO"
    # At least one edit should move potency; otherwise the trade-off is untested.
    assert any(e.potency_change for e in advice.edits)
    for edit in advice.edits:
        if edit.potency_change:
            assert edit.describe().endswith("%"), "potency effect must be in the summary"


def test_edits_are_ranked_by_risk_reduction(dirs):
    advice = suggest(SEQ, "ASO", top_n=8, **dirs)
    reductions = [e.total_toxicity_reduction for e in advice.edits]
    assert reductions == sorted(reductions, reverse=True)


def test_g_rich_cns_advice_targets_g(dirs):
    """Sanity against the published rule the CNS model reproduces."""
    advice = suggest(SEQ, "ASO", top_n=5, **dirs)
    improving = [e for e in advice.edits if e.total_toxicity_reduction > 0]
    assert improving, "expected some risk-reducing edits"
    assert sum(e.from_base == "G" for e in improving) >= len(improving) // 2


def test_a_bad_sequence_is_refused(dirs):
    with pytest.raises(ValueError, match="A/C/G/T"):
        suggest("NOTASEQUENCE123", "ASO", **dirs)


def test_advice_always_states_its_limits(dirs):
    advice = suggest(SEQ, "ASO", top_n=3, **dirs)
    joined = " ".join(advice.notes).lower()
    assert "not a coefficient" in joined
    assert advice.to_dict()["validated_endpoints"]
