import json

from streamlit.testing.v1 import AppTest

from oligotox.models.loader import load_endpoint_evidence


def _run(root):
    at = AppTest.from_file(str(root / "oligotox_app.py"), default_timeout=120)
    at.run()
    assert not at.exception
    at.button[0].click().run()
    assert not at.exception
    assert not [e.value for e in at.error]
    return at


def test_app_runs_an_analysis(root):
    _run(root)


def test_predictions_table_shows_only_validated_endpoints(root):
    """The tab a user reads first must not list an unsupported endpoint."""
    at = _run(root)
    evidence = load_endpoint_evidence(root / "artifacts/config")
    validated = {e for e, v in evidence["endpoints"].items() if v["tier"] == "validated"}

    predictions = at.dataframe[0].value
    assert set(predictions["endpoint"]) == validated

    unsupported = at.dataframe[1].value
    assert not set(unsupported["endpoint"]) & validated
    assert len(unsupported) == len(evidence["endpoints"]) - len(validated)
    # Every excluded endpoint has to say why, in the table itself.
    assert unsupported["why it is not a prediction"].notna().all()


def test_composite_caption_states_the_coverage(root):
    at = _run(root)
    captions = " ".join(c.value for c in at.caption)
    assert f"Composite covers {len(at.dataframe[0].value)} of 13 endpoints" in captions


def test_saved_run_report_is_gated(root):
    at = _run(root)
    report = json.loads(at.session_state["last_run"]["report"]["json"])
    evidence = report["evidence"]
    assert evidence["endpoints_supporting_assessment"]
    assert evidence["endpoints_not_supported"]
    for note in report["regulatory_notes"]:
        assert "Immunogenicity" not in note, "raised from an unsupported endpoint"
