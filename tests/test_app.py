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


def test_runs_are_not_written_to_shared_storage(root):
    """A public deployment must not show one visitor's sequences to the next.

    Runs used to be written to artifacts/runs/ and the Saved runs tab listed
    every file it found there.
    """
    source = (root / "oligotox_app.py").read_text()
    assert "RUNS_DIR" not in source, "the app still references a shared runs directory"
    assert 'st.session_state.setdefault("runs"' in source

    at = _run(root)
    runs = at.session_state["runs"]
    assert len(runs) == 1, "the run should live in session state"
    # Nothing on disk for this run.
    run_id = runs[0]["run_id"]
    assert not list((root / "artifacts" / "runs").glob(f"{run_id}.json"))


def test_the_summary_reports_the_backbone_that_was_actually_scored():
    """The interface must not print a parameter the model never saw.

    `_sequence_features` scores PS chemistry at 0.9, while the summary and the
    design scan both defaulted to 1.0 when the run payload carried no value.
    The app then displayed "100% phosphorothioate" for a candidate scored at
    90%. One definition, used everywhere.
    """
    # oligotox_app.py sits at the repository root, which is not importable
    # on its own: the other tests here reach the app through Streamlit's
    # AppTest rather than importing it, so nothing else puts the root on the
    # path. Locally it happened to be there and on CI it was not.
    import sys
    from pathlib import Path

    root = str(Path(__file__).resolve().parent.parent)
    if root not in sys.path:
        sys.path.insert(0, root)
    import oligotox_app as app

    assert app._ps_fraction("PS") == 0.9
    assert app._ps_fraction("mixed_PS_PO") == 0.9
    assert app._ps_fraction("PO") == 0.0

    feats = app._sequence_features("ACGTACGTACGTACGTACGT", "PS", 3.0, "siRNA")
    assert feats["ps_fraction"] == app._ps_fraction("PS"), (
        "the features the model sees must carry the same PS fraction the "
        "summary prints"
    )


def test_how_it_works_is_read_from_the_shipped_artifacts(root):
    """The trust page must not carry hand-written numbers.

    A page that explains the evidence is worth less than nothing if it can
    drift from the evidence, so every figure on it is derived at render time.
    These assert the derivation rather than the text: change the config and
    the page has to change with it.
    """
    at = _run(root)
    evidence = load_endpoint_evidence(root / "artifacts/config")
    validated = [e for e, v in evidence["endpoints"].items() if v["tier"] == "validated"]
    body = " ".join(m.value for m in at.markdown)

    assert (f"**{len(validated)} of {len(evidence['endpoints'])}** endpoints "
            "are presented") in body
    for endpoint in validated:
        assert f"`{endpoint}`" in body

    inventory = json.loads(
        (root / "reports/results/dataset_inventory.json").read_text())
    usable = [s for s in inventory["sources"] if s["usable_for_training"]]
    assert (f"**{len(usable)} source(s), "
            f"{inventory['usable_sequences']:,} sequences.**") in body


def test_the_page_says_grouped_cv_beats_the_flattering_split(root):
    """The leakage comparison is the page's main claim; it has to be counted.

    Writing "the random split scores higher" as prose would survive the
    numbers changing. This recounts it from the training summary.
    """
    at = _run(root)
    training = json.loads(
        (root / "reports/results/training_summary_after.json").read_text())
    scored = [s for s in training.values()
              if isinstance(s, dict) and s.get("pr_auc") is not None]
    inflated = [s for s in scored
                if s.get("random_split_pr_auc") is not None
                and s["random_split_pr_auc"] > s["pr_auc"]]
    body = " ".join(m.value for m in at.markdown)
    assert f"**{len(inflated)} of {len(scored)}**" in body


def test_the_external_check_is_on_the_page_and_not_softened(root):
    """The one result that argues against the app has to be visible on it."""
    at = _run(root)
    evidence = load_endpoint_evidence(root / "artifacts/config")
    external = evidence["endpoints"]["hepatocellular_injury"]["external_check"]
    warnings = " ".join(w.value for w in at.warning)
    assert f"{external['roc_auc']:.2f}" in warnings
    assert "0.5 is chance" in warnings


def test_adding_a_tab_cannot_renumber_the_others(root):
    """The panels were indexed by position and are now keyed by name."""
    source = (root / "oligotox_app.py").read_text()
    assert "with tabs[" not in source, "a positional tab index came back"
    assert "How it works" in source
