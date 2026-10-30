"""The patent label rule has to be the one the manifest claims.

A label built from a readout is a modelling decision, and the decision here
is a fold change against each study's own control. These pin the rule rather
than the row count, because "it added 444 rows" says nothing about whether
any of them is labelled correctly.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from ingest_patent_hepatotox import ALT_FOLD, LABEL, SOURCE, build  # noqa: E402


@pytest.fixture
def patents():
    """Two studies of one compound, with opposite outcomes.

    Taken from the shape of US11225664B2: the same sequence is dosed in two
    tables against different controls, and the worst dose in each is what
    labels that study.
    """
    return pd.DataFrame([
        # Study one: peaks at 192 against a control of 22, which is 8.7x.
        {"publication": "US1", "table": 28, "sequence": "ACGTACGTACGTACGTACGT",
         "isis_number": 233698, "motif": "5-10-5", "chemistry": None,
         "alt": 28.0, "control_alt": 22.0, "alt_fold_vs_control": 1.27},
        {"publication": "US1", "table": 28, "sequence": "ACGTACGTACGTACGTACGT",
         "isis_number": 233698, "motif": "5-10-5", "chemistry": None,
         "alt": 192.0, "control_alt": 22.0, "alt_fold_vs_control": 8.73},
        # Study two: never leaves the normal range.
        {"publication": "US1", "table": 64, "sequence": "ACGTACGTACGTACGTACGT",
         "isis_number": 233698, "motif": "5-10-5", "chemistry": None,
         "alt": 43.0, "control_alt": 57.0, "alt_fold_vs_control": 0.75},
    ])


def test_the_worst_dose_in_a_study_is_what_labels_it(patents):
    rows = build(patents)
    by_table = {r["notes"].split("table ")[1][:2]: r for _, r in rows.iterrows()}
    assert by_table["28"][LABEL] == 1.0, "8.7x its own control is a signal"
    assert by_table["64"][LABEL] == 0.0, "0.75x is not"


def test_two_studies_of_one_compound_stay_two_rows(patents):
    """Collapsing them would force a verdict the data does not give."""
    rows = build(patents)
    assert len(rows) == 2
    assert rows["sequence"].nunique() == 1


def test_the_label_is_a_fold_change_not_an_absolute_alt(patents):
    """Control ALT runs from 9 to 57 across these patents.

    An absolute cut would label whole studies by their baseline, so a high
    ALT against a correspondingly high control must stay negative.
    """
    high_baseline = patents.iloc[[2]].copy()
    high_baseline["alt"] = 120.0
    high_baseline["control_alt"] = 100.0
    high_baseline["alt_fold_vs_control"] = 1.2
    rows = build(high_baseline)
    assert rows.iloc[0][LABEL] == 0.0, "120 U/l at a control of 100 is not a signal"
    assert ALT_FOLD == 3.0


def test_a_cet_wing_is_not_recorded_as_a_plain_moe_gapmer(patents):
    """The patents that state their chemistry are taken at their word."""
    stated = patents.iloc[[1]].copy()
    stated["chemistry"] = "Deoxy, MOE and cEt"
    assert build(stated).iloc[0]["modality"] == "LNA_gapmer"
    assert build(patents.iloc[[1]]).iloc[0]["modality"] == "ASO_gapmer"


def test_the_unknown_chemistry_is_left_empty_rather_than_zeroed(patents):
    """Entry 13: an imputed zero reads as a measurement of zero.

    It cannot be fixed here, since the patent does not say, but it must not
    be written down as though it were known.
    """
    rows = build(patents)
    assert rows["lna_fraction"].isna().all()


def test_the_ingested_rows_are_in_the_corpus_under_their_own_source():
    """So they can be counted, audited and taken back out."""
    training = pd.read_csv(ROOT / "data/processed/training_dataset.csv",
                           low_memory=False)
    patent_rows = training[training["source"] == SOURCE]
    if patent_rows.empty:
        pytest.skip("patent rows are not currently ingested")
    assert patent_rows[LABEL].notna().all(), "every patent row carries a label"
    assert patent_rows["notes"].str.contains("vs control").all(), (
        "each row must record the measurement it was labelled from"
    )
