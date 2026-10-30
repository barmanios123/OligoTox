"""The patent extractor is pinned against a patent read by hand.

Every value below was verified against the rendered table in US11225664B2
before it was written here:

    Saline            ALT  22   AST  47
    ISIS 233698 10    ALT  28   AST  55
    ISIS 233698 25    ALT  16   AST  53
    ISIS 233698 50    ALT 192   AST 175

An extractor is only worth running if it is right, and "it produced rows" is
not evidence of that. These assert the numbers, not the row count.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from extract_patent_hepatotox import (  # noqa: E402
    _column_names, _find_column, _groups, _rejoin, _tables, compound_details,
    extract, prose_dose, sequence_map,
)

FIXTURE = ROOT / "data" / "external" / "patent_cache" / "US11225664B2.html"
pytestmark = pytest.mark.skipif(not FIXTURE.exists(),
                                reason="patent fixture not cached")


@pytest.fixture(scope="module")
def parsed():
    html = FIXTURE.read_text(encoding="utf-8", errors="replace")
    rows, rejected = extract("US11225664B2", html)
    return rows, rejected, html


def test_the_numbers_match_the_patent(parsed):
    rows, _, _ = parsed
    by_dose = {r["dose_mg_kg"]: r for r in rows
               if r["isis_number"] == 233698 and r["table"] == 28}
    assert set(by_dose) == {10.0, 25.0, 50.0}
    assert (by_dose[10.0]["alt"], by_dose[10.0]["ast"]) == (28.0, 55.0)
    assert (by_dose[25.0]["alt"], by_dose[25.0]["ast"]) == (16.0, 53.0)
    assert (by_dose[50.0]["alt"], by_dose[50.0]["ast"]) == (192.0, 175.0)
    assert all(r["control_alt"] == 22.0 for r in by_dose.values())


def test_the_sequence_is_the_one_the_patent_names(parsed):
    """ISIS 233698 is SEQ ID NO 136, and 136 is this 20-mer."""
    rows, _, _ = parsed
    seq = {r["sequence"] for r in rows if r["isis_number"] == 233698}
    assert seq == {"CTGACCTTGTGTTGAAGGGC"}
    assert {r["seq_id"] for r in rows if r["isis_number"] == 233698} == {136}


def test_the_compound_table_is_preferred_over_the_prose_route():
    """ISIS 233693 is recovered from the table that lists it directly.

    The prose says "ISIS 233693 (SEQ ID NO: 131)", but reading SEQ ID 131
    back out of the sequence listing failed, because that table ends in two
    numeric columns and taking the last one picked 355. The row itself,
    `233693 | 6 | 349 | GACATGTTCTTCACCTCCTC | 80 | 131 | 355`, carries the
    sequence next to the compound and needs no lookup at all. Values checked
    against the patent.
    """
    html = FIXTURE.read_text(encoding="utf-8", errors="replace")
    rows, _ = extract("US11225664B2", html)
    by_dose = {r["dose_mg_kg"]: r for r in rows
               if r["isis_number"] == 233693 and r["table"] == 28}
    assert set(by_dose) == {10.0, 25.0, 50.0}
    assert (by_dose[10.0]["alt"], by_dose[25.0]["alt"], by_dose[50.0]["alt"]) \
        == (23.0, 21.0, 52.0)
    assert {r["sequence"] for r in by_dose.values()} == {"GACATGTTCTTCACCTCCTC"}


def test_a_second_study_in_one_patent_is_kept_apart_from_the_first(parsed):
    """ISIS 233693 is dosed twice in US11225664B2, against different controls.

    Table 28 gives 10/25/50 mg/kg against a saline ALT of 22; table 64 gives
    12.5/25/50 against a PBS ALT of 57. Both were read off the patent. Keyed
    by dose alone the two studies overwrite each other, which is what made
    one sequence look as though it carried contradictory ALT at one dose.
    """
    rows, _, _ = parsed
    second = {r["dose_mg_kg"]: r for r in rows
              if r["isis_number"] == 233693 and r["table"] == 64}
    assert set(second) == {12.5, 25.0, 50.0}
    assert (second[50.0]["alt"], second[50.0]["ast"]) == (138.0, 152.0)
    assert (second[25.0]["alt"], second[25.0]["ast"]) == (43.0, 70.0)
    assert (second[12.5]["alt"], second[12.5]["ast"]) == (45.0, 87.0)
    assert all(r["control_alt"] == 57.0 for r in second.values())


def test_every_row_carries_the_table_it_came_from(parsed):
    """Without it two studies of one compound cannot be told apart."""
    rows, _, _ = parsed
    assert all(r.get("table") for r in rows)


def test_two_compounds_may_share_a_sequence_and_differ_by_motif():
    """This looked like a mapping bug and is not one.

    Ionis files several compounds on one base sequence that differ only in
    gapmer architecture, so `426183 | 3-14-3` and `426261 | 2-13-5` both carry
    AGCCATGGTGATCAGGAGGC. Dropping the motif column is what made them
    indistinguishable, and a corpus that cannot tell them apart would hold
    contradictory labels for one sequence.
    """
    cached = FIXTURE.parent / "US9567587B2.html"
    if not cached.exists():
        pytest.skip("US9567587B2 not cached")
    details = compound_details(_tables(
        cached.read_text(encoding="utf-8", errors="replace")))
    a, b = details.get(426183), details.get(426261)
    assert a and b
    assert a["sequence"] == b["sequence"] == "AGCCATGGTGATCAGGAGGC"
    assert a["motif"] != b["motif"], "the motif is what tells them apart"
    assert {a["motif"], b["motif"]} == {"3-14-3", "2-13-5"}


def test_dose_rows_are_rejoined_with_their_label():
    """"ISIS 233693," and "50 mg/kg" arrive as separate rows."""
    merged = _rejoin([
        ["Saline", "22", "47"],
        ["ISIS 233693,", "52", "79"],
        ["50 mg/kg"],
    ])
    assert merged[2 - 1][0] == "ISIS 233693, 50 mg/kg"
    assert len(merged) == 2


def test_a_measurement_without_a_control_is_not_emitted():
    """An ALT of 52 means nothing without the same study's baseline."""
    html = "<table><tr><td>Treatment</td><td>ALT</td></tr></table>" \
           "<table><tr><td>ISIS 233698, 50 mg/kg</td><td>192</td></tr></table>"
    rows, rejected = extract("TEST", html)
    assert rows == []
    assert rejected["no_control"] >= 1 or rejected["no_seq_id"] >= 1


def test_sequence_map_reads_the_seq_id_column(parsed):
    _, _, html = parsed
    mapping = sequence_map(_tables(html))
    assert mapping[136] == "CTGACCTTGTGTTGAAGGGC"
    assert all(set(v) <= set("ACGT") for v in mapping.values())


def test_the_dose_column_is_read_when_the_table_has_one():
    """US9567587B2 table 68 is a dose-response table with its own column.

    Checked against the patent: PBS — 52, ISIS 426325 at 8 -> 51, at 20 -> 38.
    Reading the label alone saw one unexplained ALT per compound and lost the
    second dose entirely, because the lower dose leaves the label blank.
    """
    cached = FIXTURE.parent / "US9567587B2.html"
    if not cached.exists():
        pytest.skip("US9567587B2 not cached")
    rows, _ = extract("US9567587B2",
                      cached.read_text(encoding="utf-8", errors="replace"))
    got = {r["dose_mg_kg"]: r for r in rows if r["isis_number"] == 426325}
    assert set(got) == {8.0, 20.0}
    assert (got[8.0]["alt"], got[8.0]["ast"]) == (51.0, 62.0)
    assert (got[20.0]["alt"], got[20.0]["ast"]) == (38.0, 46.0)
    assert all(r["dose_source"] == "table" for r in got.values())
    assert all(r["control_alt"] == 52.0 for r in got.values())


def test_a_label_split_across_rows_still_names_its_compound():
    """"ISIS" and "304299" print on separate rows of a dose-response table.

    US9399774B2 table 118 renders the compound over three rows, one per dose.
    Reading row[0] matched neither fragment, so the whole table was rejected
    as nameless; the group has to be assembled before the label means
    anything.
    """
    rows = _groups([
        ["PBS", "—", "48", "112"],
        ["ISIS", "25", "42", "93"],
        ["304299", "10", "37", "56"],
        ["", "4", "35", "71"],
    ], 0)
    assert [label for label, _ in rows] == ["PBS", "ISIS 304299"]
    assert len(rows[1][1]) == 3


def test_a_dose_column_is_matched_to_cells_from_the_right():
    """Data rows often carry a leading spacer cell the header does not have.

    Aligning the header to the data from the left read the dose column as the
    label and the ALT column as the dose.
    """
    names = _column_names([["", "Dose", "ALT", "AST"], ["", "(mg/kg)", "(IU/L)", "(IU/L)"]])
    assert names == ["", "Dose (mg/kg)", "ALT (IU/L)", "AST (IU/L)"]
    assert _find_column(names, r"\bdose\b") == 1
    assert _find_column(names, r"\bALT\b") == 2


def test_an_ambiguous_prose_dose_is_refused():
    """Two doses in the passage means the row's dose is unknown, not the first."""
    one = "Example 9 mice were injected with 50 mg/kg of ISIS 304538 weekly."
    two = "Example 9 mice were injected with 25 mg/kg or 50 mg/kg of ISIS 304538."
    assert prose_dose(one, len(one), {304538}) == 50.0
    assert prose_dose(two, len(two), {304538}) is None


def test_prose_is_not_borrowed_from_an_unrelated_passage():
    """The window must mention the compounds the table is about."""
    text = "Example 9 mice were injected with 50 mg/kg of ISIS 999999 weekly."
    assert prose_dose(text, len(text), {304538}) is None
