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

from extract_patent_hepatotox import _rejoin, extract, sequence_map, _tables  # noqa: E402

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
    by_dose = {r["dose_mg_kg"]: r for r in rows if r["isis_number"] == 233698}
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


def test_a_compound_whose_sequence_is_absent_is_rejected_not_guessed(parsed):
    """ISIS 233693 maps to SEQ ID 131, which the rendered patent never prints.

    The sequence listing is a separate document Google Patents does not show.
    Emitting that row with a guessed or blank sequence would put three
    fabricated measurements into the training corpus.
    """
    rows, rejected, _ = parsed
    assert not [r for r in rows if r["isis_number"] == 233693]
    assert rejected["no_sequence"] == 3


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
