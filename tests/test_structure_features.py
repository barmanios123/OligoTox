"""Structure and off-target features, and the drift they could cause."""
import json
from pathlib import Path

import pandas as pd
import pytest

from oligotox.features.accessibility import find_site, site_accessibility
from oligotox.features.offtarget import burden, offtarget_features, reverse_complement
from oligotox.features.structure_features import build_structure_features

SEQ = "CCGGCCCCTAGCGCGCGACT"


def test_real_folding_energies_are_used_not_the_placeholder():
    """A GC-derived placeholder where ViennaRNA was expected is silent drift."""
    out = build_structure_features(pd.DataFrame({"sequence": [SEQ]}), cache_path=None)
    assert out["struct_viennarna_used"].iloc[0] == 1.0
    # A real MFE for this sequence is negative; the placeholder is not a fold.
    assert out["struct_mfe_proxy"].iloc[0] < 0


def test_training_matrix_was_built_with_real_folding(root):
    """Training and scoring must agree on how struct_mfe_proxy was made."""
    matrix = root / "data/processed/feature_matrix.csv"
    if not matrix.exists():
        pytest.skip("no feature matrix in this checkout")
    frame = pd.read_csv(matrix, usecols=lambda c: c == "struct_mfe_proxy")
    if "struct_mfe_proxy" not in frame.columns:
        pytest.skip("struct_mfe_proxy not among the trained features")
    # The placeholder is non-negative by construction; real MFEs are not.
    assert (frame["struct_mfe_proxy"] < 0).any(), (
        "no negative folding energies: the matrix looks built with the placeholder"
    )


def test_reverse_complement_round_trips():
    assert reverse_complement(reverse_complement(SEQ)) == SEQ
    assert reverse_complement("AAAA") == "TTTT"
    # U is read as T, so RNA and DNA spellings agree.
    assert reverse_complement("AUGC") == reverse_complement("ATGC")


def test_offtarget_counts_real_matches_not_a_gc_proxy():
    """The old implementation returned GC content rescaled."""
    background = "N".join(["AAAA" + reverse_complement(SEQ) + "TTTT", "GGGGCCCCGGGG"])
    result = burden(SEQ, background, n_transcripts=2, megabases=len(background) / 1e6)
    assert result.full_hits == 1, "the planted full-length site must be found"

    # A sequence absent from the background must score zero, whatever its GC.
    absent = burden("TTTTTTTTTTTTTTTTTTTT", background, 2, len(background) / 1e6)
    assert absent.full_hits == 0


def test_offtarget_without_a_cache_reports_unavailable_not_zero_risk(tmp_path):
    feats = offtarget_features(SEQ, transcript_dir=tmp_path / "missing")
    assert feats["offtarget_available"] == 0.0
    assert feats["offtarget_full_hits"] == 0.0


def test_binding_site_is_found_by_reverse_complement():
    transcript = "GGGGGGGG" + reverse_complement(SEQ) + "CCCCCCCC"
    assert find_site(SEQ, transcript) == 8
    # The forward sequence is not the binding site.
    assert find_site(SEQ, "GGGG" + SEQ + "CCCC") == -1


def test_accessibility_reports_absence_rather_than_inaccessibility():
    """An unlocatable site must not look like a paired one."""
    access = site_accessibility(SEQ, "ACGT" * 50)
    assert access.site_found is False
    assert access.as_features()["site_found"] == 0.0


def test_an_open_site_scores_higher_than_a_buried_one():
    """Accessibility must track structure, not just return a constant."""
    site = reverse_complement(SEQ)
    # Flanked by unstructured A-runs the site is relatively open.
    open_context = "A" * 120 + site + "A" * 120
    # Flanked by its own complement it is sequestered in a long duplex.
    closed_context = "A" * 60 + reverse_complement(site) + site + "A" * 60
    opened = site_accessibility(SEQ, open_context)
    closed = site_accessibility(SEQ, closed_context)
    assert opened.site_found and closed.site_found
    assert opened.accessibility > closed.accessibility


def test_accessibility_result_is_recorded_with_its_coverage(root):
    path = root / "reports/results/accessibility.json"
    if not path.exists():
        pytest.skip("accessibility not evaluated in this checkout")
    payload = json.loads(path.read_text())
    # Coverage is the caveat that matters; it must never be dropped.
    assert 0.0 < payload["coverage"]["coverage"] < 1.0
    assert "intron" in payload["coverage"]["note"].lower()
