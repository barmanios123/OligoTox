"""The CNS endpoint is the one model validated on an independent assay.

These tests guard the things that would silently break it: the label
direction, the external-transfer claim, and the per-endpoint feature list.
"""
import json

import pytest

from oligotox.features.candidate import candidate_features, g_motif_features
from oligotox.models.loader import load_endpoint_evidence, load_endpoint_features
from oligotox.pipelines.score_candidate import score_candidate


@pytest.fixture
def external(root):
    return json.loads((root / "reports/results/cns_external.json").read_text())


def test_assays_share_no_sequences(external):
    """Transfer is only external if the two assays are disjoint."""
    assert external["shared_sequences"] == 0


def test_cross_assay_transfer_beats_chance_both_ways(external):
    transfers = [r for r in external["results"]
                 if r["evaluation"] == "cross-assay transfer" and r["features"].startswith("All")]
    assert len(transfers) == 2, "both directions must be reported"
    for r in transfers:
        assert r["ci_low"] > 0.5, f"{r['assay']} does not clear chance"


def test_hagedorn_label_direction_is_low_equals_toxic(external):
    """Hagedorn et al. tie toxicity to REDUCED calcium oscillation.

    The shipped oligogym_tox_combined_fixed.csv has this inverted. If someone
    re-imports that column, the label flips and this test fails.
    """
    label = external["datasets"]["Hagedorn_2022"]["label"]
    assert "<" in label, f"expected a 'less than' threshold, got {label!r}"
    hag = external["datasets"]["Hagedorn_2022"]
    rate = hag["positives"] / hag["sequences"]
    assert 0.2 < rate < 0.6, f"base rate {rate:.2f} suggests the direction flipped"


def test_cns_is_validated_and_names_its_caveat(root):
    evidence = load_endpoint_evidence(root / "artifacts/config")["endpoints"]["cns_toxicity"]
    assert evidence["tier"] == "validated"
    assert evidence["validation"] == "cross-assay external"
    # It does not beat the published G rule, and the reason has to say so.
    if not evidence["checks"]["beats_published_rule"]:
        assert "G-count rule" in evidence["reason"]


def test_cns_uses_its_own_sequence_only_features(root):
    cols = load_endpoint_features(root / "artifacts/config")["cns_toxicity"]["feature_columns"]
    assert cols, "cns_toxicity must declare its own feature list"
    # No chemistry or modality feature: the external data has no such annotation.
    for col in cols:
        assert not col.startswith("aso_"), col
        assert col not in {"ps_fraction", "lna_fraction", "ps_x_len", "lna_x_gc"}, col
    # Everything it needs is computable at scoring time.
    built = candidate_features("AUGCUAGCUAGCUAGCUAGC", "ASO", "PS", ps_fraction=1.0)
    assert [c for c in cols if c not in built] == []


def test_g_rich_sequence_scores_higher_than_g_poor(root, tmp_path):
    """Direction check against the published driver, through the real pipeline."""
    def cns(seq):
        out = score_candidate(
            candidate_features=candidate_features(seq, "ASO", "PS", ps_fraction=1.0),
            processed_dir=tmp_path / "unused",
            model_dir=root / "artifacts/models",
            config_dir=root / "artifacts/config",
        )
        return out["per_endpoint"]["cns_toxicity"]["prob"]

    assert cns("CCGGCCCCTAGCGCGCGACT") > cns("AAATCTATAATAACCACCAC")


def test_g_motif_features_are_computed_from_the_three_prime_end():
    # Nearest G to the 3' end: last base is G, so distance 0.
    assert g_motif_features("AAAAAAAAAG")["dist_3p_to_g"] == 0.0
    # No G at all: distance is the full length, not an error.
    f = g_motif_features("AAAACCCCTT")
    assert f["dist_3p_to_g"] == 10.0 and f["g_count"] == 0.0
    assert g_motif_features("AGGGA")["longest_g_run"] == 3.0
