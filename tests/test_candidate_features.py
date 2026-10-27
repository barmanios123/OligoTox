import json

from oligotox.features.candidate import candidate_features


def test_covers_every_trained_column(root):
    cfg = json.loads((root / "artifacts/config/composite_config.json").read_text())
    feats = candidate_features("TGGCAAGCATCCTGTA", "LNA_gapmer", "LNA_PS", 1.0, 0.375)
    missing = [c for c in cfg["trained_feature_columns"] if c not in feats]
    assert missing == []


def test_matches_training_definitions():
    feats = candidate_features("TGGCAAGCATCCTGTA", "LNA_gapmer", "LNA_PS", 1.0, 0.375)
    assert feats["seq_length"] == 16
    assert feats["gc_content"] == 0.5
    assert round(feats["seed_gc"], 6) == 0.714286
    assert feats["ps_x_len"] == 16
    assert feats["lna_x_gc"] == 0.1875
    assert feats["aso_gapmer_like"] == 1.0


def test_rna_and_dna_alphabets_agree():
    rna = candidate_features("AUGCUAGCUAGCUAGCUAGC")
    dna = candidate_features("ATGCTAGCTAGCTAGCTAGC")
    assert rna == dna
