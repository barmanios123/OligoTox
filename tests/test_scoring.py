import math

from oligotox.features.candidate import candidate_features
from oligotox.pipelines.score_candidate import score_candidate


def _score(root, tmp_path, seq):
    return score_candidate(
        candidate_features=candidate_features(seq, "siRNA", "PS", ps_fraction=0.9),
        processed_dir=tmp_path / "no_training_data_needed",
        model_dir=root / "artifacts/models",
        config_dir=root / "artifacts/config",
    )


def test_scores_from_fresh_clone(root, tmp_path):
    out = _score(root, tmp_path, "AUGCUAGCUAGCUAGCUAGC")
    assert len(out["per_endpoint"]) == 13
    assert 0.0 <= out["composite_score"] <= 1.0
    for info in out["per_endpoint"].values():
        assert 0.0 <= info["prob"] <= 1.0 and not math.isnan(info["prob"])


def test_scoring_is_deterministic(root, tmp_path):
    a = _score(root, tmp_path, "GACUUAUAACACUGAUGAAA")
    b = _score(root, tmp_path, "GACUUAUAACACUGAUGAAA")
    assert a["composite_score"] == b["composite_score"]


def test_different_sequences_get_different_scores(root, tmp_path):
    a = _score(root, tmp_path, "GGGGCCCCGGGGCCCCGGGG")
    b = _score(root, tmp_path, "AUAUAUAUAUAUAUAUAUAU")
    assert a["composite_score"] != b["composite_score"]
