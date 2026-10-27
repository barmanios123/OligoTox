import json

import pandas as pd

from oligotox.training.train import (
    ASSAY_READOUT_COLUMNS,
    LEAKY_PREFIXES,
    _drop_leaky_features,
    _get_groups_for_cv,
)


def test_drop_leaky_features():
    X = pd.DataFrame({
        "gc_content": [0.4, 0.6],
        "label_kupffer_activation": [0, 1],
        "endpoint_cns_toxicity": [1, 1],
        "cas_avg": [0.2, 0.9],
        "always_zero": [0.0, 0.0],
    })
    kept, dropped = _drop_leaky_features(X)
    assert list(kept.columns) == ["gc_content"]
    assert dropped["outcome_annotations"] == ["label_kupffer_activation", "endpoint_cns_toxicity"]
    assert dropped["assay_readouts"] == ["cas_avg"]
    assert dropped["constant"] == ["always_zero"]


def test_groups_treat_rna_and_dna_as_same_sequence():
    X = pd.DataFrame({"sequence": ["augc", "ATGC", "GGCC"]})
    groups = _get_groups_for_cv(X)
    assert groups[0] == groups[1] != groups[2]


def test_shipped_models_use_no_leaky_columns(root):
    cfg = json.loads((root / "artifacts/config/composite_config.json").read_text())
    cols = cfg["trained_feature_columns"]
    assert not [c for c in cols if c.startswith(LEAKY_PREFIXES)]
    assert not set(cols) & set(ASSAY_READOUT_COLUMNS)


def test_published_metrics_are_sequence_grouped(root):
    summary = json.loads((root / "reports/results/training_summary_after.json").read_text())
    for endpoint, s in summary.items():
        if s["pr_auc"] is not None:
            assert s["cv"] == "sequence_grouped", endpoint
            assert s["positive_sequences"] >= 3, endpoint
