import os
import sys
import json
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from oligotox.pipelines.build_training_dataset import run
from oligotox.config.paths import PROCESSED_DATA_DIR


def test_omics_columns_flow_into_feature_matrix():
    run(
        str(PROCESSED_DATA_DIR / "training_dataset.csv"),
        str(PROCESSED_DATA_DIR / "feature_matrix.csv"),
        str(PROCESSED_DATA_DIR / "label_matrix.csv"),
        str(PROCESSED_DATA_DIR / "dataset_metadata.json"),
    )

    X = pd.read_csv(PROCESSED_DATA_DIR / "feature_matrix.csv")
    meta = json.loads((PROCESSED_DATA_DIR / "dataset_metadata.json").read_text())

    assert "omics_modalities_loaded" in meta
    assert "omics_feature_columns" in meta

    for col in meta["omics_feature_columns"]:
        assert col in X.columns
