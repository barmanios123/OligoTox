import os
import sys
import json
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from oligotox.pipelines.build_training_dataset import run
from oligotox.config.paths import PROCESSED_DATA_DIR


def test_training_dataset_includes_omics_features():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    training_path = os.path.join(PROCESSED_DATA_DIR, "training_dataset.csv")
    feature_matrix = os.path.join(PROCESSED_DATA_DIR, "feature_matrix.csv")
    label_matrix = os.path.join(PROCESSED_DATA_DIR, "label_matrix.csv")
    metadata_path = os.path.join(PROCESSED_DATA_DIR, "dataset_metadata.json")

    info = run(
        training_dataset_path=training_path,
        feature_matrix_path=feature_matrix,
        label_matrix_path=label_matrix,
        metadata_path=metadata_path,
    )

    X = pd.read_csv(feature_matrix)
    with open(metadata_path, "r") as f:
        meta = json.load(f)

    omics_cols = [c for c in X.columns if c.startswith("omics_")]
    assert len(omics_cols) == meta.get("omics_features_added", len(omics_cols))
    print("Omics feature columns:", omics_cols[:10])
    print("Total omics features:", len(omics_cols))
    print("run() info:", info)


if __name__ == "__main__":
    test_training_dataset_includes_omics_features()
