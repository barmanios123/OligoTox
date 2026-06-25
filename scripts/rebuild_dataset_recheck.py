from pathlib import Path
from oligotox.pipelines import build_training_dataset

def main() -> None:
    build_training_dataset.run(
        training_csv_path=Path("data/processed/training_dataset_with_chem.csv"),
        feature_matrix_path=Path("data/processed/feature_matrix_with_chem_recheck.csv"),
        label_matrix_path=Path("data/processed/label_matrix_with_chem_recheck.csv"),
        dataset_metadata_path=Path("data/processed/dataset_metadata_with_chem_recheck.json"),
    )
    print("dataset recheck complete")

if __name__ == "__main__":
    main()
