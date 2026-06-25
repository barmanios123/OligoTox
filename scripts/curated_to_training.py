from __future__ import annotations

from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CURATED = ROOT / "data/curated/curated_oligo_tox_master_v3_plus_pmo_papargyri.csv"
OUT_TRAINING = ROOT / "data/processed/training_dataset.csv"


def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = (
        df.columns
        .map(str)
        .str.strip()
        .str.replace("\n", " ", regex=False)
        .str.split()
        .str.join("_")
        .str.lower()
    )
    return df


def main() -> None:
    if not CURATED.exists():
        raise FileNotFoundError(f"Curated master not found at {CURATED}")

    df = pd.read_csv(CURATED, low_memory=False)
    df = normalize_columns(df)

    if "candidate_id" not in df.columns:
        raise ValueError("Expected 'candidate_id' in curated master; unify script must be fixed.")

    OUT_TRAINING.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_TRAINING, index=False)

    print(f"Wrote base training dataset: {OUT_TRAINING} (rows={len(df)}, cols={df.shape[1]})")


if __name__ == "__main__":
    main()
