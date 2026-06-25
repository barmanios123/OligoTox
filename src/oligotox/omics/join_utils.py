from __future__ import annotations

from pathlib import Path
import pandas as pd


def safe_left_join_on_candidate(
    base_df: pd.DataFrame,
    omics_path: str | Path,
    id_column: str = "candidate_id",
    prefix: str = "omics_",
) -> pd.DataFrame:
    """
    Left-join candidate-level omics features onto base_df.
    - base_df: rows are candidates (must contain id_column).
    - omics_path: CSV with id_column + omics feature columns.
    - Returns base_df with extra columns prefixed.
    """
    omics_path = Path(omics_path)
    if not omics_path.exists():
        return base_df

    omics = pd.read_csv(omics_path)
    if omics.empty:
        return base_df

    if id_column not in base_df.columns:
        raise ValueError(f"Base DF missing id column '{id_column}' for omics join.")
    if id_column not in omics.columns:
        raise ValueError(f"Omics DF missing id column '{id_column}' for join.")

    # Drop duplicate candidate rows in omics if any
    omics = omics.drop_duplicates(subset=[id_column])

    # Avoid collisions with base_df columns
    rename = {}
    for c in omics.columns:
        if c == id_column:
            continue
        if c in base_df.columns or not c.startswith(prefix):
            rename[c] = prefix + c
    if rename:
        omics = omics.rename(columns=rename)

    merged = base_df.merge(omics, on=id_column, how="left")
    return merged
