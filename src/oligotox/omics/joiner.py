from __future__ import annotations

import pandas as pd

from oligotox.omics.schemas import OmicsTable


def join_omics_tables(
    base_df: pd.DataFrame,
    omics_tables: list[OmicsTable],
    id_column: str = "candidate_id",
) -> tuple[pd.DataFrame, list[str], dict]:
    merged = base_df.copy()
    all_feature_columns: list[str] = []
    modality_summary: dict[str, dict] = {}

    for table in omics_tables:
        merged = merged.merge(
            table.data,
            how="left",
            left_on=id_column,
            right_on=table.id_column,
            suffixes=("", f"_{table.modality}"),
        )

        if table.id_column != id_column and table.id_column in merged.columns:
            merged = merged.drop(columns=[table.id_column])

        all_feature_columns.extend(table.feature_columns)
        modality_summary[table.modality] = {
            "n_features": len(table.feature_columns),
            **table.metadata,
        }

    return merged, all_feature_columns, modality_summary
