from __future__ import annotations

from pathlib import Path
import pandas as pd


AGG_NUMERIC_DEFAULT = "mean"


def build_candidate_level_omics_table(
    omics_training_table_path: str | Path,
    output_path: str | Path,
    id_column: str = "candidate_id",
) -> pd.DataFrame:
    path = Path(omics_training_table_path)
    if not path.exists():
        raise FileNotFoundError(f"Missing omics training table: {path}")

    df = pd.read_csv(path)
    if df.empty:
        out = pd.DataFrame(columns=[id_column])
        out.to_csv(output_path, index=False)
        return out

    numeric_cols = [c for c in df.columns if c != id_column and pd.api.types.is_numeric_dtype(df[c])]
    agg = {c: AGG_NUMERIC_DEFAULT for c in numeric_cols}

    out = df.groupby(id_column, as_index=False).agg(agg)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output_path, index=False)
    return out
