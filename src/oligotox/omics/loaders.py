from __future__ import annotations

from pathlib import Path
import pandas as pd

from oligotox.omics.schemas import OmicsModalitySpec, OmicsTable
from oligotox.omics.reducers import passthrough, variance_topk


def _apply_prefix(df: pd.DataFrame, id_column: str, prefix: str | None) -> pd.DataFrame:
    if not prefix:
        return df
    rename_map = {c: f"{prefix}{c}" for c in df.columns if c != id_column}
    return df.rename(columns=rename_map)


def _reduce(df: pd.DataFrame, spec: OmicsModalitySpec) -> pd.DataFrame:
    if spec.reduce_method == "passthrough":
        return passthrough(df, spec.required_id_column)
    if spec.reduce_method == "variance_topk":
        return variance_topk(df, spec.required_id_column, k=spec.max_features or 50)
    raise ValueError(f"Unknown reduce_method={spec.reduce_method} for modality={spec.name}")


def load_omics_table(raw_dir: str | Path, spec: OmicsModalitySpec) -> OmicsTable | None:
    raw_dir = Path(raw_dir)
    if spec.file_name is None:
        return None

    path = raw_dir / spec.file_name
    if not path.exists():
        return None

    df = pd.read_csv(path)
    if spec.required_id_column not in df.columns:
        raise ValueError(
            f"{spec.name}: required id column '{spec.required_id_column}' not found in {path}"
        )

    df = _reduce(df, spec)
    df = _apply_prefix(df, spec.required_id_column, spec.prefix)

    feature_columns = [c for c in df.columns if c != spec.required_id_column]
    return OmicsTable(
        modality=spec.name,
        id_column=spec.required_id_column,
        feature_columns=feature_columns,
        data=df,
        metadata={
            "source_file": str(path),
            "n_rows": int(len(df)),
            "n_features": int(len(feature_columns)),
        },
    )
