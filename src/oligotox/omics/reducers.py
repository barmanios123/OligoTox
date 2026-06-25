from __future__ import annotations

import pandas as pd


def passthrough(df: pd.DataFrame, id_column: str) -> pd.DataFrame:
    return df.copy()


def variance_topk(df: pd.DataFrame, id_column: str, k: int = 50) -> pd.DataFrame:
    if id_column not in df.columns:
        raise ValueError(f"Missing id column: {id_column}")

    feature_cols = [c for c in df.columns if c != id_column]
    if not feature_cols:
        return df[[id_column]].copy()

    numeric = df[feature_cols].apply(pd.to_numeric, errors="coerce")
    variances = numeric.var(axis=0, skipna=True).sort_values(ascending=False)
    keep = list(variances.head(min(k, len(variances))).index)

    out = pd.concat([df[[id_column]].copy(), numeric[keep]], axis=1)
    return out
