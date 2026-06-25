from __future__ import annotations
import ast
import math
from typing import Any, Dict, Iterable, Optional

import pandas as pd

from oligotox.integrations_v03.feature_matrix import build_feature_row

V03_FEATURE_PREFIXES = ("chem_", "d2_", "struct_", "ctx_")

def _missing(x: Any) -> bool:
    try:
        return x is None or (isinstance(x, float) and math.isnan(x))
    except Exception:
        return x is None

def _coerce_listlike(x: Any):
    if _missing(x) or x == "":
        return []
    if isinstance(x, list):
        return x
    if isinstance(x, tuple):
        return list(x)
    if isinstance(x, str):
        s = x.strip()
        if not s:
            return []
        try:
            parsed = ast.literal_eval(s)
            if isinstance(parsed, list):
                return parsed
            if isinstance(parsed, tuple):
                return list(parsed)
        except Exception:
            pass
        if "," in s:
            return [p.strip() for p in s.split(",") if p.strip()]
        return [s]
    return []

def _pick(row: pd.Series, *names: str, default=None):
    for n in names:
        if n in row and not _missing(row[n]):
            return row[n]
    return default

def canonical_candidate_from_row(row: pd.Series) -> Dict[str, Any]:
    seq = _pick(row, "sequence", "seq", "oligo_sequence", "guide_sequence", default="")
    delivery = _pick(row, "delivery", "delivery_system", "delivery_modality", default="")
    route = _pick(row, "route", "administration_route", default="")
    modality = _pick(row, "modality", "oligo_type", "class", default="")
    ps_fraction = _pick(row, "ps_fraction", default=0.0)
    mod_string = _pick(row, "modification_string", "mods", "mod_string", default=None)
    lna_positions = _pick(row, "lna_positions", default=[])
    sugar_mods = _pick(row, "sugar_mods", default=[])

    return {
        "sequence": "" if _missing(seq) else str(seq),
        "delivery": "" if _missing(delivery) else str(delivery),
        "route": "" if _missing(route) else str(route),
        "modality": "" if _missing(modality) else str(modality),
        "ps_fraction": 0.0 if _missing(ps_fraction) else float(ps_fraction),
        "modification_string": None if _missing(mod_string) else str(mod_string),
        "lna_positions": _coerce_listlike(lna_positions),
        "sugar_mods": _coerce_listlike(sugar_mods),
    }

def build_v03_feature_df(df: pd.DataFrame, use_viennarna: bool = False, use_embeddings: bool = False) -> pd.DataFrame:
    rows = []
    for _, row in df.iterrows():
        cand = canonical_candidate_from_row(row)
        feats = build_feature_row(
            cand,
            use_viennarna=use_viennarna,
            use_embeddings=use_embeddings,
        )
        rows.append(feats)
    out = pd.DataFrame(rows).fillna(0.0)
    v03_cols = [c for c in out.columns if c.startswith(V03_FEATURE_PREFIXES)]
    return out[v03_cols].copy()

def merge_v03_into_canonical(
    input_csv: str = "data/processed/feature_matrix.csv",
    output_csv: str = "data/processed/feature_matrix.csv",
    backup_csv: Optional[str] = None,
    use_viennarna: bool = False,
    use_embeddings: bool = False,
):
    df = pd.read_csv(input_csv, low_memory=False)
    v03 = build_v03_feature_df(df, use_viennarna=use_viennarna, use_embeddings=use_embeddings)

    overlap = [c for c in v03.columns if c in df.columns]
    if overlap:
        df = df.drop(columns=overlap)

    merged = pd.concat([df.reset_index(drop=True), v03.reset_index(drop=True)], axis=1)

    if backup_csv:
        merged.to_csv(backup_csv, index=False)
    merged.to_csv(output_csv, index=False)

    return {
        "rows": int(merged.shape[0]),
        "cols": int(merged.shape[1]),
        "new_v03_cols": int(v03.shape[1]),
        "v03_colnames": list(v03.columns),
    }
