from __future__ import annotations

import re
import pandas as pd

from oligotox.omics.reference.gene_aliases import GENE_ALIAS_MAP


def clean_gene_symbol(x: str) -> str:
    s = str(x).strip().upper()
    s = re.sub(r"\s+", "", s)
    s = s.replace("/", "")
    s = s.replace("\\", "")
    s = s.replace("-", "")
    s = GENE_ALIAS_MAP.get(s, s)
    return s


def normalize_gene_symbols(df: pd.DataFrame, gene_col: str = "gene_symbol") -> pd.DataFrame:
    df = df.copy()
    if gene_col not in df.columns:
        raise ValueError(f"Missing gene column: {gene_col}")
    df[gene_col] = df[gene_col].astype(str).map(clean_gene_symbol)
    df = df[df[gene_col].notna() & (df[gene_col] != "")]
    return df
