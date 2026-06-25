from __future__ import annotations

import re
import pandas as pd


COLUMN_CANDIDATES = {
    "gene_symbol": [
        "gene",
        "gene_symbol",
        "symbol",
        "genesymbol",
        "gene_name",
        "hgnc_symbol",
        "external_gene_name",
    ],
    "log2fc": [
        "log2fc",
        "logfc",
        "log_fold_change",
        "log2_fold_change",
        "log2foldchange",
        "logfoldchange",
        "lfc",
        "avg_log2fc",
        "avg_logfc",
        "effect",
        "estimate",
    ],
    "pval_adj": [
        "padj",
        "adj_p",
        "adj_pval",
        "adjusted_p_value",
        "adjpvalue",
        "adj_p_value",
        "fdr",
        "qvalue",
        "pval_adj",
        "p_adj",
        "fdr_bh",
        "adj_p_val",
        "adj.p.val",
    ],
    "pval": [
        "pvalue",
        "p_value",
        "pval",
        "p",
        "p.value",
    ],
    "base_mean": [
        "basemean",
        "base_mean",
        "mean_expr",
        "avg_expr",
        "expression_mean",
    ],
}


def normalize_name(x: str) -> str:
    s = str(x).strip()
    s = re.sub(r"([a-z0-9])([A-Z])", r"\1_\2", s)
    s = s.lower()
    s = s.replace(".", "_")
    s = s.replace("-", "_")
    s = s.replace("/", "_")
    s = re.sub(r"\s+", "_", s)
    s = re.sub(r"_+", "_", s)
    s = s.strip("_")
    return s


def standardize_deg_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()

    norm_cols = {c: normalize_name(c) for c in df.columns}
    df = df.rename(columns=norm_cols)

    reverse_lookup = {}
    for target, candidates in COLUMN_CANDIDATES.items():
        for c in candidates:
            reverse_lookup[normalize_name(c)] = target

    rename = {}
    for c in df.columns:
        if c in reverse_lookup:
            rename[c] = reverse_lookup[c]

    df = df.rename(columns=rename)
    return df
