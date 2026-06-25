from __future__ import annotations

import pandas as pd


QC_PROFILES = {
    "strict": {
        "min_rows": 100,
        "min_unique_genes": 100,
        "max_missing_gene_fraction": 0.10,
        "max_missing_effect_fraction": 0.10,
    },
    "standard": {
        "min_rows": 50,
        "min_unique_genes": 50,
        "max_missing_gene_fraction": 0.20,
        "max_missing_effect_fraction": 0.20,
    },
    "relaxed": {
        "min_rows": 10,
        "min_unique_genes": 10,
        "max_missing_gene_fraction": 0.30,
        "max_missing_effect_fraction": 0.30,
    },
}


def qc_deg_table(
    df: pd.DataFrame,
    gene_col: str = "gene_symbol",
    effect_col: str = "log2fc",
    pval_adj_col: str | None = "pval_adj",
    qc_profile: str = "standard",
) -> dict:
    if qc_profile not in QC_PROFILES:
        raise ValueError(f"Unknown qc_profile='{qc_profile}'. Expected one of {list(QC_PROFILES)}")

    cfg = QC_PROFILES[qc_profile]

    out = {
        "qc_profile": qc_profile,
        "n_rows": int(len(df)),
        "n_unique_genes": 0,
        "missing_gene_fraction": 0.0,
        "missing_effect_fraction": 0.0,
        "has_pval_adj": False,
        "significant_gene_count": 0,
    }

    if gene_col in df.columns:
        out["n_unique_genes"] = int(df[gene_col].astype(str).nunique())
        out["missing_gene_fraction"] = float(df[gene_col].isna().mean())

    if effect_col in df.columns:
        out["missing_effect_fraction"] = float(pd.to_numeric(df[effect_col], errors="coerce").isna().mean())

    if pval_adj_col and pval_adj_col in df.columns:
        out["has_pval_adj"] = True
        p = pd.to_numeric(df[pval_adj_col], errors="coerce")
        out["significant_gene_count"] = int((p < 0.05).fillna(False).sum())

    out["passes_basic_qc"] = bool(
        out["n_rows"] >= cfg["min_rows"] and
        out["n_unique_genes"] >= cfg["min_unique_genes"] and
        out["missing_gene_fraction"] <= cfg["max_missing_gene_fraction"] and
        out["missing_effect_fraction"] <= cfg["max_missing_effect_fraction"]
    )
    return out
