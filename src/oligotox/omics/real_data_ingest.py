from __future__ import annotations

from pathlib import Path
from typing import Any
import json
import pandas as pd

from oligotox.omics.column_mapper import standardize_deg_columns
from oligotox.omics.symbol_cleaner import normalize_gene_symbols
from oligotox.omics.deg_qc import qc_deg_table
from oligotox.omics.deg_to_readout import deg_table_to_readout


def _require_any_effect_column(df: pd.DataFrame) -> str:
    """
    Return the name of the numeric effect column to use.
    Accepts either 'log2fc' or 'log_fc', preferring 'log2fc' when both exist.
    """
    if "log2fc" in df.columns:
        return "log2fc"
    if "log_fc" in df.columns:
        return "log_fc"
    raise ValueError(
        "Missing required standardized DEG effect column. "
        "Expected one of ['log2fc', 'log_fc'], "
        f"got columns: {list(df.columns)}"
    )


def _aggregate_duplicate_genes(
    df: pd.DataFrame,
    gene_col: str,
    effect_col: str,
    pval_adj_col: str | None = "pval_adj",
) -> pd.DataFrame:
    df = df.copy()
    df[effect_col] = pd.to_numeric(df[effect_col], errors="coerce").fillna(0.0)

    if pval_adj_col and pval_adj_col in df.columns:
        df[pval_adj_col] = pd.to_numeric(df[pval_adj_col], errors="coerce").fillna(1.0)
        df = df.sort_values([gene_col, pval_adj_col, effect_col], ascending=[True, True, False])
    else:
        df = df.sort_values([gene_col, effect_col], ascending=[True, False])

    df = df.drop_duplicates(subset=[gene_col], keep="first")
    return df


def ingest_deg_file_to_readout(
    file_path: str | Path,
    metadata: dict[str, Any],
    qc_profile: str = "standard",
    fail_on_qc: bool = True,
) -> tuple[object, dict]:
    file_path = Path(file_path)

    if file_path.suffix.lower() in [".tsv", ".txt"]:
        df = pd.read_csv(file_path, sep="\t")
    else:
        df = pd.read_csv(file_path)

    # Standardize columns and gene symbols
    df = standardize_deg_columns(df)
    if "gene_symbol" not in df.columns:
        raise ValueError(
            f"Missing required standardized DEG column 'gene_symbol'. "
            f"Available columns after mapping: {list(df.columns)}"
        )

    effect_col = _require_any_effect_column(df)

    df = normalize_gene_symbols(df, gene_col="gene_symbol")
    df = _aggregate_duplicate_genes(
        df,
        gene_col="gene_symbol",
        effect_col=effect_col,
        pval_adj_col="pval_adj",
    )

    # QC
    qc = qc_deg_table(
        df,
        gene_col="gene_symbol",
        effect_col=effect_col,
        pval_adj_col="pval_adj",
        qc_profile=qc_profile,
    )

    if not qc["passes_basic_qc"] and fail_on_qc:
        raise ValueError(f"DEG file failed basic QC: {json.dumps(qc, indent=2)}")

    # Always pass the actual effect_col into deg_table_to_readout
    readout = deg_table_to_readout(
        df,
        gene_col="gene_symbol",
        effect_col=effect_col,
        pval_col="pval_adj" if "pval_adj" in df.columns else None,
        metadata=metadata,
    )

    readout.raw_summary["qc"] = qc
    return readout, qc


def readout_to_payload(readout) -> dict:
    payload = {
        "candidate_id": readout.metadata.candidate_id,
        "modality": readout.metadata.modality,
        "tissue": readout.metadata.tissue,
        "species": readout.metadata.species,
        "system": readout.metadata.system,
        "platform": readout.metadata.platform,
        "timepoint_hr": readout.metadata.timepoint_hr,
        "dose_level": readout.metadata.dose_level,
        "dose_value": readout.metadata.dose_value,
        "dose_unit": readout.metadata.dose_unit,
        "replicate_count": readout.metadata.replicate_count,
        "coverage_score": readout.confidence.coverage_score,
        "quality_score": readout.confidence.quality_score,
        "concordance_score": readout.confidence.concordance_score,
        "translational_relevance": readout.confidence.translational_relevance,
    }
    for k, v in readout.pathway_scores.items():
        payload[f"pw_{k}"] = v
    for k, v in readout.sentinel_markers.items():
        payload[f"sentinel_{k}"] = v
    return payload
