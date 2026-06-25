from __future__ import annotations

from typing import Dict, Any, Tuple, Iterable
import pandas as pd
import numpy as np

from oligotox.omics.gene_sets import GENE_SETS
from oligotox.omics.feature_dictionary import PATHWAY_FEATURES, SENTINEL_MARKERS
from oligotox.omics.mechanistic_schema import OmicsMetadata, OmicsConfidence, OmicsReadout


def _norm_effect(row: pd.Series, effect_col: str, pval_col: str | None = None) -> float:
    v = float(row.get(effect_col, 0.0) or 0.0)
    if pval_col is not None:
        p = float(row.get(pval_col, 1.0) or 1.0)
        if p < 1e-16:
            p = 1e-16
        weight = -np.log10(p) / 8.0
        weight = max(0.0, min(1.0, weight))
        v *= (0.5 + 0.5 * weight)
    return float(max(-4.0, min(4.0, v)) / 4.0)


def _score_pathway(
    deg_df: pd.DataFrame,
    geneset: set[str],
    gene_col: str,
    effect_col: str,
    pval_col: str | None = None,
) -> Tuple[float, int]:
    df = deg_df.copy()
    df["_gene_upper"] = df[gene_col].astype(str).str.upper()

    mask = df["_gene_upper"].isin(geneset)
    sub = df.loc[mask].copy()
    n = int(len(sub))
    if n == 0:
        return 0.0, 0

    vals = sub.apply(lambda row: _norm_effect(row, effect_col, pval_col), axis=1).values
    score = float(np.mean(vals))
    return score, n


def _score_sentinel(
    deg_df: pd.DataFrame,
    gene_symbol: str,
    gene_col: str,
    effect_col: str,
    pval_col: str | None = None,
) -> float:
    df = deg_df.copy()
    df["_gene_upper"] = df[gene_col].astype(str).str.upper()
    sub = df.loc[df["_gene_upper"] == gene_symbol.upper()]
    if sub.empty:
        return 0.0
    best = sub.sort_values(by=pval_col, ascending=True).iloc[0] if pval_col and pval_col in sub.columns else sub.iloc[0]
    return _norm_effect(best, effect_col, pval_col)


def _estimate_confidence(
    deg_df: pd.DataFrame,
    gene_col: str,
    effect_col: str,
    pval_col: str | None,
) -> OmicsConfidence:
    df = deg_df.copy()
    df["_gene_upper"] = df[gene_col].astype(str).str.upper()

    all_genes = set(df["_gene_upper"].unique())
    genes_in_sets = set().union(*GENE_SETS.values())
    overlap = all_genes & genes_in_sets

    coverage = len(overlap) / max(1, len(genes_in_sets))
    coverage_score = float(max(0.0, min(1.0, coverage * 2.0)))

    if pval_col and pval_col in df.columns:
        sig = df[df[pval_col] < 0.05]
        quality = float(max(0.0, min(1.0, len(sig) / max(1, len(df)) * 5.0)))
    else:
        quality = 0.5

    abs_vals = df[effect_col].astype(float).abs().values
    if len(abs_vals) > 0:
        m = float(np.median(abs_vals))
        concordance = float(max(0.0, min(1.0, m / 2.0)))
    else:
        concordance = 0.0

    translational = (coverage_score + quality) / 2.0

    return OmicsConfidence(
        coverage_score=coverage_score,
        quality_score=quality,
        concordance_score=concordance,
        translational_relevance=translational,
    )


def deg_table_to_readout(
    deg_df: pd.DataFrame,
    gene_col: str = "gene_symbol",
    effect_col: str = "log2fc",
    pval_col: str | None = "pval_adj",
    metadata: Dict[str, Any] | None = None,
) -> OmicsReadout:
    if metadata is None:
        metadata = {}

    needed_cols = [gene_col, effect_col]
    for col in needed_cols:
        if col not in deg_df.columns:
            raise ValueError(f"Required column '{col}' not found in DEG table")

    df = deg_df.copy()
    df[effect_col] = pd.to_numeric(df[effect_col], errors="coerce").fillna(0.0)
    if pval_col and pval_col in df.columns:
        df[pval_col] = pd.to_numeric(df[pval_col], errors="coerce").fillna(1.0)

    md = OmicsMetadata(
        candidate_id=metadata.get("candidate_id"),
        modality="transcriptomics",
        tissue=metadata.get("tissue"),
        species=metadata.get("species"),
        system=metadata.get("system"),
        platform=metadata.get("platform"),
        timepoint_hr=float(metadata.get("timepoint_hr", 0.0) or 0.0),
        dose_level=metadata.get("dose_level"),
        dose_value=float(metadata.get("dose_value", 0.0) or 0.0),
        dose_unit=metadata.get("dose_unit"),
        replicate_count=int(metadata.get("replicate_count", 0) or 0),
    )

    confidence = _estimate_confidence(df, gene_col=gene_col, effect_col=effect_col, pval_col=pval_col)

    pathway_scores: Dict[str, float] = {}
    coverage_counts: Dict[str, int] = {}

    for key in PATHWAY_FEATURES:
        geneset = GENE_SETS.get(key, set())
        if not geneset:
            pathway_scores[key] = 0.0
            coverage_counts[key] = 0
            continue
        score, n = _score_pathway(df, geneset, gene_col=gene_col, effect_col=effect_col, pval_col=pval_col)
        pathway_scores[key] = score
        coverage_counts[key] = n

    sentinel_markers: Dict[str, float] = {}
    for marker in SENTINEL_MARKERS:
        sentinel_markers[marker] = _score_sentinel(df, marker, gene_col=gene_col, effect_col=effect_col, pval_col=pval_col)

    raw_summary = {
        "pathway_gene_counts": coverage_counts,
        "n_genes": int(len(df)),
    }

    return OmicsReadout(
        metadata=md,
        pathway_scores=pathway_scores,
        sentinel_markers=sentinel_markers,
        confidence=confidence,
        raw_summary=raw_summary,
    )
