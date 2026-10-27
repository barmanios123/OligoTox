from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from oligotox.config.paths import RAW_DATA_DIR
from oligotox.omics import OMICS_REGISTRY, load_omics_table, join_omics_tables
from oligotox.omics.join_utils import safe_left_join_on_candidate
from oligotox.features.aso_rich_features import build_aso_rich_features
from oligotox.features.candidate import recompute_sequence_features

try:
    from features.aso_optimizer_v03 import batch_aso_features as batch_aso_features_v03
except Exception:
    batch_aso_features_v03 = None

from oligotox.features.structure_features import build_structure_features
from oligotox.features.embedding_features import build_embedding_features

try:
    from features.sequence_embeddings_v03 import SequenceEmbedder as SequenceEmbedder_v03
except Exception:
    SequenceEmbedder_v03 = None


ID_COLUMNS = ["candidate_id", "source", "notes", "series_id"]

LABEL_COLUMNS = [
    "hepatocellular_injury",
    "hepatic_vacuolation",
    "renal_tubular_injury",
    "renal_accumulation",
    "thrombocytopenia",
    "coagulation_effect",
    "complement_activation",
    "innate_immune_activation",
    "cytokine_release",
    "kupffer_cell_activation",
    "cns_toxicity",
    "offtarget_gene_silencing",
    "injection_site_reaction",
]

CURATED_ENDPOINT_MAP = {
    "hepatocellular_injury": "endpoint_hepatocellular_injury",
    "hepatic_vacuolation": "endpoint_hepatic_vacuolation",
    "renal_tubular_injury": "endpoint_renal_tubular_injury",
    "renal_accumulation": "endpoint_renal_accumulation",
    "thrombocytopenia": "endpoint_thrombocytopenia",
    "coagulation_effect": "endpoint_coagulation_effect",
    "complement_activation": "endpoint_complement_activation",
    "innate_immune_activation": "endpoint_innate_immune_activation",
    "cytokine_release": "endpoint_cytokine_release",
    "kupffer_cell_activation": "endpoint_kupffer_cell_activation",
    "cns_toxicity": "endpoint_cns_neurotoxicity",
    "offtarget_gene_silencing": "endpoint_offtarget_gene_silencing",
    "injection_site_reaction": "endpoint_injection_site_reaction",
}


def _json_dump(path: str | Path, obj: dict[str, Any]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2))


def _load_base_training_table(training_csv_path: str | Path) -> pd.DataFrame:
    training_csv_path = Path(training_csv_path)

    if training_csv_path.exists():
        return pd.read_csv(training_csv_path)

    raise FileNotFoundError(
        f"Base training dataset not found at {training_csv_path}. "
        "Provide data/processed/training_dataset.csv or update the loader."
    )


def _detect_feature_columns(df: pd.DataFrame) -> list[str]:
    protected = set(ID_COLUMNS + LABEL_COLUMNS + list(CURATED_ENDPOINT_MAP.values()))
    protected.discard("sequence")
    return [c for c in df.columns if c not in protected]


def _maybe_join_candidate_level_omics(
    base_df: pd.DataFrame,
    omics_features_path: str | Path = "data/processed/omics/candidate_level_omics_features.csv",
    id_column: str = "candidate_id",
) -> tuple[pd.DataFrame, list[str]]:
    """
    Join candidate-level omics features produced by the batch DEG ingestion pipeline.

    Returns:
      - updated dataframe
      - list of added omics column names
    """
    before_cols = set(base_df.columns)

    joined = safe_left_join_on_candidate(
        base_df=base_df,
        omics_path=omics_features_path,
        id_column=id_column,
        prefix="omics_",
    )

    added_cols = [c for c in joined.columns if c not in before_cols and c.startswith("omics_")]
    return joined, added_cols


def _normalize_label_series(series: pd.Series) -> pd.Series:
    s = pd.to_numeric(series, errors="coerce")
    return s.map(lambda x: 1.0 if x == 1 else (0.0 if x == 0 else np.nan))


def _coerce_feature_output(
    obj: Any,
    n_rows: int,
    default_index: pd.Index,
) -> pd.DataFrame:
    if obj is None:
        return pd.DataFrame(index=default_index)

    if isinstance(obj, pd.DataFrame):
        df = obj.copy()
    elif isinstance(obj, list):
        try:
            df = pd.DataFrame(obj)
        except Exception:
            return pd.DataFrame(index=default_index)
    elif isinstance(obj, dict):
        try:
            df = pd.DataFrame(obj)
        except Exception:
            return pd.DataFrame(index=default_index)
    else:
        try:
            df = pd.DataFrame(obj)
        except Exception:
            return pd.DataFrame(index=default_index)

    if df.empty:
        return pd.DataFrame(index=default_index)

    df = df.reset_index(drop=True)

    if len(df) != n_rows:
        return pd.DataFrame(index=default_index)

    return df


def run(
    training_csv_path: str | Path | None = None,
    feature_matrix_path: str | Path | None = None,
    label_matrix_path: str | Path | None = None,
    dataset_metadata_path: str | Path | None = None,
    training_dataset_path: str | Path | None = None,
    metadata_path: str | Path | None = None,
) -> dict[str, Any]:
    """
    Build feature/label matrices and dataset metadata for training.

    Backward-compatible aliases:
    - training_dataset_path -> training_csv_path
    - metadata_path -> dataset_metadata_path
    """
    if training_csv_path is None:
        training_csv_path = training_dataset_path
    if dataset_metadata_path is None:
        dataset_metadata_path = metadata_path

    if training_csv_path is None:
        raise ValueError("training_csv_path (or training_dataset_path) is required")
    if feature_matrix_path is None:
        raise ValueError("feature_matrix_path is required")
    if label_matrix_path is None:
        raise ValueError("label_matrix_path is required")
    if dataset_metadata_path is None:
        raise ValueError("dataset_metadata_path (or metadata_path) is required")

    training_csv_path = Path(training_csv_path)
    feature_matrix_path = Path(feature_matrix_path)
    label_matrix_path = Path(label_matrix_path)
    dataset_metadata_path = Path(dataset_metadata_path)

    # 1) Load base training table
    df = _load_base_training_table(training_csv_path)
    if "sequence" in df.columns:
        df = recompute_sequence_features(df)

    if "candidate_id" not in df.columns:
        raise ValueError("training dataset must contain candidate_id")

    # 2) Attach v0.3 ASO chemistry features at the base level (df)
    aso_chemistry_features_v03 = pd.DataFrame(index=df.index)
    if batch_aso_features_v03 is not None:
        try:
            raw_aso = batch_aso_features_v03(df)
            aso_chemistry_features_v03 = _coerce_feature_output(
                raw_aso,
                n_rows=len(df),
                default_index=df.index,
            )
        except Exception:
            aso_chemistry_features_v03 = pd.DataFrame(index=df.index)

    if not aso_chemistry_features_v03.empty:
        overlap = [c for c in aso_chemistry_features_v03.columns if c in df.columns]
        if overlap:
            aso_chemistry_features_v03 = aso_chemistry_features_v03.drop(columns=overlap)
        df = pd.concat(
            [df.reset_index(drop=True), aso_chemistry_features_v03.reset_index(drop=True)],
            axis=1,
        )

    # 3) Omics joins (use df with chem features)
    available_omics = []
    for _, spec in OMICS_REGISTRY.items():
        table = load_omics_table(RAW_DATA_DIR, spec)
        if table is not None:
            available_omics.append(table)

    df_joined, registry_omics_feature_columns, omics_summary = join_omics_tables(
        base_df=df,
        omics_tables=available_omics,
        id_column="candidate_id",
    )

    df_joined, candidate_level_omics_columns = _maybe_join_candidate_level_omics(
        base_df=df_joined,
        omics_features_path="data/processed/omics/candidate_level_omics_features.csv",
        id_column="candidate_id",
    )

    # 4) ASO-rich sequence features
    rich_aso_features = build_aso_rich_features(df_joined)
    rich_aso_features = _coerce_feature_output(
        rich_aso_features,
        n_rows=len(df_joined),
        default_index=df_joined.index,
    )
    if not rich_aso_features.empty:
        overlap = [c for c in rich_aso_features.columns if c in df_joined.columns]
        if overlap:
            rich_aso_features = rich_aso_features.drop(columns=overlap)
        df_joined = pd.concat(
            [df_joined.reset_index(drop=True), rich_aso_features.reset_index(drop=True)],
            axis=1,
        )

    # 5) ViennaRNA-based structure features
    structure_features = build_structure_features(
        df_joined,
        cache_path="data/processed/structure_features_cache.csv",
        use_viennarna=True,
    )
    structure_features = _coerce_feature_output(
        structure_features,
        n_rows=len(df_joined),
        default_index=df_joined.index,
    )
    if not structure_features.empty:
        overlap = [c for c in structure_features.columns if c in df_joined.columns]
        if overlap:
            structure_features = structure_features.drop(columns=overlap)
        df_joined = pd.concat(
            [df_joined.reset_index(drop=True), structure_features.reset_index(drop=True)],
            axis=1,
        )

    # 6) Baseline embedding_features (stub_hash_rna)
    embedding_features = build_embedding_features(
        df_joined,
        cache_path="data/processed/embedding_features_cache.csv",
        model_name="stub_hash_rna",
        dim=16,
    )
    embedding_features = _coerce_feature_output(
        embedding_features,
        n_rows=len(df_joined),
        default_index=df_joined.index,
    )
    if not embedding_features.empty:
        overlap = [c for c in embedding_features.columns if c in df_joined.columns]
        if overlap:
            embedding_features = embedding_features.drop(columns=overlap)
        df_joined = pd.concat(
            [df_joined.reset_index(drop=True), embedding_features.reset_index(drop=True)],
            axis=1,
        )

    # 7) v0.3 sequence embeddings: disabled until we align API (no .transform in local SequenceEmbedder)
    sequence_embedding_features_v03 = pd.DataFrame(index=df_joined.index)

    # 8) Labels
    for label_col, curated_col in CURATED_ENDPOINT_MAP.items():
        if curated_col in df_joined.columns:
            df_joined[label_col] = _normalize_label_series(df_joined[curated_col])
        elif label_col in df_joined.columns:
            df_joined[label_col] = _normalize_label_series(df_joined[label_col])
        else:
            df_joined[label_col] = np.nan

    feature_columns = _detect_feature_columns(df_joined)
    missing_labels = [c for c in LABEL_COLUMNS if c not in df_joined.columns]
    if missing_labels:
        raise ValueError(f"Missing expected label columns: {missing_labels}")

    X = df_joined[feature_columns].copy()
    Y = df_joined[LABEL_COLUMNS].copy()

    feature_matrix_path.parent.mkdir(parents=True, exist_ok=True)
    label_matrix_path.parent.mkdir(parents=True, exist_ok=True)
    dataset_metadata_path.parent.mkdir(parents=True, exist_ok=True)

    X.to_csv(feature_matrix_path, index=False)
    Y.to_csv(label_matrix_path, index=False)

    all_omics_feature_columns = sorted(
        set(registry_omics_feature_columns) | set(candidate_level_omics_columns)
    )

    metadata = {
        "training_csv": str(training_csv_path),
        "rows": int(len(df_joined)),
        "n_features": int(len(feature_columns)),
        "n_labels": int(len(LABEL_COLUMNS)),
        "feature_columns": feature_columns,
        "label_columns": LABEL_COLUMNS,
        "id_columns": [c for c in ID_COLUMNS if c in df_joined.columns],
        "omics_modalities_loaded": list(omics_summary.keys()),
        "omics_feature_columns": all_omics_feature_columns,
        "omics_feature_columns_registry": registry_omics_feature_columns,
        "omics_feature_columns_candidate_level": candidate_level_omics_columns,
        "omics_features_added": int(len(candidate_level_omics_columns)),
        "omics_summary": omics_summary,
        "candidate_level_omics_joined": bool(len(candidate_level_omics_columns) > 0),
        "candidate_level_omics_source": "data/processed/omics/candidate_level_omics_features.csv",
        "aso_rich_feature_columns": list(rich_aso_features.columns) if not rich_aso_features.empty else [],
        "aso_rich_features_added": int(rich_aso_features.shape[1]) if not rich_aso_features.empty else 0,
        "aso_optimizer_v03_feature_columns": list(aso_chemistry_features_v03.columns) if not aso_chemistry_features_v03.empty else [],
        "aso_optimizer_v03_features_added": int(aso_chemistry_features_v03.shape[1]) if not aso_chemistry_features_v03.empty else 0,
        "structure_feature_columns": list(structure_features.columns) if not structure_features.empty else [],
        "structure_features_added": int(structure_features.shape[1]) if not structure_features.empty else 0,
        "embedding_feature_columns": list(embedding_features.columns) if not embedding_features.empty else [],
        "embedding_features_added": int(embedding_features.shape[1]) if not embedding_features.empty else 0,
        "embedding_v03_feature_columns": list(sequence_embedding_features_v03.columns) if not sequence_embedding_features_v03.empty else [],
        "embedding_v03_features_added": int(sequence_embedding_features_v03.shape[1]) if not sequence_embedding_features_v03.empty else 0,
        "viennarna_requested": True,
    }
    _json_dump(dataset_metadata_path, metadata)

    return metadata