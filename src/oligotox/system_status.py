from __future__ import annotations
import json
from pathlib import Path
from typing import Dict, Any
import pandas as pd

def _safe_read_csv(path: Path):
    if not path.exists():
        return None
    try:
        return pd.read_csv(path, low_memory=False)
    except Exception:
        return None

def get_system_status(repo_root=".") -> Dict[str, Any]:
    root = Path(repo_root)
    training_path = root / "data/processed/training_dataset.csv"
    with_chem_path = root / "data/processed/training_dataset_with_chem.csv"
    audit_path = root / "data/processed/training_dataset_rebuilt_audit.csv"
    mapping_path = root / "data/metadata/candidate_to_oligo_map.csv"
    training_df = _safe_read_csv(training_path)
    with_chem_df = _safe_read_csv(with_chem_path)
    audit_df = _safe_read_csv(audit_path)
    mapping_df = _safe_read_csv(mapping_path)
    total_rows = None
    rows_with_sequence = None
    rows_missing_sequence = None
    if training_df is not None and "sequence" in training_df.columns:
        seq = training_df["sequence"].astype(str).str.strip().replace({"": pd.NA, "nan": pd.NA, "None": pd.NA})
        total_rows = int(len(training_df))
        rows_with_sequence = int(seq.notna().sum())
        rows_missing_sequence = int(total_rows - rows_with_sequence)
    active_sequence_subset_rows = int(len(with_chem_df)) if with_chem_df is not None else None
    mapping_rows = int(len(mapping_df)) if mapping_df is not None else None
    labeled_rows = None
    if audit_df is not None and "endpoint" in audit_df.columns and "n_labeled" in audit_df.columns:
        try:
            labeled_rows = int(audit_df["n_labeled"].fillna(0).max())
        except Exception:
            labeled_rows = None
    trust_banner = None
    if total_rows is not None and rows_with_sequence is not None:
        trust_banner = (f"Pipeline operational, but only {rows_with_sequence} of {total_rows} training rows currently have valid sequence. "
                        f"Usable sequence-aware model subset remains limited.")
    return {
        "training_dataset_path": str(training_path),
        "with_chem_dataset_path": str(with_chem_path),
        "audit_path": str(audit_path),
        "mapping_path": str(mapping_path),
        "total_rows": total_rows,
        "rows_with_sequence": rows_with_sequence,
        "rows_missing_sequence": rows_missing_sequence,
        "active_sequence_subset_rows": active_sequence_subset_rows,
        "mapping_rows": mapping_rows,
        "labeled_rows": labeled_rows,
        "trust_banner": trust_banner,
    }
