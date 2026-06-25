from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

import pandas as pd


@dataclass
class ValidationResult:
    ok: bool
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    normalized_value: Optional[str] = None
    stats: dict = field(default_factory=dict)


def validate_sequence_input(seq: str) -> ValidationResult:
    if seq is None:
        return ValidationResult(ok=False, errors=["Sequence is required."])
    normalized = str(seq).strip().upper().replace(" ", "").replace("U", "T")
    allowed = set("ACGTN")
    bad = sorted({c for c in normalized if c not in allowed})
    if not normalized:
        return ValidationResult(ok=False, errors=["Sequence is empty after normalization."])
    if bad:
        return ValidationResult(ok=False, errors=[f"Invalid characters in sequence: {bad}"])
    if len(normalized) < 8:
        return ValidationResult(ok=False, errors=["Sequence is too short; expected at least 8 bases."])
    return ValidationResult(
        ok=True,
        normalized_value=normalized,
        stats={"length": len(normalized), "gc_fraction": (normalized.count("G") + normalized.count("C")) / len(normalized)},
    )


def validate_mapping_file(path: str | Path) -> ValidationResult:
    p = Path(path)
    required = ["candidate_id"]
    optional = [
        "oligo_id",
        "registry_id",
        "manual_sequence",
        "manual_modification",
        "mapping_notes",
        "curation_status",
    ]
    if not p.exists():
        return ValidationResult(ok=False, errors=[f"Mapping file not found: {p}"])
    df = pd.read_csv(p, low_memory=False)
    missing_required = [c for c in required if c not in df.columns]
    if missing_required:
        return ValidationResult(ok=False, errors=[f"Missing required columns: {missing_required}"])
    warnings = []
    for c in optional:
        if c not in df.columns:
            warnings.append(f"Optional column missing: {c}")
    dupes = int(df["candidate_id"].astype(str).duplicated().sum())
    rows_with_any = 0
    if len(df):
        mask = False
        for c in ["oligo_id", "registry_id", "manual_sequence"]:
            if c in df.columns:
                mask = mask | df[c].notna()
        rows_with_any = int(mask.sum()) if hasattr(mask, "sum") else 0
    return ValidationResult(
        ok=True,
        warnings=warnings + ([f"Duplicate candidate_id rows: {dupes}"] if dupes else []),
        stats={"rows": int(len(df)), "rows_with_any_mapping_info": int(rows_with_any)},
    )


def validate_file_exists(path: str | Path, expected_cols: list[str] | None = None) -> ValidationResult:
    p = Path(path)
    if not p.exists():
        return ValidationResult(ok=False, errors=[f"File not found: {p}"])
    if expected_cols:
        df = pd.read_csv(p, nrows=5, low_memory=False)
        missing = [c for c in expected_cols if c not in df.columns]
        if missing:
            return ValidationResult(ok=False, errors=[f"Missing expected columns in {p}: {missing}"])
    return ValidationResult(ok=True)
