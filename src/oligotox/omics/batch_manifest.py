from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import pandas as pd


REQUIRED_COLUMNS = [
    "candidate_id",
    "study_id",
    "comparison_label",
    "file_path",
    "tissue",
    "species",
    "system",
    "platform",
    "timepoint_hr",
    "dose_level",
    "dose_value",
    "dose_unit",
    "replicate_count",
    "qc_profile",
]


@dataclass
class ManifestRow:
    candidate_id: str
    study_id: str
    comparison_label: str
    file_path: str
    tissue: str
    species: str
    system: str
    platform: str
    timepoint_hr: float
    dose_level: str
    dose_value: float
    dose_unit: str
    replicate_count: int
    qc_profile: str


def load_manifest(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Manifest not found: {path}")

    df = pd.read_csv(path)
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"Manifest missing required columns: {missing}")

    return df


def manifest_to_rows(df: pd.DataFrame) -> list[ManifestRow]:
    rows = []
    for _, r in df.iterrows():
        rows.append(
            ManifestRow(
                candidate_id=str(r["candidate_id"]),
                study_id=str(r["study_id"]),
                comparison_label=str(r["comparison_label"]),
                file_path=str(r["file_path"]),
                tissue=str(r["tissue"]),
                species=str(r["species"]),
                system=str(r["system"]),
                platform=str(r["platform"]),
                timepoint_hr=float(r["timepoint_hr"]),
                dose_level=str(r["dose_level"]),
                dose_value=float(r["dose_value"]),
                dose_unit=str(r["dose_unit"]),
                replicate_count=int(r["replicate_count"]),
                qc_profile=str(r["qc_profile"]),
            )
        )
    return rows
