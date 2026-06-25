from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


@dataclass
class OmicsModalitySpec:
    name: str
    required_id_column: str = "candidate_id"
    file_name: str | None = None
    prefix: str | None = None
    reduce_method: str = "passthrough"
    max_features: int | None = None
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass
class OmicsTable:
    modality: str
    id_column: str
    feature_columns: list[str]
    data: Any
    metadata: dict[str, Any] = field(default_factory=dict)
