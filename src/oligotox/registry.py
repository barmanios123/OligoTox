from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> Optional[str]:
    if not path.exists() or not path.is_file():
        return None
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


@dataclass
class ArtifactRecord:
    name: str
    path: str
    kind: str
    checksum_sha256: Optional[str] = None
    size_bytes: Optional[int] = None


@dataclass
class ModelRecord:
    model_name: str
    endpoint: str
    artifact_path: str
    metrics: Dict[str, Any] = field(default_factory=dict)
    feature_schema_version: Optional[str] = None
    training_dataset_hash: Optional[str] = None
    created_at: str = field(default_factory=utc_now_iso)
    source_run_id: Optional[str] = None
    is_active: bool = False
    artifact_kind: str = "model_binary"


@dataclass
class PipelineRun:
    run_id: str
    job_name: str
    status: str
    trigger: str = "manual"
    started_at: str = field(default_factory=utc_now_iso)
    finished_at: Optional[str] = None
    parameters: Dict[str, Any] = field(default_factory=dict)
    logs_path: Optional[str] = None
    artifacts: List[ArtifactRecord] = field(default_factory=list)
    metrics: Dict[str, Any] = field(default_factory=dict)
    error: Optional[str] = None


class Registry:
    def __init__(self, base_dir: str | Path = "artifacts"):
        self.base_dir = Path(base_dir)
        self.runs_dir = self.base_dir / "runs"
        self.models_dir = self.base_dir / "models"
        self.registry_dir = self.base_dir / "registry"
        self.runs_dir.mkdir(parents=True, exist_ok=True)
        self.models_dir.mkdir(parents=True, exist_ok=True)
        self.registry_dir.mkdir(parents=True, exist_ok=True)
        self.runs_index_path = self.registry_dir / "runs_index.json"
        self.models_index_path = self.registry_dir / "models_index.json"

    def _load_index(self, path: Path) -> list[dict]:
        if not path.exists():
            return []
        try:
            payload = json.loads(path.read_text())
            return payload if isinstance(payload, list) else []
        except Exception:
            return []

    def _save_index(self, path: Path, rows: list[dict]) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(rows, indent=2, default=str))

    def save_run(self, run: PipelineRun) -> Path:
        run_dir = self.runs_dir / run.run_id
        run_dir.mkdir(parents=True, exist_ok=True)
        run_path = run_dir / "run.json"
        run_path.write_text(json.dumps(asdict(run), indent=2, default=str))
        rows = [r for r in self._load_index(self.runs_index_path) if r.get("run_id") != run.run_id]
        rows.append(asdict(run))
        self._save_index(self.runs_index_path, rows)
        return run_path

    def load_run(self, run_id: str) -> Optional[dict]:
        run_path = self.runs_dir / run_id / "run.json"
        if not run_path.exists():
            return None
        try:
            return json.loads(run_path.read_text())
        except Exception:
            return None

    def list_runs(self) -> list[dict]:
        rows = self._load_index(self.runs_index_path)
        if not rows:
            for path in sorted(self.runs_dir.glob("*/run.json")):
                try:
                    rows.append(json.loads(path.read_text()))
                except Exception:
                    pass
            if rows:
                self._save_index(self.runs_index_path, rows)
        return sorted(rows, key=lambda x: x.get("started_at", ""), reverse=True)

    def register_artifact(self, path: str | Path, kind: str, name: Optional[str] = None) -> ArtifactRecord:
        p = Path(path)
        return ArtifactRecord(
            name=name or p.name,
            path=str(p),
            kind=kind,
            checksum_sha256=sha256_file(p),
            size_bytes=p.stat().st_size if p.exists() and p.is_file() else None,
        )

    def register_model(self, record: ModelRecord) -> None:
        rows = self._load_index(self.models_index_path)
        if record.is_active:
            for r in rows:
                if r.get("artifact_kind") == "model_binary":
                    r["is_active"] = False
        rows = [r for r in rows if r.get("model_name") != record.model_name]
        rows.append(asdict(record))
        self._save_index(self.models_index_path, rows)

    def list_models(self) -> list[dict]:
        rows = self._load_index(self.models_index_path)
        return sorted(rows, key=lambda x: x.get("created_at", ""), reverse=True)

    def get_active_model(self) -> Optional[dict]:
        rows = self.list_models()
        binaries = [r for r in rows if r.get("artifact_kind") == "model_binary"]
        active = [r for r in binaries if r.get("is_active")]
        return active[0] if active else (binaries[0] if binaries else None)

    def list_artifacts_for_run(self, run_id: str) -> list[dict]:
        run = self.load_run(run_id)
        if not run:
            return []
        return run.get("artifacts", [])