from __future__ import annotations

import json
import os
import shutil
import subprocess
import uuid
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional

from oligotox.registry import ModelRecord, PipelineRun, Registry
from oligotox.validators import validate_mapping_file, validate_sequence_input

REPO_ROOT = Path(os.getenv("OLIGOTOX_REPO_ROOT", Path.cwd()))
ARTIFACTS_DIR = REPO_ROOT / "artifacts"
REGISTRY = Registry(ARTIFACTS_DIR)

CANONICAL_REBUILD_STEPS = [
    ("mapping_validate", ["python3", "scripts/validate_candidate_mapping_file.py"]),
    ("dataset_rebuild", ["python3", "scripts/rebuild_base_training_dataset.py"]),
    ("add_chem_features", ["python3", "scripts/add_chem_features.py"]),
    ("build_dataset", ["python3", "scripts/build_dataset.py"]),
    ("rebuild_recheck", ["python3", "scripts/rebuild_dataset_recheck.py"]),
    ("audit_sequence_coverage", ["python3", "scripts/audit_sequence_coverage.py"]),
    ("train_baseline", ["python3", "scripts/train_baseline.py"]),
]

MODEL_BINARY_EXTS = {".pkl", ".joblib", ".pt", ".pth", ".onnx", ".bin", ".sav"}
REPORT_EXTS = {".csv", ".json", ".txt", ".md", ".html", ".parquet"}


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def _run_id(prefix: str) -> str:
    return f"{prefix}_{datetime.now().strftime('%Y%m%d_%H%M%S')}_{uuid.uuid4().hex[:8]}"


def _env() -> dict:
    env = os.environ.copy()
    env["PYTHONPATH"] = f"{Path.home()}:{REPO_ROOT}:{REPO_ROOT / 'src'}"
    return env


def _write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, default=str))


def _safe_read_text(path: Path) -> str | None:
    if not path.exists():
        return None
    try:
        return path.read_text()
    except Exception:
        return None


def _execute_step(run_dir: Path, step_name: str, command: List[str]) -> dict:
    log_path = run_dir / f"{step_name}.log"
    with log_path.open("w") as log:
        proc = subprocess.run(
            command,
            cwd=REPO_ROOT,
            env=_env(),
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True,
        )
    return {"step": step_name, "command": command, "returncode": proc.returncode, "log_path": str(log_path)}


def run_mapping_validate(trigger: str = "manual") -> dict:
    run_id = _run_id("mapping_validate")
    run_dir = ARTIFACTS_DIR / "runs" / run_id
    run_dir.mkdir(parents=True, exist_ok=True)

    run = PipelineRun(
        run_id=run_id,
        job_name="mapping_validate",
        status="running",
        trigger=trigger,
        logs_path=str(run_dir / "mapping_validate.log"),
    )
    REGISTRY.save_run(run)

    validation = validate_mapping_file(REPO_ROOT / "data/metadata/candidate_to_oligo_map.csv")
    payload = {"validation": asdict(validation)}
    _write_json(run_dir / "metrics.json", payload)
    (run_dir / "mapping_validate.log").write_text(json.dumps(payload, indent=2))

    run.status = "succeeded" if validation.ok else "failed"
    run.finished_at = utc_now_iso()
    run.metrics = payload
    run.artifacts.append(REGISTRY.register_artifact(run_dir / "metrics.json", kind="validation_metrics"))
    run.artifacts.append(REGISTRY.register_artifact(run_dir / "mapping_validate.log", kind="validation_log"))

    mapping_file = REPO_ROOT / "data/metadata/candidate_to_oligo_map.csv"
    if mapping_file.exists():
        run.artifacts.append(REGISTRY.register_artifact(mapping_file, kind="mapping_snapshot"))

    REGISTRY.save_run(run)
    return REGISTRY.load_run(run_id)


def run_full_rebuild(trigger: str = "manual", parameters: Optional[Dict[str, Any]] = None) -> dict:
    run_id = _run_id("dataset_rebuild_full")
    run_dir = ARTIFACTS_DIR / "runs" / run_id
    run_dir.mkdir(parents=True, exist_ok=True)

    run = PipelineRun(
        run_id=run_id,
        job_name="dataset_rebuild_full",
        status="running",
        trigger=trigger,
        parameters=parameters or {},
        logs_path=str(run_dir / "master.log"),
    )
    REGISTRY.save_run(run)

    step_results = []
    overall_ok = True
    for step_name, command in CANONICAL_REBUILD_STEPS:
        result = _execute_step(run_dir, step_name, command)
        step_results.append(result)
        if result["returncode"] != 0:
            overall_ok = False
            run.error = f"Step failed: {step_name}"
            break

    outputs = {
        "step_results": step_results,
        "expected_outputs": {
            "training_dataset": str(REPO_ROOT / "data/processed/training_dataset.csv"),
            "with_chem_dataset": str(REPO_ROOT / "data/processed/training_dataset_with_chem.csv"),
            "audit_file": str(REPO_ROOT / "data/processed/training_dataset_rebuilt_audit.csv"),
        },
    }

    run.status = "succeeded" if overall_ok else "failed"
    run.finished_at = utc_now_iso()
    run.metrics = {"n_steps": len(step_results), "n_steps_succeeded": sum(r["returncode"] == 0 for r in step_results)}
    _write_json(run_dir / "outputs.json", outputs)
    _write_json(run_dir / "metrics.json", run.metrics)

    for key, value in outputs["expected_outputs"].items():
        p = Path(value)
        if p.exists():
            run.artifacts.append(REGISTRY.register_artifact(p, kind=key))

    combined_logs = []
    for r in step_results:
        txt = _safe_read_text(Path(r["log_path"]))
        if txt:
            combined_logs.append(f"\n===== {r['step']} =====\n{txt}")
            run.artifacts.append(REGISTRY.register_artifact(r["log_path"], kind="step_log"))
    if combined_logs:
        (run_dir / "master.log").write_text("\n".join(combined_logs))
        run.artifacts.append(REGISTRY.register_artifact(run_dir / "master.log", kind="master_log"))

    run.artifacts.append(REGISTRY.register_artifact(run_dir / "outputs.json", kind="outputs_manifest"))
    run.artifacts.append(REGISTRY.register_artifact(run_dir / "metrics.json", kind="run_metrics"))
    REGISTRY.save_run(run)
    return REGISTRY.load_run(run_id)


def run_train_model_v4(trigger: str = "manual") -> dict:
    run_id = _run_id("train_model_v4")
    run_dir = ARTIFACTS_DIR / "runs" / run_id
    model_dir = ARTIFACTS_DIR / "models" / run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    model_dir.mkdir(parents=True, exist_ok=True)

    run = PipelineRun(
        run_id=run_id,
        job_name="train_model_v4",
        status="running",
        trigger=trigger,
        logs_path=str(run_dir / "train_model_v4.log"),
    )
    REGISTRY.save_run(run)

    result = _execute_step(run_dir, "train_model_v4", ["python3", "scripts/train_baseline.py"])
    run.finished_at = utc_now_iso()
    run.artifacts.append(REGISTRY.register_artifact(result["log_path"], kind="training_log"))

    if result["returncode"] != 0:
        run.status = "failed"
        run.error = "Model training script failed."
        REGISTRY.save_run(run)
        return REGISTRY.load_run(run_id)

    discovered = []
    stable_model_drop = REPO_ROOT / "artifacts" / "models" / "latest_train_baseline"

    if stable_model_drop.exists():
        for p in stable_model_drop.rglob("*"):
            if p.is_file():
                discovered.append(p)

    for folder in [REPO_ROOT / "artifacts", REPO_ROOT / "reports", REPO_ROOT / "models"]:
        if folder.exists():
            for p in folder.rglob("*"):
                if p.is_file() and p not in discovered and any(tok in p.name.lower() for tok in ["model", "baseline", "logreg", "eval", "metric"]):
                    discovered.append(p)

    copied_binaries = []
    copied_reports = []
    copied_other = []
    seen = set()

    for p in discovered:
        if str(p) in seen:
            continue
        seen.add(str(p))
        dest = model_dir / p.name
        try:
            if p.resolve() != dest.resolve():
                shutil.copy2(p, dest)
        except Exception:
            try:
                if p.exists():
                    dest.write_bytes(p.read_bytes())
            except Exception:
                continue

        ext = dest.suffix.lower()
        if ext in MODEL_BINARY_EXTS:
            copied_binaries.append(dest)
            run.artifacts.append(REGISTRY.register_artifact(dest, kind="model_binary"))
        elif ext in REPORT_EXTS:
            copied_reports.append(dest)
            run.artifacts.append(REGISTRY.register_artifact(dest, kind="model_report"))
        else:
            copied_other.append(dest)
            run.artifacts.append(REGISTRY.register_artifact(dest, kind="model_related_file"))

    active_assigned = False
    for p in sorted(copied_binaries):
        REGISTRY.register_model(
            ModelRecord(
                model_name=f"{run_id}:{p.name}",
                endpoint=p.stem,
                artifact_path=str(p),
                metrics={"source_file": p.name},
                feature_schema_version="v4",
                training_dataset_hash=None,
                source_run_id=run_id,
                is_active=not active_assigned,
                artifact_kind="model_binary",
            )
        )
        active_assigned = True

    for p in sorted(copied_reports):
        REGISTRY.register_model(
            ModelRecord(
                model_name=f"{run_id}:{p.name}",
                endpoint="report",
                artifact_path=str(p),
                metrics={"source_file": p.name},
                feature_schema_version="v4",
                training_dataset_hash=None,
                source_run_id=run_id,
                is_active=False,
                artifact_kind="model_report",
            )
        )

    run.status = "succeeded"
    run.metrics = {
        "copied_model_binaries": len(copied_binaries),
        "copied_model_reports": len(copied_reports),
        "copied_other_files": len(copied_other),
        "active_model_assigned": active_assigned,
        "stable_model_drop": str(stable_model_drop),
    }
    _write_json(run_dir / "metrics.json", run.metrics)
    run.artifacts.append(REGISTRY.register_artifact(run_dir / "metrics.json", kind="training_metrics"))
    REGISTRY.save_run(run)
    return REGISTRY.load_run(run_id)


def run_predict_sequence(sequence: str, trigger: str = "manual") -> dict:
    from oligotox import ml_predict

    validation = validate_sequence_input(sequence)
    if not validation.ok:
        return {"ok": False, "errors": validation.errors}

    run_id = _run_id("predict_sequence")
    run_dir = ARTIFACTS_DIR / "runs" / run_id
    run_dir.mkdir(parents=True, exist_ok=True)

    run = PipelineRun(
        run_id=run_id,
        job_name="predict_sequence",
        status="running",
        trigger=trigger,
        parameters={"sequence": validation.normalized_value},
    )
    REGISTRY.save_run(run)

    if hasattr(ml_predict, "predict_sequence"):
        prediction = ml_predict.predict_sequence(validation.normalized_value)
    elif hasattr(ml_predict, "main_predict"):
        prediction = ml_predict.main_predict(validation.normalized_value)
    else:
        raise AttributeError("ml_predict.py must expose predict_sequence() or main_predict().")

    if hasattr(ml_predict, "get_top_features"):
        try:
            top_features = ml_predict.get_top_features(validation.normalized_value)
        except Exception as e:
            top_features = {"warning": str(e)}
    else:
        top_features = {"warning": "get_top_features not available"}

    active_model = REGISTRY.get_active_model()
    payload = {
        "input": validation.normalized_value,
        "validation": asdict(validation),
        "prediction": prediction,
        "top_features": top_features,
        "run_id": run_id,
        "active_model": active_model,
    }

    _write_json(run_dir / "prediction.json", payload)
    run.status = "succeeded"
    run.finished_at = utc_now_iso()
    run.metrics = {"sequence_length": validation.stats.get("length")}
    run.artifacts.append(REGISTRY.register_artifact(run_dir / "prediction.json", kind="prediction"))
    REGISTRY.save_run(run)
    return {"ok": True, "run_id": run_id, "payload": payload}
