from __future__ import annotations
import json
from pathlib import Path
from typing import List, Dict, Any
from oligotox.jobs import REGISTRY

def list_recent_prediction_runs(limit: int = 20) -> List[Dict[str, Any]]:
    runs = REGISTRY.list_runs()
    preds = [r for r in runs if r.get("job_name") == "predict_sequence"]
    return preds[:limit]

def load_prediction_payload(run_id: str):
    pred_path = Path("artifacts") / "runs" / run_id / "prediction.json"
    if not pred_path.exists():
        return None
    try:
        return json.loads(pred_path.read_text())
    except Exception:
        return None

def latest_rebuild_run():
    runs = REGISTRY.list_runs()
    rebuilds = [r for r in runs if r.get("job_name") == "dataset_rebuild_full"]
    return rebuilds[0] if rebuilds else None

def latest_training_run():
    runs = REGISTRY.list_runs()
    trainings = [r for r in runs if r.get("job_name") == "train_model_v4"]
    return trainings[0] if trainings else None
