from __future__ import annotations

from fastapi import FastAPI, HTTPException
from pydantic import BaseModel

from oligotox.jobs import REGISTRY, run_full_rebuild, run_mapping_validate, run_predict_sequence

app = FastAPI(title="OligoTox API", version="0.1.0")


class PredictRequest(BaseModel):
    sequence: str


@app.get("/health")
def health():
    return {"ok": True, "service": "oligotox-api"}


@app.get("/models")
def list_models():
    return {"models": REGISTRY.list_models()}


@app.get("/runs")
def list_runs():
    return {"runs": REGISTRY.list_runs()}


@app.get("/runs/{run_id}")
def get_run(run_id: str):
    run = REGISTRY.load_run(run_id)
    if not run:
        raise HTTPException(status_code=404, detail="Run not found")
    return run


@app.post("/jobs/mapping-validate")
def api_mapping_validate():
    return run_mapping_validate(trigger="api")


@app.post("/jobs/rebuild")
def api_rebuild():
    return run_full_rebuild(trigger="api")


@app.post("/predict")
def api_predict(body: PredictRequest):
    result = run_predict_sequence(body.sequence, trigger="api")
    if not result.get("ok"):
        raise HTTPException(status_code=400, detail=result.get("errors"))
    return result
