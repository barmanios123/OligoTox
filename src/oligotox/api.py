from __future__ import annotations

from datetime import datetime, timezone
from typing import Any, Dict, List, Literal, Optional

from fastapi import FastAPI, HTTPException, Query
from pydantic import BaseModel, Field, field_validator

from oligotox.jobs import REGISTRY, run_full_rebuild, run_mapping_validate, run_predict_sequence

APP_NAME = "OligoTox API"
APP_VERSION = "0.2.0"


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def normalize_sequence(seq: str) -> str:
    return seq.upper().replace(" ", "").replace("\n", "").replace("\t", "").replace("T", "U")


def validate_sequence(seq: str) -> list[str]:
    errors: list[str] = []
    if not seq:
        return ["Sequence is required."]
    allowed = set("ACGUN")
    bad = sorted(set(ch for ch in seq if ch not in allowed))
    if bad:
        errors.append(
            f"Invalid characters found: {', '.join(bad)}. Allowed characters are A, C, G, U, N, and DNA T (auto-normalized to U)."
        )
    if len(seq) < 8:
        errors.append("Sequence is too short; minimum supported length is 8 nt.")
    if len(seq) > 80:
        errors.append("Sequence is too long; maximum supported length is 80 nt for this API route.")
    return errors


def registry_models_safe() -> list[Any]:
    try:
        models = REGISTRY.list_models()
        return models if isinstance(models, list) else list(models or [])
    except Exception:
        return []


def registry_runs_safe() -> list[Any]:
    try:
        runs = REGISTRY.list_runs()
        return runs if isinstance(runs, list) else list(runs or [])
    except Exception:
        return []


app = FastAPI(
    title=APP_NAME,
    version=APP_VERSION,
    description=(
        "API surface for OligoTox sequence prediction, run tracking, and operational rebuild jobs. "
        "This service exposes the currently wired job backend through a cleaner production-style interface."
    ),
)


class HealthResponse(BaseModel):
    ok: bool
    service: str
    version: str
    timestamp_utc: str
    model_count: int
    run_count: int


class PredictRequest(BaseModel):
    sequence: str = Field(..., description="Oligonucleotide sequence. DNA T is auto-normalized to U.")
    modality: Optional[str] = Field(default=None, description="Optional modality context, e.g. ASO or siRNA.")
    chemistry: Optional[str] = Field(default=None, description="Optional chemistry context, e.g. PS, LNA, 2'-MOE.")
    conjugate: Optional[str] = Field(default=None, description="Optional conjugate context, e.g. GalNAc.")
    trigger: str = Field(default="api", description="Trigger source for run tracking.")

    @field_validator("sequence")
    @classmethod
    def check_sequence(cls, value: str) -> str:
        seq = normalize_sequence(value)
        errors = validate_sequence(seq)
        if errors:
            raise ValueError(" ".join(errors))
        return seq


class PredictResponse(BaseModel):
    ok: bool
    input_sequence: str
    normalized_sequence: str
    sequence_length: int
    trigger: str
    timestamp_utc: str
    result: Dict[str, Any]


class ModelsResponse(BaseModel):
    models: List[Any]
    count: int


class RunsResponse(BaseModel):
    runs: List[Any]
    count: int


class JobResponse(BaseModel):
    ok: bool = True
    job: str
    trigger: str
    timestamp_utc: str
    payload: Dict[str, Any]


class ErrorResponse(BaseModel):
    ok: bool = False
    error: str
    detail: Any | None = None
    timestamp_utc: str


@app.get("/", tags=["meta"])
def root() -> Dict[str, Any]:
    return {
        "ok": True,
        "service": APP_NAME,
        "version": APP_VERSION,
        "docs": "/docs",
        "redoc": "/redoc",
        "timestamp_utc": utc_now_iso(),
    }


@app.get("/health", response_model=HealthResponse, tags=["meta"])
def health() -> HealthResponse:
    models = registry_models_safe()
    runs = registry_runs_safe()
    return HealthResponse(
        ok=True,
        service="oligotox-api",
        version=APP_VERSION,
        timestamp_utc=utc_now_iso(),
        model_count=len(models),
        run_count=len(runs),
    )


@app.get("/models", response_model=ModelsResponse, tags=["registry"])
def list_models() -> ModelsResponse:
    models = registry_models_safe()
    return ModelsResponse(models=models, count=len(models))


@app.get("/runs", response_model=RunsResponse, tags=["registry"])
def list_runs(
    limit: int = Query(default=100, ge=1, le=1000, description="Maximum number of runs to return."),
) -> RunsResponse:
    runs = registry_runs_safe()
    return RunsResponse(runs=runs[:limit], count=min(len(runs), limit))


@app.get("/runs/{run_id}", tags=["registry"])
def get_run(run_id: str) -> Dict[str, Any]:
    try:
        run = REGISTRY.load_run(run_id)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to load run: {e}") from e

    if not run:
        raise HTTPException(status_code=404, detail="Run not found")
    return run


@app.post("/jobs/mapping-validate", response_model=JobResponse, tags=["jobs"])
def api_mapping_validate(
    trigger: Literal["api", "manual", "system"] = Query(default="api"),
) -> JobResponse:
    try:
        payload = run_mapping_validate(trigger=trigger)
        return JobResponse(
            ok=True,
            job="mapping-validate",
            trigger=trigger,
            timestamp_utc=utc_now_iso(),
            payload=payload if isinstance(payload, dict) else {"result": payload},
        )
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Mapping validation failed: {e}") from e


@app.post("/jobs/rebuild", response_model=JobResponse, tags=["jobs"])
def api_rebuild(
    trigger: Literal["api", "manual", "system"] = Query(default="api"),
) -> JobResponse:
    try:
        payload = run_full_rebuild(trigger=trigger)
        return JobResponse(
            ok=True,
            job="rebuild",
            trigger=trigger,
            timestamp_utc=utc_now_iso(),
            payload=payload if isinstance(payload, dict) else {"result": payload},
        )
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Full rebuild failed: {e}") from e


@app.post("/predict", response_model=PredictResponse, tags=["prediction"])
def api_predict(body: PredictRequest) -> PredictResponse:
    seq = normalize_sequence(body.sequence)

    try:
        raw_result = run_predict_sequence(seq, trigger=body.trigger)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Prediction execution failed: {e}") from e

    if not isinstance(raw_result, dict):
        raise HTTPException(
            status_code=500,
            detail="Prediction backend returned a non-dictionary response.",
        )

    if not raw_result.get("ok"):
        raise HTTPException(
            status_code=400,
            detail=raw_result.get("errors", "Prediction failed"),
        )

    enriched_result = {
        **raw_result,
        "request_context": {
            "modality": body.modality,
            "chemistry": body.chemistry,
            "conjugate": body.conjugate,
        },
    }

    return PredictResponse(
        ok=True,
        input_sequence=body.sequence,
        normalized_sequence=seq,
        sequence_length=len(seq),
        trigger=body.trigger,
        timestamp_utc=utc_now_iso(),
        result=enriched_result,
    )


@app.get("/predict/{sequence}", response_model=PredictResponse, tags=["prediction"])
def api_predict_get(
    sequence: str,
    trigger: Literal["api", "manual", "system"] = Query(default="api"),
) -> PredictResponse:
    seq = normalize_sequence(sequence)
    errors = validate_sequence(seq)
    if errors:
        raise HTTPException(status_code=400, detail=errors)

    try:
        raw_result = run_predict_sequence(seq, trigger=trigger)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Prediction execution failed: {e}") from e

    if not isinstance(raw_result, dict):
        raise HTTPException(
            status_code=500,
            detail="Prediction backend returned a non-dictionary response.",
        )

    if not raw_result.get("ok"):
        raise HTTPException(
            status_code=400,
            detail=raw_result.get("errors", "Prediction failed"),
        )

    return PredictResponse(
        ok=True,
        input_sequence=sequence,
        normalized_sequence=seq,
        sequence_length=len(seq),
        trigger=trigger,
        timestamp_utc=utc_now_iso(),
        result=raw_result,
    )


@app.get("/system/status", tags=["meta"])
def system_status() -> Dict[str, Any]:
    return {
        "ok": True,
        "service": APP_NAME,
        "version": APP_VERSION,
        "timestamp_utc": utc_now_iso(),
        "registry": {
            "model_count": len(registry_models_safe()),
            "run_count": len(registry_runs_safe()),
        },
        "capabilities": {
            "prediction": True,
            "mapping_validate": True,
            "full_rebuild": True,
            "run_registry": True,
        },
    }