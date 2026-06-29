from __future__ import annotations

from dataclasses import dataclass, asdict
from datetime import datetime
from typing import Any, Dict, List


@dataclass
class OligotoxRun:
    run_id: str
    timestamp: str
    sequence: str
    modality: str
    chemistry: str
    conjugate: str
    inputs: Dict[str, Any]
    result: Dict[str, Any]
    regulatory_md: str


def _now_iso() -> str:
    return datetime.now().isoformat(timespec="seconds")


def make_run_id(seq: str) -> str:
    return f"run-{_now_iso().replace(':','').replace('-','')}-{len(seq)}"


def log_run(
    st,
    sequence: str,
    modality: str,
    chemistry: str,
    conjugate: str,
    inputs: Dict[str, Any],
    result: Dict[str, Any],
    regulatory_md: str,
) -> None:
    runs: List[Dict[str, Any]] = st.session_state.get("runs", [])
    run = OligotoxRun(
        run_id=make_run_id(sequence),
        timestamp=_now_iso(),
        sequence=sequence,
        modality=modality,
        chemistry=chemistry,
        conjugate=conjugate,
        inputs=inputs,
        result=result,
        regulatory_md=regulatory_md,
    )
    runs.append(asdict(run))
    st.session_state["runs"] = runs


def get_runs(st) -> List[Dict[str, Any]]:
    return st.session_state.get("runs", [])


def get_last_run(st) -> Dict[str, Any] | None:
    runs = st.session_state.get("runs", [])
    return runs[-1] if runs else None
