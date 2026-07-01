from __future__ import annotations

from dataclasses import dataclass, asdict
from datetime import datetime
from typing import Any, Dict, List, Optional


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
    ts = _now_iso().replace(":", "").replace("-", "")
    return f"run-{ts}-{len(seq)}"


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
    st.session_state["last_run"] = asdict(run)


def get_runs(st) -> List[Dict[str, Any]]:
    return st.session_state.get("runs", [])


def get_last_run(st) -> Optional[Dict[str, Any]]:
    last = st.session_state.get("last_run")
    if last:
        return last
    runs = st.session_state.get("runs", [])
    return runs[-1] if runs else None


def clear_runs(st) -> None:
    st.session_state["runs"] = []
    st.session_state["last_run"] = None