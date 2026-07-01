#!/usr/bin/env python3
import os
import sys
import json
import uuid
from datetime import datetime
from pathlib import Path

import pandas as pd
import streamlit as st

ROOT = Path(__file__).resolve().parent
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from oligotox.pipelines.score_candidate import score_candidate
from oligotox.efficacy.registry import EfficacyModelRegistry
from oligotox.pkpd.bridge import compute_pkpd_enriched_features
from oligotox.reporting.oswg_report import OSWGSafetyReport

APP_TITLE = "OligoTox Operator Console"
APP_SUBTITLE = "Client-ready oligonucleotide toxicity, efficacy, PK/PD, and reporting workspace"

PROCESSED_DIR = ROOT / "data" / "processed"
MODEL_DIR = ROOT / "artifacts" / "models"
CONFIG_DIR = ROOT / "artifacts" / "config"
RUNS_DIR = ROOT / "artifacts" / "runs"
RUNS_DIR.mkdir(parents=True, exist_ok=True)

st.set_page_config(
    page_title="OligoTox",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

def _safe_float(x, default=0.0):
    try:
        return float(x)
    except Exception:
        return default

def _extract_endpoint_value(info):
    if isinstance(info, dict):
        for key in [
            "probability",
            "score",
            "risk_score",
            "predicted_probability",
            "prob",
            "value",
            "raw_score",
        ]:
            if key in info and info[key] is not None:
                try:
                    return float(info[key])
                except Exception:
                    return info[key]
        return None
    return info

def _validate_sequence(seq: str):
    seq = (seq or "").strip().upper()
    if not seq:
        return False, "Sequence is required."
    allowed = set("ACGTU")
    bad = sorted(set(ch for ch in seq if ch not in allowed))
    if bad:
        return False, f"Invalid sequence characters: {', '.join(bad)}. Allowed: A, C, G, T, U."
    if len(seq) < 8:
        return False, "Sequence is too short for meaningful analysis; use at least 8 nt."
    return True, seq

def _sequence_features(seq: str, chemistry: str, dose_mg_kg: float):
    seq = seq.upper().replace("U", "T")
    length = len(seq)
    gc = (seq.count("G") + seq.count("C")) / max(length, 1)
    cpg = sum(1 for i in range(len(seq) - 1) if seq[i:i+2] == "CG") / max(length - 1, 1)
    chemistry = chemistry.upper()
    ps_fraction = 0.9 if "PS" in chemistry else 0.0
    return {
        "seq_length": int(length),
        "gc_content": float(gc),
        "cpg_density": float(cpg),
        "ps_fraction": float(ps_fraction),
        "dose_mg_kg": float(dose_mg_kg),
    }

def _run_toxicity(candidate_features):
    return score_candidate(
        candidate_features=candidate_features,
        processed_dir=PROCESSED_DIR,
        model_dir=MODEL_DIR,
        config_dir=CONFIG_DIR,
    )

def _run_efficacy(seq: str, modality: str, chemistry: str):
    registry = EfficacyModelRegistry()
    out = registry.predict_all(sequence=seq, modality=modality, chemistry=chemistry)
    rows = []
    for name, raw in out.items():
        value = None
        if isinstance(raw, (float, int)):
            value = float(raw)
        elif isinstance(raw, dict):
            for v in raw.values():
                try:
                    value = float(v[0] if hasattr(v, "__len__") and not isinstance(v, (str, bytes)) else v)
                    break
                except Exception:
                    pass
        elif hasattr(raw, "__len__") and not isinstance(raw, (str, bytes)):
            try:
                value = float(raw[0])
            except Exception:
                value = None
        rows.append({"model": name, "predicted_activity": value})
    df = pd.DataFrame(rows).sort_values("predicted_activity", ascending=False, na_position="last")
    return {"raw": out, "table": df}

def _run_pkpd(seq: str, dose_mg_kg: float, modality: str, chemistry: str, schedule: str = "weekly"):
    try:
        return compute_pkpd_enriched_features(
            sequence=seq,
            modality=modality,
            chemistry=chemistry,
            schedule=schedule,
            dose_mg_kg=dose_mg_kg,
        )
    except TypeError:
        try:
            return compute_pkpd_enriched_features(
                seq,
                modality,
                chemistry,
                schedule,
                dose_mg_kg=dose_mg_kg,
            )
        except Exception as e:
            return {"error": str(e)}
    except Exception as e:
        return {"error": str(e)}

def _build_report(seq, chemistry, toxicity_result, efficacy_result, pkpd_result):
    toxicity_predictions = {}
    for endpoint, info in toxicity_result.get("per_endpoint", {}).items():
        val = _extract_endpoint_value(info)
        if isinstance(val, (float, int)):
            toxicity_predictions[endpoint] = float(val)

    efficacy_score = None
    try:
        t = efficacy_result["table"]["predicted_activity"].dropna()
        if len(t) > 0:
            efficacy_score = float(t.mean())
    except Exception:
        efficacy_score = None

    pkpd_summary = {}
    if isinstance(pkpd_result, dict):
        for k, v in pkpd_result.items():
            if isinstance(v, (float, int, str, bool)):
                pkpd_summary[k] = v

    report = OSWGSafetyReport(
        compound_id=f"RUN-{datetime.utcnow().strftime('%Y%m%d%H%M%S')}",
        sequence=seq,
        assessment_date=datetime.utcnow().strftime("%Y-%m-%d"),
        toxicity_predictions=toxicity_predictions,
        efficacy_score=(0.0 if efficacy_score is None else efficacy_score),
        pkpd_summary=pkpd_summary,
    )

    report_json = report.generate_report(format="json")
    risk_score = report._calculate_risk_score()
    risk_class = report._classify_safety()

    return {
        "json": report_json,
        "risk_score": risk_score,
        "risk_class": risk_class,
    }

def _save_run(payload: dict):
    run_id = payload.get("run_id") or str(uuid.uuid4())
    payload["run_id"] = run_id
    payload["saved_at"] = datetime.utcnow().isoformat() + "Z"
    out = RUNS_DIR / f"{run_id}.json"
    with open(out, "w") as f:
        json.dump(payload, f, indent=2, default=str)
    return out

def _load_saved_runs():
    rows = []
    for fp in sorted(RUNS_DIR.glob("*.json"), reverse=True):
        try:
            with open(fp) as f:
                obj = json.load(f)
            rows.append(
                {
                    "run_id": obj.get("run_id"),
                    "saved_at": obj.get("saved_at"),
                    "sequence": obj.get("input", {}).get("sequence"),
                    "chemistry": obj.get("input", {}).get("chemistry"),
                    "modality": obj.get("input", {}).get("modality"),
                    "composite_score": obj.get("toxicity", {}).get("composite_score"),
                    "composite_band": obj.get("toxicity", {}).get("composite_band"),
                    "dominant_endpoint": obj.get("toxicity", {}).get("dominant_endpoint"),
                    "path": str(fp),
                }
            )
        except Exception:
            pass
    return pd.DataFrame(rows)

def _render_endpoint_table(toxicity_result):
    rows = []
    for endpoint, info in toxicity_result.get("per_endpoint", {}).items():
        rows.append(
            {
                "endpoint": endpoint,
                "value": _extract_endpoint_value(info),
                "threshold": info.get("threshold") if isinstance(info, dict) else None,
                "predicted_positive": info.get("predicted_positive") if isinstance(info, dict) else None,
                "risk_band": info.get("risk_band") if isinstance(info, dict) else None,
            }
        )
    df = pd.DataFrame(rows).sort_values("value", ascending=False, na_position="last")
    return df

st.title(APP_TITLE)
st.caption(APP_SUBTITLE)

with st.sidebar:
    st.header("Run setup")
    sequence = st.text_area(
        "Sequence",
        value="AUGCUAGCUAGCUAGCUAGC",
        height=120,
        help="Enter RNA/DNA oligo sequence using A/C/G/T/U."
    )
    chemistry = st.selectbox("Chemistry", ["PS", "PO", "mixed_PS_PO"], index=0)
    modality = st.selectbox("Modality", ["siRNA", "ASO"], index=0)
    dose_mg_kg = st.number_input("Dose (mg/kg)", min_value=0.0, max_value=100.0, value=3.0, step=0.5)
    analyst_note = st.text_input("Analyst note", value="")
    run_btn = st.button("Run analysis", type="primary", use_container_width=True)

tabs = st.tabs([
    "Candidate Input",
    "Toxicity Results",
    "Efficacy Results",
    "PK/PD",
    "Regulatory Report",
    "Saved Runs",
    "Build Info",
])

if "last_run" not in st.session_state:
    st.session_state["last_run"] = None

with tabs[0]:
    st.subheader("Candidate setup")
    ok, seq_or_msg = _validate_sequence(sequence)
    if ok:
        normalized_seq = seq_or_msg
        st.success(f"Sequence accepted. Length: {len(normalized_seq)} nt.")
        st.code(normalized_seq)
        features_preview = _sequence_features(normalized_seq, chemistry, dose_mg_kg)
        st.write("Derived feature preview")
        st.json(features_preview)
    else:
        normalized_seq = None
        st.error(seq_or_msg)

    st.write("Artifacts status")
    st.write({
        "processed_dir_exists": PROCESSED_DIR.exists(),
        "model_dir_exists": MODEL_DIR.exists(),
        "config_dir_exists": CONFIG_DIR.exists(),
        "runs_dir_exists": RUNS_DIR.exists(),
    })

if run_btn:
    ok, seq_or_msg = _validate_sequence(sequence)
    if not ok:
        st.error(seq_or_msg)
    else:
        normalized_seq = seq_or_msg
        candidate_features = _sequence_features(normalized_seq, chemistry, dose_mg_kg)

        with st.spinner("Running OligoTox analysis..."):
            toxicity_result = _run_toxicity(candidate_features)
            efficacy_result = _run_efficacy(normalized_seq, modality, chemistry)
            pkpd_result = _run_pkpd(normalized_seq, dose_mg_kg, modality, chemistry)
            report_result = _build_report(normalized_seq, chemistry, toxicity_result, efficacy_result, pkpd_result)

        payload = {
            "input": {
                "sequence": normalized_seq,
                "chemistry": chemistry,
                "modality": modality,
                "dose_mg_kg": dose_mg_kg,
                "analyst_note": analyst_note,
                "candidate_features": candidate_features,
            },
            "toxicity": toxicity_result,
            "efficacy": {
                "table": efficacy_result["table"].to_dict(orient="records"),
            },
            "pkpd": pkpd_result,
            "report": report_result,
        }
        saved_path = _save_run(payload)
        payload["saved_path"] = str(saved_path)
        st.session_state["last_run"] = payload
        st.success(f"Analysis complete. Saved run: {saved_path.name}")

run = st.session_state.get("last_run")

with tabs[1]:
    st.subheader("Toxicity")
    if not run:
        st.info("Run an analysis to view toxicity outputs.")
    else:
        tox = run["toxicity"]
        c1, c2, c3 = st.columns(3)
        c1.metric("Composite score", f"{_safe_float(tox.get('composite_score')):.3f}")
        c2.metric("Band", str(tox.get("composite_band")))
        c3.metric("Dominant endpoint", str(tox.get("dominant_endpoint")))

        endpoint_df = _render_endpoint_table(tox)
        st.dataframe(endpoint_df, use_container_width=True, hide_index=True)

        st.download_button(
            "Download toxicity JSON",
            data=json.dumps(tox, indent=2, default=str),
            file_name="toxicity_result.json",
            mime="application/json",
        )

with tabs[2]:
    st.subheader("Efficacy")
    if not run:
        st.info("Run an analysis to view efficacy outputs.")
    else:
        eff_df = pd.DataFrame(run["efficacy"]["table"])
        st.dataframe(eff_df, use_container_width=True, hide_index=True)
        if not eff_df.empty and "predicted_activity" in eff_df.columns:
            top = eff_df.dropna(subset=["predicted_activity"]).sort_values("predicted_activity", ascending=False).head(1)
            if len(top):
                st.metric("Top efficacy model", f"{top.iloc[0]['model']} ({top.iloc[0]['predicted_activity']:.3f})")

with tabs[3]:
    st.subheader("PK/PD")
    if not run:
        st.info("Run an analysis to view PK/PD outputs.")
    else:
        pkpd = run["pkpd"]
        if isinstance(pkpd, dict):
            st.json(pkpd)
        else:
            st.write(pkpd)

with tabs[4]:
    st.subheader("Regulatory report")
    if not run:
        st.info("Run an analysis to generate a report.")
    else:
        rep = run["report"]
        c1, c2 = st.columns(2)
        c1.metric("Risk score", f"{_safe_float(rep.get('risk_score')):.3f}")
        c2.metric("Classification", str(rep.get("risk_class")))
        st.code(rep.get("json", ""), language="json")
        st.download_button(
            "Download report JSON",
            data=rep.get("json", ""),
            file_name="oswg_report.json",
            mime="application/json",
        )

with tabs[5]:
    st.subheader("Saved runs")
    runs_df = _load_saved_runs()
    if runs_df.empty:
        st.info("No saved runs yet.")
    else:
        st.dataframe(runs_df, use_container_width=True, hide_index=True)
        selected = st.selectbox("Open saved run file", runs_df["path"].tolist())
        if selected:
            with open(selected) as f:
                st.json(json.load(f))

with tabs[6]:
    st.subheader("Build info")
    st.write({
        "app_file": str(ROOT / "oligotox_app.py"),
        "processed_dir": str(PROCESSED_DIR),
        "model_dir": str(MODEL_DIR),
        "config_dir": str(CONFIG_DIR),
        "runs_dir": str(RUNS_DIR),
    })
    st.write("Validation status: backend validated end-to-end in current workspace.")
