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

from oligotox.features.candidate import candidate_features as build_candidate_features
from oligotox.pipelines.score_candidate import score_candidate
from oligotox.pkpd.bridge import compute_pkpd_enriched_features
from oligotox.efficacy.potency import score_efficacy
from oligotox.reporting.oswg_report import OSWGSafetyReport

APP_TITLE = "OligoTox Operator Console"
APP_SUBTITLE = (
    "Research prototype. Most toxicity endpoints rest on 74 unique sequences and are "
    "shown as unsupported; CNS neurotoxicity and siRNA potency are trained on external "
    "data. Each number carries the evidence behind it."
)

# Endpoints are presented by how well the evaluation supports them, not by how
# confident the model sounds: see scripts/derive_evidence_tiers.py.
TIER_LABELS = {
    "validated": "Validated",
    "weak": "Weak signal",
    "source_confounded": "Source-confounded",
    "insufficient_data": "Insufficient data",
}

PROCESSED_DIR = ROOT / "data" / "processed"
MODEL_DIR = ROOT / "artifacts" / "models"
CONFIG_DIR = ROOT / "artifacts" / "config"

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

def _sequence_features(seq: str, chemistry: str, dose_mg_kg: float, modality: str = "siRNA"):
    ps_fraction = 0.9 if "PS" in chemistry.upper() else 0.0
    return build_candidate_features(
        seq,
        modality=modality,
        chemistry=chemistry,
        ps_fraction=ps_fraction,
        dose_mg_kg=dose_mg_kg,
    )

def _run_toxicity(candidate_features):
    return score_candidate(
        candidate_features=candidate_features,
        processed_dir=PROCESSED_DIR,
        model_dir=MODEL_DIR,
        config_dir=CONFIG_DIR,
    )

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

def _build_report(seq, chemistry, toxicity_result, pkpd_result):
    toxicity_predictions = {}
    for endpoint, info in toxicity_result.get("per_endpoint", {}).items():
        val = _extract_endpoint_value(info)
        if isinstance(val, (float, int)):
            toxicity_predictions[endpoint] = float(val)

    pkpd_summary = {}
    if isinstance(pkpd_result, dict):
        for k, v in pkpd_result.items():
            if isinstance(v, (float, int, str, bool)):
                pkpd_summary[k] = v

    evidence_tiers = {
        endpoint: info.get("evidence_tier")
        for endpoint, info in toxicity_result.get("per_endpoint", {}).items()
        if isinstance(info, dict)
    }

    report = OSWGSafetyReport(
        compound_id=f"RUN-{datetime.utcnow().strftime('%Y%m%d%H%M%S')}",
        sequence=seq,
        assessment_date=datetime.utcnow().strftime("%Y-%m-%d"),
        toxicity_predictions=toxicity_predictions,
        pkpd_summary=pkpd_summary,
        evidence_tiers=evidence_tiers,
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
    """Keep the run in this visitor's session only.

    Runs used to be written to artifacts/runs/ and the Saved runs tab listed
    every file it found there. On a shared deployment that is one visitor's
    sequences shown to the next, so runs now live in session state and leave
    with the session. Use the download buttons to keep one.
    """
    run_id = payload.get("run_id") or str(uuid.uuid4())
    payload["run_id"] = run_id
    payload["saved_at"] = datetime.utcnow().isoformat() + "Z"
    st.session_state.setdefault("runs", []).insert(0, payload)
    return run_id

def _load_saved_runs():
    """This session's runs, newest first."""
    rows = []
    for obj in st.session_state.get("runs", []):
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
            }
        )
    return pd.DataFrame(rows)


def _render_endpoint_table(toxicity_result, tiers=None):
    rows = []
    for endpoint, info in toxicity_result.get("per_endpoint", {}).items():
        if tiers is not None and info.get("evidence_tier") not in tiers:
            continue
        rows.append(
            {
                "endpoint": endpoint,
                "value": _extract_endpoint_value(info),
                # A number that is not calibrated is a rank, not a percentage.
                "reads as": ("probability" if info.get("probabilities_calibrated")
                             else "ranking only"),
                "threshold": info.get("threshold") if isinstance(info, dict) else None,
                "predicted_positive": info.get("predicted_positive") if isinstance(info, dict) else None,
                "risk_band": info.get("risk_band") if isinstance(info, dict) else None,
            }
        )
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    return df.sort_values("value", ascending=False, na_position="last")


def _render_unsupported_table(toxicity_result):
    """Endpoints the evaluation doesn't support, with why."""
    rows = []
    for endpoint, info in toxicity_result.get("per_endpoint", {}).items():
        tier = info.get("evidence_tier")
        if tier == "validated":
            continue
        rows.append(
            {
                "endpoint": endpoint,
                "tier": TIER_LABELS.get(tier, tier),
                "model output": _extract_endpoint_value(info),
                "held-out PR-AUC": info.get("cv_pr_auc"),
                "base rate": info.get("base_rate"),
                "positive sequences": info.get("positive_sequences"),
                "why it is not a prediction": info.get("evidence_note"),
            }
        )
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    order = {label: i for i, label in enumerate(TIER_LABELS.values())}
    return df.sort_values(
        ["tier", "positive sequences"],
        ascending=[True, False],
        key=lambda s: s.map(order) if s.name == "tier" else s,
    )

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
    run_btn = st.button("Run analysis", type="primary", width="stretch")

tabs = st.tabs([
    "Candidate Input",
    "Toxicity Results",
    "PK/PD",
    "Regulatory Report",
    "Efficacy",
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
        features_preview = _sequence_features(normalized_seq, chemistry, dose_mg_kg, modality)
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
    })

if run_btn:
    ok, seq_or_msg = _validate_sequence(sequence)
    if not ok:
        st.error(seq_or_msg)
    else:
        normalized_seq = seq_or_msg
        candidate_features = _sequence_features(normalized_seq, chemistry, dose_mg_kg, modality)

        with st.spinner("Running OligoTox analysis..."):
            toxicity_result = _run_toxicity(candidate_features)
            pkpd_result = _run_pkpd(normalized_seq, dose_mg_kg, modality, chemistry)
            report_result = _build_report(normalized_seq, chemistry, toxicity_result, pkpd_result)

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
            "pkpd": pkpd_result,
            "report": report_result,
        }
        run_id = _save_run(payload)
        st.session_state["last_run"] = payload
        st.success(f"Analysis complete. Run {run_id} kept in this session.")

run = st.session_state.get("last_run")

with tabs[1]:
    st.subheader("Toxicity")
    if not run:
        st.info("Run an analysis to view toxicity outputs.")
    else:
        tox = run["toxicity"]
        basis = tox.get("composite_basis") or []
        total = len(tox.get("per_endpoint", {}))

        if tox.get("composite_score") is None:
            st.warning(
                "No endpoint in this build passes the evaluation, so no composite score is "
                "reported. The per-endpoint outputs below are not predictions."
            )
        else:
            c1, c2, c3 = st.columns(3)
            c1.metric("Composite score", f"{_safe_float(tox.get('composite_score')):.3f}")
            c2.metric("Band", str(tox.get("composite_band")))
            c3.metric("Highest validated endpoint", str(tox.get("dominant_endpoint")))
            st.caption(
                f"Composite covers {len(basis)} of {total} endpoints "
                f"({', '.join(e.replace('_', ' ') for e in basis)}). The other "
                f"{total - len(basis)} are excluded because the baselines don't support them."
            )

        st.markdown("**Predictions**")
        validated_df = _render_endpoint_table(tox, tiers={"validated"})
        if validated_df.empty:
            st.info("No validated endpoints in this build.")
        else:
            st.dataframe(validated_df, width="stretch", hide_index=True)

        unsupported_df = _render_unsupported_table(tox)
        if not unsupported_df.empty:
            with st.expander(
                f"Not supported by the evaluation ({len(unsupported_df)} endpoints) — "
                "model output shown, but do not read these as predictions"
            ):
                st.caption(
                    "These models run, but each one either loses to a baseline that sees only "
                    "the source dataset or one crude sequence number, or rests on fewer than 10 "
                    "positive sequences. They are excluded from the composite score."
                )
                st.dataframe(unsupported_df, width="stretch", hide_index=True)

        st.download_button(
            "Download toxicity JSON",
            data=json.dumps(tox, indent=2, default=str),
            file_name="toxicity_result.json",
            mime="application/json",
        )

with tabs[2]:
    st.subheader("PK/PD")
    if not run:
        st.info("Run an analysis to view PK/PD outputs.")
    else:
        pkpd = run["pkpd"]
        if isinstance(pkpd, dict):
            st.json(pkpd)
        else:
            st.write(pkpd)

with tabs[3]:
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


with tabs[4]:
    st.subheader("Efficacy")
    if not run:
        st.info("Run an analysis to view efficacy outputs.")
    else:
        seq = run["input"]["sequence"]
        modality = run["input"]["modality"]
        eff = score_efficacy(seq, modality, MODEL_DIR, CONFIG_DIR)
        if not eff["available"]:
            st.info(eff["reason"])
        else:
            for name, pred in eff["predictions"].items():
                label = name.replace("_", " ")
                if pred["tier"] == "validated":
                    st.metric(label, f"{pred['score']:.1f}", help=pred["units"])
                else:
                    st.metric(f"{label} (weak evidence)", f"{pred['score']:.2f}", help=pred["units"])
                st.caption(f"{pred['readout']}. {pred['evidence']}")
                if pred["caveat"]:
                    st.caption(f"Caveat: {pred['caveat']}")
                st.divider()
            st.caption(
                "These are raw model readouts in the units of the training data, not "
                "probabilities, and they are not comparable across modalities."
            )

with tabs[5]:
    st.subheader("Saved runs")
    st.caption(
        "Runs from this browser session only. They are not written to the server "
        "and are gone when you close the tab, so nothing you enter here is visible "
        "to anyone else. Download a run to keep it."
    )
    runs_df = _load_saved_runs()
    if runs_df.empty:
        st.info("No runs in this session yet.")
    else:
        st.dataframe(runs_df, width="stretch", hide_index=True)
        session_runs = st.session_state.get("runs", [])
        labels = {f"{r['saved_at']}  {r['input']['sequence'][:24]}": r for r in session_runs}
        selected = st.selectbox("Open a run from this session", list(labels))
        if selected:
            chosen = labels[selected]
            st.download_button(
                "Download this run as JSON",
                data=json.dumps(chosen, indent=2, default=str),
                file_name=f"oligotox_run_{chosen['run_id']}.json",
                mime="application/json",
            )
            st.json(chosen)

with tabs[6]:
    st.subheader("Build info")
    st.write({
        "app_file": str(ROOT / "oligotox_app.py"),
        "processed_dir": str(PROCESSED_DIR),
        "model_dir": str(MODEL_DIR),
        "config_dir": str(CONFIG_DIR),
        "run_storage": "session only, not written to disk",
    })
    st.write("Validation status: backend validated end-to-end in current workspace.")
