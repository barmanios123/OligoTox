from __future__ import annotations

from pathlib import Path
import pandas as pd
import streamlit as st

from oligotox.jobs import REGISTRY, run_full_rebuild, run_mapping_validate, run_train_model_v4
from oligotox.system_status import get_system_status
from oligotox.ui_state import latest_rebuild_run, latest_training_run, list_recent_prediction_runs
from oligotox.validators import validate_mapping_file

st.set_page_config(page_title="OligoTox Operator", layout="wide")
st.title("OligoTox Operator Console")

status = get_system_status(".")
active_model = REGISTRY.get_active_model()

if status.get("trust_banner"):
    st.warning(status["trust_banner"])

tabs = st.tabs(["Overview", "Curation", "Runs", "Models", "Predictions"])

with tabs[0]:
    st.subheader("Overview")
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Total training rows", status.get("total_rows") or "NA")
    c2.metric("Rows with sequence", status.get("rows_with_sequence") or "NA")
    c3.metric("Rows missing sequence", status.get("rows_missing_sequence") or "NA")
    c4.metric("Active sequence subset", status.get("active_sequence_subset_rows") or "NA")

    st.info(
        "System status: pipeline is operational, but current model-building evidence remains constrained by limited "
        "sequence-linked coverage in the labeled dataset."
    )

    latest_rebuild = latest_rebuild_run()
    latest_train = latest_training_run()

    left, right = st.columns(2)
    with left:
        st.markdown("### Latest rebuild")
        if latest_rebuild:
            st.json(latest_rebuild)
        else:
            st.caption("No rebuild runs yet.")
    with right:
        st.markdown("### Active model binary")
        if active_model:
            st.json(active_model)
        else:
            st.caption("No active model binary registered yet.")

    mapping_path = Path("data/metadata/candidate_to_oligo_map.csv")
    if mapping_path.exists():
        result = validate_mapping_file(mapping_path)
        a, b, c = st.columns(3)
        a.metric("Mapping rows", result.stats.get("rows", 0))
        b.metric("Rows with mapping info", result.stats.get("rows_with_any_mapping_info", 0))
        c.metric("Rows unresolved", result.stats.get("rows_unresolved", 0))

with tabs[1]:
    st.subheader("Curation")
    mapping_path = Path("data/metadata/candidate_to_oligo_map.csv")
    if mapping_path.exists():
        df = pd.read_csv(mapping_path, low_memory=False)
        edited = st.data_editor(df, num_rows="dynamic", use_container_width=True)
        c1, c2 = st.columns(2)
        with c1:
            if st.button("Save mapping CSV"):
                edited.to_csv(mapping_path, index=False)
                st.success(f"Saved {mapping_path}")
        with c2:
            if st.button("Validate mapping file"):
                result = run_mapping_validate(trigger="ui")
                st.json(result)
    else:
        st.error("Mapping file not found.")

with tabs[2]:
    st.subheader("Runs")
    c1, c2 = st.columns(2)
    with c1:
        if st.button("Run full rebuild pipeline"):
            result = run_full_rebuild(trigger="ui")
            st.success(f"Rebuild complete: {result['run_id']}")
            st.json(result)
    with c2:
        if st.button("Train / register v4 model"):
            result = run_train_model_v4(trigger="ui")
            if result.get("status") == "succeeded":
                st.success(f"Training complete: {result['run_id']}")
            else:
                st.error(f"Training failed: {result.get('error')}")
            st.json(result)

    for run in REGISTRY.list_runs()[:30]:
        with st.expander(f"{run['run_id']} — {run['job_name']} — {run['status']}"):
            st.json(run)
            for artifact in run.get("artifacts", []):
                p = artifact.get("path")
                if p and Path(p).exists():
                    try:
                        with open(p, "rb") as f:
                            st.download_button(
                                label=f"Download {artifact.get('kind', 'artifact')}: {artifact.get('name', Path(p).name)}",
                                data=f.read(),
                                file_name=Path(p).name,
                                key=f"{run['run_id']}::{p}",
                            )
                    except Exception:
                        pass

with tabs[3]:
    st.subheader("Models")
    models = REGISTRY.list_models()
    if models:
        df = pd.DataFrame(models)
        binaries = df[df["artifact_kind"] == "model_binary"] if "artifact_kind" in df.columns else df.iloc[0:0]
        reports = df[df["artifact_kind"] == "model_report"] if "artifact_kind" in df.columns else df.iloc[0:0]

        st.markdown("### Model binaries")
        if not binaries.empty:
            st.dataframe(binaries, use_container_width=True)
            if active_model and active_model.get("artifact_path") and Path(active_model["artifact_path"]).exists():
                with open(active_model["artifact_path"], "rb") as f:
                    st.download_button(
                        "Download active model binary",
                        data=f.read(),
                        file_name=Path(active_model["artifact_path"]).name,
                        key="active_model_download",
                    )
        else:
            st.info("No model binaries registered yet.")

        st.markdown("### Reports and metrics")
        if not reports.empty:
            st.dataframe(reports, use_container_width=True)
            for _, row in reports.iterrows():
                p = Path(row["artifact_path"])
                if p.exists():
                    with open(p, "rb") as f:
                        st.download_button(
                            f"Download report: {p.name}",
                            data=f.read(),
                            file_name=p.name,
                            key=f"report::{p}",
                        )
        else:
            st.info("No model reports registered yet.")
    else:
        st.info("No models registered yet.")

with tabs[4]:
    st.subheader("Recent predictions")
    preds = list_recent_prediction_runs(limit=20)
    if preds:
        st.dataframe(pd.DataFrame(preds), use_container_width=True)
    else:
        st.info("No prediction runs yet. Use the ML prediction page to generate one.")
