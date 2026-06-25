"""
ML Toxicity Prediction Page  —  OligoTox v4
==========================================
Streamlit page that runs sequence-based RF models to predict
multi-endpoint toxicity risk for a user-supplied nucleotide sequence.

Navigate to this page from the main app sidebar.
"""
import sys, os
from pathlib import Path

# Ensure repo root on path
REPO = Path(__file__).resolve().parent.parent
for p in [str(REPO), str(REPO / "src")]:
    if p not in sys.path:
        sys.path.insert(0, p)

import pandas as pd
import streamlit as st

try:
    from oligotox.ml_predict import predict_sequence, get_top_features
    ML_AVAILABLE = True
except ImportError as e:
    ML_AVAILABLE = False
    ML_ERR = str(e)

st.set_page_config(page_title="ML Toxicity Prediction", page_icon="🧬", layout="wide")
st.title("🧬 ML Toxicity Prediction  —  OligoTox v4")

st.markdown("""
Enter a nucleotide sequence (ASO/siRNA, RNA or DNA) and receive per-endpoint
toxicity probability scores from the v4 Random Forest models trained on:
- **Papargyri 2020** LNA gapmer screening data (n=796 gapmers, hepatocellular injury + off-target)
- **MOE_Neurotox** in vivo ASO CNS safety data (n=2437 ASOs, FOB score ≥3 threshold)
- **Hagedorn 2022** in vitro Ca²⁺ oscillation neurotoxicity (n=1825 LNA ASOs)

CV-AUROC ranges: 0.70–0.91 across endpoints. Sequence-only prediction (no chemistry annotations required).
""")

if not ML_AVAILABLE:
    st.error(f"Could not import oligotox.ml_predict: {ML_ERR}")
    st.stop()

# ---- Input ----
col1, col2 = st.columns([3, 1])
with col1:
    seq_input = st.text_area(
        "Nucleotide sequence (5' to 3')",
        placeholder="e.g. TGGCATGATGCCTCGTAGACATGCT",
        height=100,
        help="RNA or DNA, upper or lower case. Spaces and dashes are stripped."
    )
with col2:
    st.markdown("### Example sequences")
    examples = {
        "Mipomersen (ASO)": "TCCGTCATCGCTCCTCAGGG",
        "Eteplirsen (PMO)": "CTCCAACATCAAGGAAGATGGCATTTCTAG",
        "Inclisiran (siRNA)": "AAGCUGCCGUGUCCUGGAAUU",
    }
    for name, seq in examples.items():
        if st.button(name, use_container_width=True):
            seq_input = seq
            st.session_state["seq_input"] = seq

if "seq_input" in st.session_state and not seq_input:
    seq_input = st.session_state["seq_input"]

run_btn = st.button("🔬 Run Prediction", type="primary", use_container_width=False)

if run_btn or seq_input:
    seq_clean = seq_input.strip().upper().replace(" ","").replace("-","")
    if len(seq_clean) < 5:
        st.warning("Please enter a sequence of at least 5 nucleotides.")
        st.stop()

    st.markdown(f"**Analysing:** `{seq_clean}` ({len(seq_clean)} nt)")

    with st.spinner("Running RF models..."):
        results = predict_sequence(seq_clean)

    if results and "error" in results[0]:
        st.error(results[0]["error"])
        st.stop()

    # ---- Results table ----
    df = pd.DataFrame(results)
    df["prob_%"] = (df["prob"] * 100).round(1)

    st.subheader("Per-Endpoint Risk Scores")

    def colour_risk(val):
        colours = {"High": "background-color: #ffcccc",
                   "Medium": "background-color: #fff3cd",
                   "Low": "background-color: #d4edda"}
        return colours.get(val, "")

    # Only select columns that exist in the results DataFrame
    base_cols = ["endpoint", "prob_%", "risk_level"]
    opt_cols = [c for c in ["n_train", "n_pos_train"] if c in df.columns]
    sel_cols = base_cols + opt_cols
    display = df[sel_cols].copy()
    col_map = {"endpoint": "Endpoint", "prob_%": "Probability (%)", "risk_level": "Risk Level", "n_train": "Training N", "n_pos_train": "Positive N"}
    display.columns = [col_map[c] for c in sel_cols]
    try:
        styled = display.style.map(colour_risk, subset=["Risk Level"])
    except AttributeError:
        styled = display.style.applymap(colour_risk, subset=["Risk Level"])
    st.dataframe(styled, hide_index=True)

    # ---- Bar chart ----
    st.subheader("Risk Score Overview")
    chart_df = df.set_index("endpoint")["prob_%"]
    st.bar_chart(chart_df, height=300)

    # ---- Feature importance drill-down ----
    st.subheader("🔍 Feature Importance Drill-Down")
    ep_choices = df["endpoint"].tolist()
    selected_ep = st.selectbox("Select endpoint to inspect:", ep_choices)
    fi = get_top_features(seq_clean, selected_ep, n=15)
    if fi:
        fi_df = pd.DataFrame(fi)
        fi_df.columns = ["Feature", "Importance", "Value"]
        st.dataframe(fi_df, hide_index=True)
        st.bar_chart(fi_df.set_index("Feature")["Importance"], height=300)
    else:
        st.info("Feature importances not available for this endpoint (model type does not support it).")

    # ---- Download ----
    csv = df.to_csv(index=False)
    st.download_button(
        label="⬇️ Download Results CSV",
        data=csv,
        file_name=f"oligotox_ml_prediction_{seq_clean[:12]}.csv",
        mime="text/csv"
    )

    # ---- Model performance note ----
    with st.expander("📊 Model Performance (CV-AUROC)"):
        perf = pd.read_csv(REPO / "reports" / "oligotox_v4_final_eval.csv")
        perf_display = perf[["endpoint","source","n_pos","model","auroc","ap"]].copy()
        perf_display.columns = ["Endpoint","Source","N pos","Model","AUROC","AP"]
        st.dataframe(perf_display, hide_index=True)

# Updated: v4 pipeline
