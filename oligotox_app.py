#!/usr/bin/env python3
"""
OligoTox Operator Console & Sequence Analysis Platform
Production-ready app wired to run_oligotox inference backend.
"""
import streamlit as st
import pandas as pd
import json
from pathlib import Path
from datetime import datetime
import sys

# PKPD imports
try:
    from src.pkpd.galnac_pkpd import simulate_galnac_pkpd
    PKPD_AVAILABLE = True
except ImportError:
    PKPD_AVAILABLE = False
    def simulate_galnac_pkpd(*args, **kwargs):
        return None

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / "src"))

from oligotox.config.paths import (
    PROCESSED_DATA_DIR,
    MODEL_ARTIFACTS_DIR,
    METRICS_ARTIFACTS_DIR,
)
from oligotox.inference.api import run_oligotox
# Simple sequence normalization inline


def normalize_sequence(seq: str) -> str:
    """Simple sequence normalization."""
    return seq.upper().replace(" ", "").replace("\n", "").replace("T", "U")
st.set_page_config(
    page_title="OligoTox Platform",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

st.title("🧬 OligoTox: Oligonucleotide Toxicity Prediction Platform")
st.caption("Mechanism-aware toxicity prediction and safety decision support for ASOs, siRNAs, and RNA-targeting therapeutics")

# Tabs
tabs = st.tabs(["Sequence Analysis", "Overview", "Curation", "Runs", "Models", "Predictions", "Evaluation"])

# ============================================================================
# TAB 0: SEQUENCE ANALYSIS (PRODUCTION INFERENCE)
# ============================================================================
with tabs[0]:
    st.subheader("Sequence Analysis")
    st.caption("Enter an oligonucleotide sequence for full toxicity risk profiling using the trained OligoTox model.")
    
    col1, col2 = st.columns([2, 1])
    
    with col1:
        seq_input = st.text_area(
            "Input Sequence (RNA/DNA)",
            height=120,
            placeholder="e.g., AUGCUAGCUAGCUAGCUAGC",
            help="Paste ASO, siRNA, or other oligo sequence (ACGTU alphabet)"
        )
        
    with col2:
        modality = st.selectbox("Modality", ["ASO", "siRNA", "lncRNA-targeting", "Other"], index=0)
        chemistry = st.selectbox("Chemistry", ["PS", "PO", "2'-OMe", "2'-MOE", "LNA", "Mixed"], index=0)
        conjugate = st.selectbox("Conjugate", ["None", "GalNAc", "Cholesterol", "Other"], index=0)
        
    run_btn = st.button("🔬 Run Full OligoTox Prediction", type="primary", use_container_width=True)
    
    if run_btn:
        if not seq_input or not seq_input.strip():
            st.error("Please enter a valid sequence.")
        else:
            try:
                seq_clean = normalize_sequence(seq_input)
                
                with st.spinner("Running OligoTox inference backend..."):
                    # Build input payload
                    inp = {
                        "sequence": seq_clean,
                        "modality": modality,
                        "chemistry": chemistry,
                        "conjugate": conjugate if conjugate != "None" else "",
                        "length": len(seq_clean),
                    }
                    
                    # Call production inference
                    result = run_oligotox(inp, omics_payload=None)
                    
                st.success("✅ Prediction complete!")
                
                # Display results
                st.markdown("### Toxicity Risk Profile")
                st.markdown(f"**Sequence:** `{seq_clean}` ({len(seq_clean)} nt)")
                st.markdown(f"**Modality:** {modality} | **Chemistry:** {chemistry} | **Conjugate:** {conjugate}")
                
                st.markdown("---")
                st.markdown("#### Per-Endpoint Risk Scores")
                
                if "per_endpoint" in result:
                    df_endpoints = pd.DataFrame(result["per_endpoint"]).T
                    df_endpoints = df_endpoints.reset_index().rename(columns={"index": "Endpoint"})
                    st.dataframe(df_endpoints, use_container_width=True)
                    
                    # Downloadable JSON
                    st.download_button(
                        label="📥 Download Full Result (JSON)",
                        data=json.dumps(result, indent=2),
                        file_name=f"oligotox_result_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json",
                        mime="application/json"
                    )
                
                if "dominant_endpoint" in result:
                    st.info(f"🎯 **Dominant Risk:** {result['dominant_endpoint']}")
                
                st.markdown("---")
                st.markdown("#### Mechanistic Attribution (D1/D2/D3)")
                st.markdown("""
                - **D1 (On-target):** Exaggerated pharmacology from intended target engagement  
                - **D2 (Hybridization off-target):** Unintended RNA binding via sequence complementarity  
                - **D3 (Chemistry/delivery):** Chemical-class or delivery-vehicle liabilities  
                """)
                
                # Display any mechanistic notes if available
                if "omics" in result and "adjustment_summary" in result["omics"]:
                    st.markdown(f"**Mechanistic Notes:** {result['omics']['adjustment_summary']}")
            
            except Exception as e:
                st.error(f"Prediction failed: {str(e)}")
                st.exception(e)

# ============================================================================
# TAB 1: OVERVIEW
# ============================================================================
with tabs[1]:
    st.subheader("Platform Overview")
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.metric("Training Rows", "1576")
    with col2:
        st.metric("Active Endpoints", "15+")
    with col3:
        st.metric("Model Status", "Trained")
    
    st.markdown("---")
    st.markdown("""
    ### OligoTox Platform
    
    **Mission:** Mechanism-aware toxicity prediction for oligonucleotide therapeutics.
    
    **Core Framework: D1/D2/D3**
    - **D1:** On-target pharmacology
    - **D2:** Hybridization-dependent off-target effects
    - **D3:** Chemistry and delivery liabilities
    
    **Key Features:**
    - Endpoint-specific risk scores
    - Uncertainty and confidence indicators
    - Mechanistic attribution
    - Assay and biomarker recommendations
    - Candidate comparison support
    """)

# ============================================================================
# TAB 2: CURATION
# ============================================================================
with tabs[2]:
    st.subheader("Data Curation Status")
    
    try:
        dataset_path = PROCESSED_DATA_DIR / "training_dataset.csv"
        if dataset_path.exists():
            df = pd.read_csv(dataset_path, nrows=10)
            st.success(f"✅ Training dataset found: {dataset_path}")
            st.dataframe(df.head())
            
            # Show sequence completeness
            if "sequence" in df.columns:
                full_df = pd.read_csv(dataset_path)
                n_total = len(full_df)
                n_with_seq = full_df["sequence"].notna().sum()
                st.metric("Sequence Coverage", f"{n_with_seq}/{n_total}")
        else:
            st.warning("Training dataset not found.")
    except Exception as e:
        st.error(f"Error loading dataset: {e}")

# ============================================================================
# TAB 3: RUNS
# ============================================================================
with tabs[3]:
    st.subheader("Pipeline Runs")
    st.info("This tab shows rebuild, training, and batch prediction run history.")
    st.markdown("*(Run tracking not yet implemented in this version)*")

# ============================================================================
# TAB 4: MODELS
# ============================================================================
with tabs[4]:
    st.subheader("Model Registry")
    
    try:
        if MODEL_ARTIFACTS_DIR.exists():
            models = list(MODEL_ARTIFACTS_DIR.glob("*.pkl")) + list(MODEL_ARTIFACTS_DIR.glob("*.joblib"))
            if models:
                st.success(f"Found {len(models)} model artifact(s)")
                for m in models:
                    st.markdown(f"- `{m.name}` ({m.stat().st_size // 1024} KB)")
            else:
                st.info("No model binaries found. Train a model first.")
        else:
            st.warning("Model artifacts directory does not exist.")
    except Exception as e:
        st.error(f"Error reading model registry: {e}")

# ============================================================================
# TAB 5: PREDICTIONS
# ============================================================================
with tabs[5]:
    st.subheader("Recent Predictions")
    st.info("This tab tracks recent prediction runs and batch analyses.")
    st.markdown("*(Prediction logging not yet implemented in this version)*")

# ============================================================================
# TAB 6: EVALUATION (NEW - CLINICAL PREDICTION METRICS)
# ============================================================================
with tabs[6]:
    st.subheader("📊 Model Evaluation & Clinical Metrics")
    st.caption("Rigorous evaluation with leakage controls, calibration, and uncertainty quantification")
    
    st.markdown("""
    ### Evaluation Philosophy
    
    This platform follows best practices for clinical prediction:
    
    1. **Leakage-Resistant Splitting:** DataSAIL-based train/val/test splits accounting for:
       - Sequence similarity
       - Batch effects
       - Temporal structure
       - Study/assay families
    
    2. **Clinically Meaningful Metrics:**
       - AUROC and AUPRC for discrimination
       - Brier score for calibration quality
       - Calibration curves (reliability diagrams)
       - Per-endpoint performance tracking
    
    3. **Uncertainty Quantification:**
       - Confidence bands (Low/Medium/High)
       - Coverage estimation
       - Out-of-distribution detection
    
    4. **Endpoint-Specific Evaluation:**
       - Risk band thresholds per endpoint
       - D1/D2/D3 mechanistic decomposition
       - Translational biomarker alignment
    """)
    
    st.markdown("---")
    
    # Load metrics if available
    try:
        metrics_file = METRICS_ARTIFACTS_DIR / "training_summary.json"
        if metrics_file.exists():
            with open(metrics_file, "r") as f:
                metrics = json.load(f)
            
            st.success("✅ Model evaluation metrics found")
            
            st.markdown("### Global Metrics")
            col1, col2, col3 = st.columns(3)
            
            # Extract aggregate metrics if available
            if "composite_auroc" in metrics:
                col1.metric("Composite AUROC", f"{metrics['composite_auroc']:.3f}")
            if "composite_brier" in metrics:
                col2.metric("Composite Brier Score", f"{metrics['composite_brier']:.3f}")
            if "n_endpoints" in metrics:
                col3.metric("Evaluated Endpoints", metrics["n_endpoints"])
            
            st.markdown("---")
            st.markdown("### Per-Endpoint Performance")
            
            if "per_endpoint" in metrics and isinstance(metrics["per_endpoint"], dict):
                endpoint_rows = []
                for ep, vals in metrics["per_endpoint"].items():
                    endpoint_rows.append({
                        "Endpoint": ep,
                        "AUROC": vals.get("auroc", "N/A"),
                        "AUPRC": vals.get("auprc", "N/A"),
                        "Brier": vals.get("brier", "N/A"),
                        "Threshold": vals.get("threshold", "N/A"),
                    })
                df_metrics = pd.DataFrame(endpoint_rows)
                st.dataframe(df_metrics, use_container_width=True)
            
            st.markdown("---")
            st.markdown("### Leakage Diagnostic")
            
            datasail_file = METRICS_ARTIFACTS_DIR / "datasail_split_diagnostic.json"
            if datasail_file.exists():
                with open(datasail_file, "r") as f:
                    datasail_info = json.load(f)
                st.json(datasail_info)
            else:
                st.info("DataSAIL leakage diagnostic not found. Run training pipeline to generate.")
            
        else:
            st.warning("⚠️ No evaluation metrics found. Train the model to generate metrics.")
            
            st.markdown("""\n            **To generate evaluation metrics:**\n            1. Ensure `data/processed/training_dataset.csv` exists\n            2. Run: `python train.py`\n            3. Metrics will be saved to `artifacts/metrics/training_summary.json`\n            """)
    
    except Exception as e:
        st.error(f"Error loading evaluation metrics: {e}")
        st.exception(e)
    
    st.markdown("---")
    st.markdown("### Calibration & Uncertainty")
    st.markdown("""
    **Calibration:** Does a predicted 70% risk correspond to 70% observed frequency?  
    **Uncertainty:** How confident is the model? Does it know when it doesn't know?
    
    These are critical for clinical decision support, where miscalibrated predictions can lead to poor therapeutic decisions.
    
    *(Calibration plots and uncertainty visualizations will be added in the next iteration)*
    """)

