import streamlit as st
import sys
sys.path.insert(0, 'src')
from oligotox.inference.api import run_oligotox
import pandas as pd
import json

st.set_page_config(page_title="OligoTox Platform", layout="wide")

st.title("🧬 OligoTox Platform - One-of-a-Kind Safety & Efficacy Prediction")

st.markdown("""
### Multi-Endpoint Toxicity + Efficacy + PKPD Integration
**The only platform combining:**
- 13 Toxicity Endpoint Models ✅
- RNA Foundation Models (Nucleotide Transformer, RiNALMo, BiRNA-BERT) ✅
- 7 Efficacy Models (Coming: OligoFormer, ENsiRNA, siDPT, etc.) 🚧
- GalNAc PKPD Simulator 🚧
- OSWG/FDA-Aligned Reporting 🚧
""")

# Input section
st.header("🧬 1. Input Sequence")
col1, col2, col3 = st.columns(3)

with col1:
    sequence = st.text_area("RNA/DNA Sequence", "AUGCUAGCUAGCUAGCUAGC", height=100)

with col2:
    modality = st.selectbox("Modality", ["ASO", "siRNA", "mRNA", "PMO"])
    
with col3:
    chemistry = st.selectbox("Chemistry", ["PS", "MOE", "2'-OMe", "LNA", "PMO"])
    conjugate = st.selectbox("Conjugate", ["None", "GalNAC", "Other"])

if st.button("🔬 Run Full OligoTox Prediction", type="primary"):
    with st.spinner("Running toxicity prediction..."):
        try:
            # Run prediction
            result = run_oligotox(sequence, modality=modality.lower())
            
            st.success("✅ Prediction Complete!")
            
            # Display results
            st.header("📊 2. Toxicity Risk Profile")
            
            # Extract predictions
            predictions = result.get('predictions', {})
            
            # Create dataframe
            df_results = pd.DataFrame([
                {"Endpoint": k, "Probability": v, "Risk": "High" if v > 0.7 else "Medium" if v > 0.3 else "Low"}
                for k, v in predictions.items()
            ]).sort_values('Probability', ascending=False)
            
            st.dataframe(df_results, use_container_width=True)
            
            # Visualize top risks
            import plotly.express as px
            fig = px.bar(df_results.head(10), x='Probability', y='Endpoint', 
                        orientation='h', title="Top 10 Toxicity Risks",
                        color='Risk', color_discrete_map={'High':'red', 'Medium':'orange', 'Low':'green'})
            st.plotly_chart(fig, use_container_width=True)
            
            # Download button
            st.download_button(
                label="📥 Download Results (JSON)",
                data=json.dumps(result, indent=2),
                file_name=f"oligotox_results_{sequence[:10]}.json",
                mime="application/json"
            )
            
        except Exception as e:
            st.error(f"❌ Error: {str(e)}")
            st.exception(e)

st.markdown("---")
st.markdown("""
### 🚧 Under Development:
- **Efficacy Models**: OligoFormer, ENsiRNA, siDPT, OligoGraph, siRNADiscovery, ASOptimizer, PFRED
- **PKPD Simulator**: GalNAc tissue distribution, 18-compartment PBPK model
- **Benchmarking**: OligoGym dataset integration
- **Reporting**: OSWG/FDA-aligned safety reports
""")
