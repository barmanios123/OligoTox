import re

# Read the current oligotox_app.py
with open('oligotox_app.py', 'r') as f:
    content = f.read()

# 1. Add PK/PD import after other oligotox imports
import_section = '''from oligotox.inference.api import run_oligotox
# Simple sequence normalization inline'''

new_import_section = '''from oligotox.inference.api import run_oligotox
from oligotox.pkpd import adjust_toxicity_risk_with_pkpd  # Phase 2A PK/PD bridge
# Simple sequence normalization inline'''

content = content.replace(import_section, new_import_section)

# 2. Add PK/PD toggle and inputs in col2 (after conjugate selectbox, before run button)
conjugate_section = '''        conjugate = st.selectbox("Conjugate", ["None", "GalNAc", "Cholesterol", "Other"],
                index=0)
    
    run_btn = st.button("🔬 Run Full OligoTox Prediction", type="primary",
                        use_container_width=True)'''

new_conjugate_section = '''        conjugate = st.selectbox("Conjugate", ["None", "GalNAc", "Cholesterol", "Other"],
                index=0)
        
        # PK/PD Integration (Phase 2A)
        st.markdown("---")
        enable_pkpd = st.checkbox(
            "🧬 Enable PK/PD-Aware Scoring",
            value=False,
            help="Adjust toxicity predictions using GalNAc PK/PD modeling (Phase 2A)"
        )
        
        dose_mg_kg = 3.0
        schedule = "Q4W"
        if enable_pkpd:
            if conjugate == "GalNAc":
                dose_mg_kg = st.number_input(
                    "Dose (mg/kg)",
                    min_value=0.1,
                    max_value=10.0,
                    value=3.0,
                    step=0.5,
                    help="GalNAc-siRNA dose for PK/PD modeling"
                )
                schedule = st.selectbox(
                    "Dosing Schedule",
                    ["Q4W", "Q8W", "Q12W", "QW", "Q2W"],
                    index=0,
                    help="Dosing frequency for PK/PD exposure calculation"
                )
            else:
                st.warning("⚠️ PK/PD modeling currently optimized for GalNAc conjugates only.")
    
    run_btn = st.button("🔬 Run Full OligoTox Prediction", type="primary",
                        use_container_width=True)'''

content = content.replace(conjugate_section, new_conjugate_section)

# 3. Integrate PK/PD computation after OligoTox prediction (after run_oligotox call)
result_section = '''                # Call production inference
                result = run_oligotox(inp, omics_payload=None)
                
                st.success("✅ Prediction complete!")'''

new_result_section = '''                # Call production inference
                result = run_oligotox(inp, omics_payload=None)
                
                # Apply PK/PD adjustment if enabled (Phase 2A)
                if enable_pkpd and conjugate == "GalNAc":
                    result = adjust_toxicity_risk_with_pkpd(
                        oligotox_result=result,
                        sequence=seq_clean,
                        modality="siRNA",
                        chemistry=chemistry,
                        dose_mg_kg=dose_mg_kg,
                        schedule=schedule
                    )
                    st.success("✅ Prediction complete! (PK/PD-adjusted)")
                else:
                    st.success("✅ Prediction complete!")'''

content = content.replace(result_section, new_result_section)

# 4. Add PK/PD metrics display in results section (after risk scores display)
risk_display = '''                st.markdown("---")
                st.markdown("#### Per-Endpoint Risk Scores")'''

new_risk_display = '''                # Display PK/PD metrics if available
                if enable_pkpd and conjugate == "GalNAc" and "pkpd_metrics" in result:
                    st.markdown("---")
                    st.markdown("#### 🧬 PK/PD Metrics (Phase 2A Bridge)")
                    pkpd_metrics = result["pkpd_metrics"]
                    col_pk1, col_pk2, col_pk3 = st.columns(3)
                    with col_pk1:
                        st.metric(
                            "Liver AUC",
                            f"{pkpd_metrics['liver_auc_pmol_g_h']:.1f} pmol·h/g",
                            help="Total liver exposure over dosing interval"
                        )
                    with col_pk2:
                        st.metric(
                            "Liver Cmax",
                            f"{pkpd_metrics['liver_cmax_pmol_g']:.1f} pmol/g",
                            help="Peak liver concentration"
                        )
                    with col_pk3:
                        st.metric(
                            "Exposure Score",
                            f"{pkpd_metrics['exposure_score']:.2f}",
                            help="Weighted composite exposure metric for risk adjustment"
                        )
                    st.caption(
                        "PK/PD metrics computed using a 2-compartment GalNAc-siRNA model. "
                        "Exposure score feeds into toxicity risk adjustment (multiplicative)." )
                
                st.markdown("---")
                st.markdown("#### Per-Endpoint Risk Scores")'''

content = content.replace(risk_display, new_risk_display)

# Write the updated content
with open('oligotox_app.py', 'w') as f:
    f.write(content)

print("✅ PK/PD UI integration added to oligotox_app.py")
print("Changes made:")
print("  1. Added PK/PD bridge import")
print("  2. Added PK/PD toggle, dose, and schedule inputs in Sequence Analysis tab")
print("  3. Integrated PK/PD adjustment after OligoTox prediction")
print("  4. Added PK/PD metrics display in results section")
