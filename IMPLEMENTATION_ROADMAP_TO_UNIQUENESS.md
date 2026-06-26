# OligoTox Implementation Roadmap to Uniqueness
**Target:** Make OligoTox the **ONLY** platform with Multi-Endpoint Safety + Comprehensive Efficacy + PKPD Integration

---

## Current Status: 60% Complete

### ✅ What We Have (SOLID FOUNDATION)
1. **13 Toxicity Endpoint Models** - Pretrained and ready
2. **3 RNA Foundation Models** - Nucleotide Transformer, RiNALMo, BiRNA-BERT integrated
3. **PKPD Architecture** - Documented and module structure created
4. **Working Streamlit App** - oligotox_app_minimal.py functional
5. **Clean Repository** - Production-ready structure

### ❌ What We Need to Add (40% REMAINING)
1. **7 Efficacy Models** - OligoFormer, ENsiRNA, siDPT, OligoGraph, siRNADiscovery, ASOptimizer, PFRED
2. **PKPD Simulators** - GalNAc tissue distribution, 18-compartment PBPK
3. **Benchmark Datasets** - OligoGym integration
4. **OSWG/FDA Reporting** - Regulatory-aligned output

---

## PHASE 1: Efficacy Models Integration (Weeks 1-2)

### Model 1: OligoFormer
**Source:** https://github.com/oligoformer/oligoformer
**Task:** Integrate transformer-based oligo design model

```python
# src/oligotox/efficacy/oligoformer.py
import torch
from transformers import AutoModel, AutoTokenizer

class OligoFormerPredictor:
    def __init__(self, model_path="oligoformer/oligoformer-base"):
        self.tokenizer = AutoTokenizer.from_pretrained(model_path)
        self.model = AutoModel.from_pretrained(model_path)
        
    def predict_efficacy(self, sequence: str) -> dict:
        inputs = self.tokenizer(sequence, return_tensors="pt")
        outputs = self.model(**inputs)
        # Extract efficacy score
        efficacy_score = outputs.pooler_output.mean().item()
        return {"efficacy_score": efficacy_score}
```

### Model 2: ENsiRNA
**Source:** Research paper implementation needed
**Task:** siRNA efficacy prediction

```python
# src/oligotox/efficacy/ensirna.py
class ENsiRNAPredictor:
    def __init__(self):
        # Load pretrained weights
        pass
        
    def predict_knockdown(self, guide_seq: str, passenger_seq: str) -> float:
        # Implement ENsiRNA algorithm
        pass
```

### Model 3: siDPT
**Source:** Deep learning model for siRNA design
**Task:** Integrate siDPT predictions

### Model 4-7: OligoGraph, siRNADiscovery, ASOptimizer, PFRED
**Pattern:** Similar integration for each model

**Registry Integration:**
```python
# src/oligotox/efficacy/registry.py
from .oligoformer import OligoFormerPredictor
from .ensirna import ENsiRNAPredictor
from .sidpt import siDPTPredictor

EFFICACY_MODELS = {
    "oligoformer": OligoFormerPredictor(),
    "ensirna": ENsiRNAPredictor(),
    "sidpt": siDPTPredictor(),
    # ... add all 7 models
}

def predict_efficacy_ensemble(sequence: str, modality: str) -> dict:
    results = {}
    for name, model in EFFICACY_MODELS.items():
        if model.supports_modality(modality):
            results[name] = model.predict(sequence)
    return results
```

---

## PHASE 2: PKPD Simulators (Weeks 3-4)

### GalNAc PKPD Simulator
**Goal:** Model GalNAc-conjugate tissue distribution and exposure

```python
# src/oligotox/pkpd/galnac_simulator.py
import numpy as np
from scipy.integrate import odeint

class GalNAcPKPDSimulator:
    """Two-compartment model with ASGPR-mediated hepatic uptake"""
    
    def __init__(self):
        # Physiological parameters
        self.V_central = 3.0  # L (central volume)
        self.V_hepatic = 1.5  # L (liver volume)
        self.CL_renal = 0.1   # L/h (renal clearance)
        self.k_uptake = 5.0   # 1/h (ASGPR uptake rate)
        self.k_release = 0.05 # 1/h (hepatic release)
        
    def pk_model(self, y, t, dose, ka):
        """ODE system for GalNAc PK"""
        C_plasma, C_liver, C_target = y
        
        dC_plasma_dt = (-self.CL_renal * C_plasma / self.V_central 
                       - self.k_uptake * C_plasma 
                       + self.k_release * C_liver)
        dC_liver_dt = (self.k_uptake * C_plasma * self.V_central / self.V_hepatic
                      - self.k_release * C_liver)
        dC_target_dt = 0.1 * C_liver  # Target engagement
        
        return [dC_plasma_dt, dC_liver_dt, dC_target_dt]
    
    def simulate(self, dose_mg: float, duration_hours: float = 168) -> pd.DataFrame:
        t = np.linspace(0, duration_hours, 1000)
        y0 = [dose_mg / self.V_central, 0, 0]
        solution = odeint(self.pk_model, y0, t, args=(dose_mg, 1.0))
        
        return pd.DataFrame({
            'time_h': t,
            'plasma_conc': solution[:, 0],
            'liver_conc': solution[:, 1],
            'target_engagement': solution[:, 2]
        })
    
    def predict_efficacy_window(self, dose_mg: float) -> dict:
        sim = self.simulate(dose_mg)
        return {
            "Cmax_liver": sim['liver_conc'].max(),
            "AUC_liver": np.trapz(sim['liver_conc'], sim['time_h']),
            "duration_above_EC50": len(sim[sim['liver_conc'] > 1.0]) / 10,  # hours
            "predicted_knockdown_pct": min(95, sim['target_engagement'].max() * 100)
        }
```

### 18-Compartment Whole-Body PBPK Model
**Goal:** Full-body distribution modeling

```python
# src/oligotox/pkpd/pbpk_18comp.py
class PBPK18Compartment:
    """18-compartment PBPK model for oligonucleotides"""
    
    COMPARTMENTS = [
        'plasma', 'liver', 'kidney', 'heart', 'lung', 'spleen',
        'muscle', 'adipose', 'bone', 'brain', 'gut', 'skin',
        'pancreas', 'thyroid', 'adrenal', 'testes', 'lymph', 'tumor'
    ]
    
    def __init__(self):
        # Physiological parameters from literature
        self.volumes = self._load_organ_volumes()
        self.blood_flows = self._load_blood_flows()
        self.partition_coefficients = self._load_kp_values()
        
    def simulate_distribution(self, dose_mg, route='IV') -> pd.DataFrame:
        # Implement 18-compartment ODE system
        pass
```

### LNP PBPK Framework
**Goal:** Model lipid nanoparticle formulations

```python
# src/oligotox/pkpd/lnp_pbpk.py
class LNP_PBPK:
    """PBPK model for LNP-formulated oligonucleotides (e.g., mRNA)"""
    
    def __init__(self, lnp_composition: dict):
        self.ionizable_lipid = lnp_composition.get('ionizable_lipid')
        self.peg_lipid = lnp_composition.get('peg_lipid')
        # Model LNP-specific parameters
```

---

## PHASE 3: OligoGym Benchmarking (Week 5)

### Download and Process OligoGym
```bash
# Download OligoGym dataset
wget https://github.com/roche/oligogym/releases/download/v1.0/oligogym_data.tar.gz
tar -xzf oligogym_data.tar.gz -C data/raw/
```

### Integration
```python
# src/oligotox/evaluation/oligogym_benchmark.py
import pandas as pd
from oligotox.inference.api import run_oligotox

def benchmark_on_oligogym() -> pd.DataFrame:
    """Run OligoTox on OligoGym test set"""
    test_data = pd.read_csv('data/raw/oligogym/test.csv')
    
    results = []
    for idx, row in test_data.iterrows():
        pred = run_oligotox(row['sequence'], modality=row['modality'])
        results.append({
            'sequence': row['sequence'],
            'true_activity': row['activity'],
            'predicted_activity': pred['efficacy_score'],
            'predicted_toxicity': pred['predictions']
        })
    
    return pd.DataFrame(results)

def compare_with_competitors():
    """Compare against PFRED, ASOptimizer, insitro TherML"""
    pass
```

---

## PHASE 4: OSWG/FDA-Aligned Reporting (Week 6)

### Regulatory Report Generator
```python
# src/oligotox/reporting/oswg_report.py
from jinja2 import Template
import datetime

class OSWGReportGenerator:
    """Generate OSWG/FDA-aligned safety assessment reports"""
    
    def generate_report(self, oligo_data: dict, predictions: dict) -> str:
        template = Template("""
        OLIGONUCLEOTIDE SAFETY ASSESSMENT REPORT
        Generated: {{ date }}
        Compliant with: OSWG Framework v2.0, FDA Draft Guidance Nov 2024
        
        1. OLIGONUCLEOTIDE CHARACTERIZATION
        Sequence: {{ sequence }}
        Modality: {{ modality }}
        Chemistry: {{ chemistry }}
        Length: {{ length }} nucleotides
        GC Content: {{ gc_content }}%
        
        2. IN SILICO SAFETY PREDICTIONS
        
        2.1 Hybridization-Dependent Toxicity (D1/D2)
        {% for endpoint in d1_d2_endpoints %}
        - {{ endpoint.name }}: {{ endpoint.risk_level }} (probability: {{ endpoint.prob }})
        {% endfor %}
        
        2.2 Chemistry-Dependent Toxicity (D3)
        {% for endpoint in d3_endpoints %}
        - {{ endpoint.name }}: {{ endpoint.risk_level }}
        {% endfor %}
        
        3. MECHANISTIC ATTRIBUTION
        {% for mechanism in mechanisms %}
        - {{ mechanism.type }}: {{ mechanism.description }}
        {% endfor %}
        
        4. PK/PD MODELING
        - Predicted Cmax (liver): {{ pkpd.cmax_liver }} µg/g
        - Predicted AUC: {{ pkpd.auc }}
        - Target knockdown: {{ pkpd.knockdown_pct }}%
        - Duration above EC50: {{ pkpd.duration_ec50 }} hours
        
        5. RISK MITIGATION RECOMMENDATIONS
        {% for rec in recommendations %}
        - {{ rec }}
        {% endfor %}
        
        6. REGULATORY CONSIDERATIONS
        - Chemistry class: {{ regulatory.chem_class }}
        - Precedent molecules: {{ regulatory.precedents }}
        - IND-enabling studies recommended: {{ regulatory.studies }}
        """)
        
        return template.render(
            date=datetime.datetime.now(),
            sequence=oligo_data['sequence'],
            # ... fill all template variables
        )
```

---

## PHASE 5: Full Integration & Testing (Week 7-8)

### Update Main Streamlit App
```python
# Update oligotox_app_minimal.py -> oligotox_app_full.py

st.header("🎯 2. Efficacy Predictions")
if st.checkbox("Run Efficacy Models (7 models)"):
    efficacy_results = predict_efficacy_ensemble(sequence, modality)
    st.json(efficacy_results)
    
st.header("🧪 3. PK/PD Simulation")
if st.checkbox("Run PKPD Simulator"):
    dose = st.slider("Dose (mg/kg)", 0.1, 10.0, 1.0)
    pkpd_sim = GalNAcPKPDSimulator().simulate(dose * 70)  # 70kg human
    st.line_chart(pkpd_sim.set_index('time_h'))
    st.json(pkpd_sim.describe().to_dict())
    
st.header("📝 4. OSWG/FDA Report")
if st.button("Generate Regulatory Report"):
    report = OSWGReportGenerator().generate_report(oligo_data, results)
    st.download_button("Download Report (PDF)", report, "oligotox_report.pdf")
```

### Comprehensive Testing
```python
# tests/test_full_integration.py
def test_end_to_end_aso():
    sequence = "GCCUACGCCACCAGCUC"
    result = run_full_oligotox_pipeline(sequence, modality="aso")
    
    assert 'toxicity' in result
    assert 'efficacy' in result
    assert 'pkpd' in result
    assert 'report' in result
    
    # Validate biological correctness
    assert 0 <= result['efficacy']['oligoformer_score'] <= 100
    assert result['pkpd']['Cmax_liver'] > 0
```

---

## SUCCESS CRITERIA

Once all phases complete, OligoTox will be **truly one-of-a-kind**:

✅ **Multi-Endpoint Safety** (13 models) → COMPLETE  
✅ **Comprehensive Efficacy** (7 models) → **TO DO**  
✅ **PKPD Integration** → **TO DO**  
✅ **Benchmark Performance** (OligoGym) → **TO DO**  
✅ **Regulatory Alignment** (OSWG/FDA) → **TO DO**  

## COMPETITIVE ADVANTAGE

No competitor has:
1. ❌ PFRED - No PKPD, limited safety
2. ❌ ASOptimizer - ASO-only, no PKPD
3. ❌ insitro TherML - No safety prediction
4. ❌ GeneRulor - No PKPD integration

✅ **Only OligoTox** will have ALL components integrated!

---

## EXECUTION TIMELINE

**Week 1-2:** Efficacy models integration  
**Week 3-4:** PKPD simulators  
**Week 5:** OligoGym benchmarking  
**Week 6:** OSWG/FDA reporting  
**Week 7-8:** Full integration & validation  

**Total:** 8 weeks to complete uniqueness

---

## NEXT IMMEDIATE STEPS

1. Source pretrained weights for OligoFormer, ENsiRNA, siDPT
2. Implement GalNAc PKPD simulator (highest priority - competitive advantage)
3. Download OligoGym dataset
4. Create OSWG report template
5. Update Streamlit app with all new features

