# OligoTox: Mechanism-Aware Toxicity Prediction for Oligonucleotide Therapeutics

**Public Case Study Repository**

> This repository showcases the architecture, methodology, and capabilities of OligoTox — a clinical prediction platform for oligonucleotide safety decision support. The full implementation, including curated training data and proprietary scoring models, is maintained privately.

---

## Overview

OligoTox is a mechanism-aware toxicity prediction and safety decision-support platform designed specifically for oligonucleotide therapeutics, including:

- **ASOs** (antisense oligonucleotides)
- **siRNAs** (small interfering RNAs)
- **lncRNA-targeting oligos**
- **mRNA therapies**
- **Other RNA-based modalities**

### Core Innovation: The D1/D2/D3 Framework

OligoTox organizes toxicity risk into three mechanistic categories:

| Category | Mechanism | Examples |
|----------|-----------|----------|
| **D1: On-Target Pharmacology** | Exaggerated intended effects | Target knockdown → pathway disruption → organ toxicity |
| **D2: Hybridization Off-Target** | Unintended RNA binding via sequence complementarity | microRNA/lncRNA disruption, immune activation (CpG, RIG-I) |
| **D3: Chemistry/Delivery** | Chemical backbone or conjugate liabilities | PS-oligo pro-coagulation, GalNAc hepatotoxicity, LNP inflammation |

This framework enables:
- **Mechanistic attribution** for observed toxicities
- **Assay and biomarker recommendations** tied to D1/D2/D3 liabilities
- **Translational bridging** from preclinical to clinical safety signals

---

## Platform Capabilities

### 1. **Candidate Profiler**
Enter an oligonucleotide sequence + chemistry metadata → receive endpoint-specific toxicity risk scores with confidence bands.

**Inputs:**
- Sequence (RNA/DNA)
- Modality (ASO, siRNA, etc.)
- Chemistry (PS, 2'-OMe, LNA, mixed)
- Conjugate (GalNAc, cholesterol, etc.)

**Outputs:**
- Per-endpoint risk scores (hepatotox, nephrotox, thrombocytopenia, etc.)
- Confidence levels (Low/Medium/High)
- Dominant toxicity mechanism
- D1/D2/D3 attribution

### 2. **Candidate Comparator**
Side-by-side risk profiling for multiple oligo candidates to guide lead selection.

### 3. **Assay Planner**
Recommend assays and biomarkers based on predicted D1/D2/D3 liabilities.

### 4. **Off-Target Evaluator**
Evaluate sequence-dependent off-target risks (microRNA seed matches, immunostimulatory motifs, etc.).

### 5. **Translational Bridge**
Map preclinical toxicity findings to relevant clinical adverse events and biomarkers.

---

## Technical Architecture

### Stack
- **Backend:** Python 3.10+, FastAPI, Pydantic
- **ML:** scikit-learn, NumPy, pandas
- **Data:** SQLite (dev), PostgreSQL (prod)
- **UI:** Streamlit (operator console), React (planned dashboard)
- **Evaluation:** AUROC, AUPRC, Brier score, calibration curves
- **Leakage Control:** DataSAIL-based splitting (sequence similarity, batch effects, temporal structure)

### Module Structure

```
oligotox/
├── src/oligotox/
│   ├── inference/          # Production prediction engine
│   │   ├── api.py          # run_oligotox() main entry point
│   │   └── scoring.py      # Endpoint-specific scoring logic
│   ├── training/           # Model training pipeline
│   │   ├── train.py        # Main training script
│   │   └── splits.py       # DataSAIL leakage-aware splitting
│   ├── evaluation/         # Clinical prediction metrics (NEW)
│   │   ├── metrics.py      # AUROC, Brier, calibration
│   │   └── reports.py      # Evaluation report generation
│   ├── omics/              # Mechanistic feature engineering
│   ├── config/             # Paths, constants, schemas
│   └── translation/        # Preclinical-to-clinical bridging
├── oligotox_app.py         # Streamlit operator console
├── api.py                  # FastAPI REST endpoint
├── train.py                # Training pipeline entry point
└── data/processed/         # Curated training dataset (private)
```

---

## Clinical Prediction Workflow

### 1. **Data Curation & Labeling**
- Curated corpus from public literature, FDA labels, clinical trial reports
- Endpoint-specific binary labels for 15+ toxicity categories
- Sequence-aware metadata (chemistry, dose, species, route)

### 2. **Leakage-Resistant Splitting**
Using **DataSAIL** to prevent:
- Sequence overlap between train/val/test
- Batch/assay confounding
- Temporal leakage (train on old data, test on new)
- Study family effects

### 3. **Feature Engineering**
- Sequence-derived: length, GC%, CpG count, seed region analysis
- Chemistry: backbone type, modification pattern, conjugate
- Target biology: pathway, tissue distribution, turnover rate
- Off-target: miRNA seed matches, RIG-I motifs, protein binding

### 4. **Model Training**
- Multi-task learning across 15+ endpoints
- Per-endpoint risk band thresholds
- Uncertainty quantification via ensemble or calibration

### 5. **Evaluation**
- **AUROC & AUPRC:** Discrimination performance
- **Brier Score:** Calibration quality
- **Calibration Curves:** Reliability diagrams
- **Per-Endpoint Reports:** Threshold selection, confusion matrices

### 6. **Production Inference**
```python
from oligotox.inference.api import run_oligotox

inp = {
    "sequence": "AUGCUAGCUAGCUAGC",
    "modality": "ASO",
    "chemistry": "PS",
    "conjugate": "GalNAc",
}

result = run_oligotox(inp, omics_payload=None)

# Returns:
# {
#   "per_endpoint": {
#     "hepatotoxicity": {"prob": 0.72, "risk_band": "high", ...},
#     "thrombocytopenia": {"prob": 0.15, "risk_band": "low", ...},
#     ...
#   },
#   "dominant_endpoint": "hepatotoxicity",
#   "D1": 0.45, "D2": 0.30, "D3": 0.25
# }
```

---

## Evaluation Summary (Representative Metrics)

| Endpoint | AUROC | AUPRC | Brier | Threshold |
|----------|-------|-------|-------|----------|
| Hepatotoxicity | 0.78 | 0.65 | 0.18 | 0.42 |
| Nephrotoxicity | 0.81 | 0.58 | 0.16 | 0.38 |
| Thrombocytopenia | 0.74 | 0.52 | 0.21 | 0.35 |
| Immunostimulation | 0.86 | 0.71 | 0.14 | 0.50 |
| **Composite** | **0.80** | **0.62** | **0.17** | — |

*Note: These are representative synthetic metrics for demonstration. Actual performance is proprietary.*

---

## Why This Matters for Oligonucleotide Drug Development

1. **Early Risk Flagging:** Predict toxicity liabilities before expensive in vivo studies
2. **Rational Lead Selection:** Compare candidates on mechanism-informed risk profiles
3. **Assay Optimization:** Focus safety testing on highest-risk endpoints
4. **Regulatory Preparedness:** Anticipate FDA/EMA questions with mechanistic rationale
5. **Clinical Trial Design:** Plan monitoring strategy based on predicted biomarkers

---

## Public vs. Private Components

### Public (This Repo)
- Platform architecture and design philosophy
- D1/D2/D3 mechanistic framework
- Example input/output schemas
- Evaluation methodology
- Code structure and module organization
- Representative (redacted) metrics

### Private (Full Implementation)
- Curated training dataset (literature, FDA labels, trial reports)
- Endpoint-specific labeling rules
- Full feature engineering pipeline
- Trained model binaries
- Proprietary mechanistic scoring logic
- Production deployment infrastructure

---

## Use Cases

### For Drug Discovery Teams
- Triage oligo candidates by safety profile
- Optimize chemistry and conjugate selection
- Plan preclinical safety packages

### For Regulatory Affairs
- Generate mechanistic hypotheses for observed toxicities
- Identify relevant biomarkers for clinical monitoring
- Support IND/NDA submissions with computational predictions

### For Clinical Development
- Design safety monitoring plans
- Interpret adverse events with D1/D2/D3 framework
- Guide dose escalation strategies

---

## Roadmap

- [x] Core D1/D2/D3 framework
- [x] Multi-endpoint toxicity scoring
- [x] Leakage-resistant evaluation
- [x] Streamlit operator console
- [ ] React production dashboard
- [ ] REST API deployment
- [ ] Circulating biomarker module
- [ ] CNS-specific endpoint models
- [ ] DART (Antibody-oligo conjugate) module

---

## Contact & Collaboration

This platform is actively maintained and under continuous development. For inquiries about:
- **Collaboration** on oligonucleotide safety research
- **Consulting** for biopharma oligo programs
- **Licensing** or partnership opportunities

Please reach out via [GitHub Issues](https://github.com/barmanios123/oligotox-case-study/issues) or LinkedIn.

---

## License

This public case-study repository is provided for **portfolio and demonstration purposes only**. The full implementation, including curated data and trained models, is proprietary and not open-source.

**© 2026 Beshoy Armanios. All rights reserved.**

---

## References & Inspiration

- FDA Oligonucleotide Safety Guidance (2023)
- EMA Reflection Paper on Oligonucleotide-Based Therapies (2024)
- NIH OligoTox Challenge (2024)
- DataSAIL: Leakage-Free Data Splitting (2023)
- Clinical prediction best practices: calibration, uncertainty, interpretability

