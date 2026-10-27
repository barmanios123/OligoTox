# OligoTox Training Pipeline Architecture
**Complete End-to-End Documentation**

> **Correction notice (2026-10-27). This document describes the pipeline as it
> was before the label leakage fix, and parts of it document the bug itself.**
>
> Do not follow it as a specification. The three things to know:
>
> 1. **The random 80/20 `train_test_split` shown under "Training Process" is
>    exactly the leakage this project later fixed.** Many rows are replicate
>    measurements of the same oligo, so a random row-level split puts the same
>    sequence in train and test. Folds now come from `StratifiedGroupKFold`
>    grouped on the normalized sequence.
> 2. **Models are not RandomForest.** They are
>    `CalibratedClassifierCV(LogisticRegression(...), method="sigmoid", cv=2)`.
>    The hyperparameter block below is dead code.
> 3. **"54 features" and the per-endpoint prevalences quoted below are wrong.**
>    55 features survive filtering, and hepatocellular injury's base rate is
>    0.126, not 0.0247. 23 outcome-annotation columns and 2 post-dose assay
>    readouts were removed as leakage.
>
> `docs/provenance/DATA_AND_MODEL_PROVENANCE.md` has the corrected pipeline and
> feature breakdown; the root `README.md` has the honest metrics. The feature
> extraction and inference sections here are still broadly accurate.

---

## 🎯 Overview: What You're Training

You're training **13 independent binary classifiers** that predict whether an oligonucleotide will cause specific toxicity endpoints. Each model learns patterns from **54 chemical/structural features** extracted from **1,614 unique oligonucleotides** with known toxicity outcomes.

**Think of it as:**
```
Oligonucleotide Sequence → Feature Extraction → 19 Models → Toxicity Risk Profile
```

---

## 📊 The Data: Input → Labels

### Input Data: `data/processed/training_dataset.csv`
**Dimensions:** 1,614 rows × 54 feature columns

#### Feature Categories (X):

**1. Chemistry Features (10 features)**
- `length`: Nucleotide count (13-30 bases)
- `gc_content`: G+C percentage (0.35-0.80)
- `ps_fraction`: Phosphorothioate backbone ratio (0-1.0)
- `chemical_mod_2ome`: 2'-O-methyl modifications
- `chemical_mod_moe`: 2'-methoxyethyl modifications
- `chemical_mod_lna`: Locked nucleic acid bridges
- `chemical_mod_fana`: 2'-fluoro-arabino modifications
- `gapmer_indicator`: Binary (gapmer vs uniform)
- `conjugate_type`: GalNAc, peptide, lipid, etc.
- `formulation`: Lipid nanoparticle vs naked

**2. Sequence Motifs (8 features)**
- `cpg_island_count`: CpG dinucleotide frequency (immune activation)
- `purine_run_max`: Longest consecutive A/G stretch
- `pyrimidine_run_max`: Longest C/U stretch
- `homopolymer_max`: Any base repeat length
- `gc_skew`: (G-C)/(G+C) strand bias
- `at_skew`: (A-T)/(A+T) strand bias
- `dinucleotide_diversity`: Shannon entropy of 2-mers
- `triplet_complexity`: Unique 3-mer count

**3. Structure Features (6 features)**
- `secondary_structure_mfe`: Minimum free energy (kcal/mol)
- `hairpin_propensity`: Self-complementarity score
- `duplex_stability_tm`: Melting temperature (°C)
- `bulge_loop_count`: Mismatches in structure
- `stem_length_avg`: Average base-paired region
- `structure_ensemble_diversity`: Boltzmann ensemble entropy

**4. Target Context (5 features)**
- `target_gene_expression`: FPKM in hepatocytes
- `target_conservation_score`: PhyloP score (human→mouse)
- `target_splice_site_proximity`: Distance to nearest exon junction
- `target_utr_location`: 5'UTR, CDS, 3'UTR, intergenic
- `target_secondary_structure`: mRNA accessibility (ΔG)

**5. Delivery Parameters (7 features)**
- `dose_mg_kg`: Administered dose
- `route`: IV, SC, oral, inhaled
- `delivery_vehicle`: LNP, naked, conjugate
- `dosing_frequency`: QD, Q3D, QW, etc.
- `treatment_duration_days`: Study length
- `species`: Mouse, rat, NHP, human
- `tissue_distribution_liver`: % of dose in liver

**6. Pharmacokinetics (8 features)**
- `cmax_plasma`: Peak plasma concentration (µg/mL)
- `auc_plasma`: Area under curve (µg·h/mL)
- `half_life_plasma`: t½ in circulation (hours)
- `clearance_ml_min_kg`: Systemic clearance
- `volume_distribution`: Vd (L/kg)
- `protein_binding_fraction`: Albumin/AAG binding
- `renal_excretion_fraction`: % excreted in urine
- `metabolite_formation`: Active metabolite presence

**7. Omics Integration (10 features)**
- `gene_expression_pathway_*`: NF-κB, inflammasome, apoptosis, etc.
- `proteomics_cytokine_score`: IL-6, TNFα, IFNγ levels
- `metabolomics_lipid_perturbation`: Triglyceride/cholesterol shifts
- `transcriptomics_stress_response`: ER stress, oxidative stress genes

---

### Label Matrix: `data/processed/label_matrix.csv`
**Dimensions:** 1,614 rows × 13 endpoints columns

#### The 13 toxicity endpoints (Y):

```python
endpoints = [
    'hepatocellular_injury',          # ALT/AST elevation, necrosis
    'hepatic_vacuolation',            # Lipid droplet accumulation
    'renal_tubular_injury',           # Proximal tubule damage
    'renal_accumulation',             # Kidney cortex oligonucleotide buildup
    'thrombocytopenia',               # Platelet count decrease
    'coagulation_effect',             # PT/aPTT prolongation
    'complement_activation',          # C3a, C5a release
    'innate_immune_activation',       # Macrophage/NK cell response
    'cytokine_release',               # Cytokine storm (CRS)
    'kupffer_cell_activation',        # Liver macrophage activation
    'cns_toxicity',                   # Neurotoxicity (rare)
    'offtarget_gene_silencing',       # Unintended mRNA knockdown
    'injection_site_reaction',        # Local inflammation
    'body_weight_loss',               # >10% weight loss
    'mortality',                      # Dose-limiting toxicity
    'proinflammatory_cytokines',      # IL-1β, IL-6, TNFα
    'tlr_activation',                 # TLR7/8/9 immune stimulation
    'hybridization_related',          # On-target exaggerated pharmacology
    'developmental_toxicity'          # Embryo-fetal effects
]
```

**Encoding:**
- `1` = Toxic event observed (positive)
- `0` = No toxicity / Safe (negative)

**Example Row:**
```
oligo_id,hepatocellular_injury,hepatic_vacuolation,kupffer_cell_activation,...
ASO_001,0,1,1,...  # Safe liver, but vacuolation + KC activation
```

---

## 🤖 Model Architecture: Random Forest Per Endpoint

### Why Random Forest?

1. **Non-linear Feature Interactions:** Captures complex relationships
   - Example: `gc_content × ps_fraction × cpg_count` → TLR activation
2. **Handles Imbalanced Data:** Class weighting compensates for rare toxicities
3. **Feature Importance:** Identifies which features drive predictions
4. **Robust to Overfitting:** Ensemble of 100 trees averages out noise

### Hyperparameters (Fixed Across All Endpoints):

```python
# SUPERSEDED. The shipped model, from _make_model() in
# src/oligotox/training/train.py:
from sklearn.calibration import CalibratedClassifierCV
from sklearn.linear_model import LogisticRegression

base = LogisticRegression(
    max_iter=2000,
    class_weight="balanced",   # upweight the rare positive class
    solver="liblinear",
)
model = CalibratedClassifierCV(base, method="sigmoid", cv=2)
```

With 74 unique sequences, a calibrated linear model is the defensible choice.
Gradient boosting was evaluated as a baseline and matches or beats it on
several endpoints, but every confidence interval overlaps, so there is no
evidence to prefer either (see the root README's Baselines section).

**Class Weighting Formula:**
```
weight_toxic = n_total / (2 × n_toxic)
weight_safe = n_total / (2 × n_safe)
```

Example for `hepatocellular_injury` (2.47% prevalence):
- `weight_toxic` = 1614 / (2 × 40) = 20.18
- `weight_safe` = 1614 / (2 × 1574) = 0.51
→ Model penalizes misclassifying toxic cases **40× more** than safe cases

---

## ⚙️ Training Process: Step-by-Step

### **For Each of the 13 endpoints:**

#### Step 1: Data Loading
```python
import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedGroupKFold, cross_val_predict

# Load features (X)
X = pd.read_csv('data/processed/training_dataset.csv').drop('oligo_id', axis=1)
# 1614 rows × 54 columns

# Load labels (Y) for ONE endpoint at a time
y = pd.read_csv('data/processed/label_matrix.csv')['hepatocellular_injury']
# 1614 binary labels (0/1)
```

#### Step 2: Sequence-grouped cross-validation

**The random split previously documented here was the leakage.** The 1,614
rows cover only 74 unique sequences, so `train_test_split` on rows put
replicate measurements of the same oligo on both sides. Hepatocellular injury
scored 0.93 PR-AUC that way and 0.27 once the split was fixed.

```python
from sklearn.model_selection import StratifiedGroupKFold, cross_val_predict

# All rows of one sequence stay on the same side of every split.
# RNA and DNA spellings of a sequence are the same group.
groups = raw["sequence"].str.upper().str.replace("U", "T", regex=False)

cv = StratifiedGroupKFold(n_splits=n_splits)
prob = cross_val_predict(_make_model(), X, y, cv=cv, groups=groups,
                         method="predict_proba")[:, 1]
```

There is no held-out test set. With 4-27 positive sequences per endpoint, a
single split would be noise, so every published number is out-of-fold under
sequence-grouped CV, reported with its base rate and a 95% bootstrap CI over
sequences.

#### Step 3: Model Training
```python
# Superseded: see Step 2. Fitting now happens inside
# cross_val_predict, once per fold, on the grouped folds.
```

**What Happens Internally:**
1. **Bootstrap Sampling:** Create 100 random subsets of 1,291 oligos
2. **Feature Selection:** Each tree sees √54 ≈ 7 random features
3. **Tree Building:** Recursively split on features that best separate toxic/safe
   - Example split: `IF cpg_count > 5 AND ps_fraction > 0.8 → PREDICT TOXIC`
4. **Ensemble:** Average predictions from all 100 trees

#### Step 4: Evaluation on Test Set
```python
from sklearn.metrics import roc_auc_score, classification_report

# Predict probabilities (0-1 range)
y_pred_proba = model.predict_proba(X_test)[:, 1]  # Column 1 = toxic probability

# Calculate AUROC (Area Under ROC Curve)
auroc = roc_auc_score(y_test, y_pred_proba)
print(f"AUROC: {auroc:.4f}")  # Goal: >0.80 for good model
```

**AUROC Interpretation:**
- **1.0** = Perfect classifier (no errors)
- **0.80-0.99** = Excellent
- **0.70-0.79** = Good
- **0.50-0.69** = Poor (barely better than random)
- **<0.50** = Failed (worse than random guessing)

#### Step 5: Save Model
```python
import joblib

model_path = f'models/hepatocellular_injury_model.pkl'
joblib.dump(model, model_path)
```

**Repeat Steps 1-5 for all 13 endpoints → 19 `.pkl` files**

---

## 📈 Performance Metrics: Understanding Your Output

### What the Final Report Shows:

```
Index | Endpoint                     | AUROC  | Severity | Base Rate
------|------------------------------|--------|----------|----------
0     | hepatocellular_injury        | 0.0123 | High     | 0.0247
9     | kupffer_cell_activation      | 0.1111 | Moderate | 0.2222
```

### Column Breakdown:

1. **AUROC (Primary Metric):**
   - **0.0123 for hepatocellular_injury** → **MODEL FAILED**
   - Why? Only 40/1614 positive examples (2.47%)
   - Test set had only **8 toxic cases** → Not enough signal to learn

2. **Base Rate (Data Characteristic):**
   - Proportion of toxic cases in training data
   - **0.0247 = 2.47%** → Extreme class imbalance
   - **0.2222 = 22.22%** → More balanced, easier to learn

3. **Severity Classification:**
   - **High (<5%):** Rare, severe toxicities (hard to predict)
   - **Moderate (5-25%):** Intermediate frequency
   - **Low (>25%):** Common events (easier to predict)

### Why Are Most AUROC Scores So Low?

**Root Causes:**

1. **Extreme Class Imbalance:**
   - `cns_toxicity`: Only 12/1614 positive cases (0.73%)
   - Test set might have 0-2 positive examples → Impossible to evaluate

2. **Insufficient Features for Rare Events:**
   - 54 features may not capture rare toxicity mechanisms
   - Need mechanistic features: protein binding, cellular uptake, etc.

3. **Small Dataset:**
   - 1,614 oligos is small for ML (ideal: >10K per endpoint)
   - Transfer learning from larger datasets could help

4. **Noise in Labels:**
   - Historical data from multiple studies → Inconsistent toxicity criteria
   - One mislabeled case out of 12 = 8% label noise

---

## 🔮 Prediction at Inference Time

### When You Input a New Sequence in Streamlit:

#### Step 1: Feature Extraction
```python
sequence = "UAGGGCCAAGUGAUCAUGAA"  # User input

# Extract all 54 features
features = extract_features(sequence)
# → [length=20, gc_content=0.50, ps_fraction=1.0, cpg_count=2, ...]
```

#### Step 2: Load All 19 Models
```python
import joblib

models = {}
for endpoint in endpoints:
    model_path = f'models/{endpoint}_model.pkl'
    models[endpoint] = joblib.load(model_path)
```

#### Step 3: Predict for Each Endpoint
```python
predictions = {}
for endpoint, model in models.items():
    prob_toxic = model.predict_proba([features])[0, 1]  # Index [0,1] = toxic probability
    
    predictions[endpoint] = {
        'probability': prob_toxic,
        'prediction': 'TOXIC' if prob_toxic > 0.5 else 'Safe',
        'severity': severity_map[endpoint],
        'confidence': max(prob_toxic, 1-prob_toxic)  # Distance from 0.5
    }
```

#### Step 4: Generate Risk Profile
```python
# Sort by risk (highest first)
sorted_risks = sorted(
    predictions.items(),
    key=lambda x: x[1]['probability'],
    reverse=True
)

# Display
for endpoint, result in sorted_risks[:5]:  # Top 5 risks
    print(f"{endpoint}: {result['probability']:.1%} risk ({result['severity']} severity)")
```

**Example Output:**
```
TOP 5 PREDICTED RISKS:
1. kupffer_cell_activation: 78% risk (Moderate severity) ⚠️
2. thrombocytopenia: 65% risk (Low severity)
3. hepatic_vacuolation: 54% risk (Moderate severity)
4. complement_activation: 42% risk (Moderate severity)
5. injection_site_reaction: 38% risk (High severity)

OVERALL SAFETY SCORE: 22/100 (HIGH RISK)
```

---

## 🛠️ How Decisions Are Made: Feature Importance

### What Drives Predictions?

After training, you can inspect which features matter most:

```python
import matplotlib.pyplot as plt

# Get feature importances
model = joblib.load('models/kupffer_cell_activation_model.pkl')
importances = model.feature_importances_
feature_names = X.columns

# Sort and plot
indices = np.argsort(importances)[::-1][:10]  # Top 10
plt.barh(feature_names[indices], importances[indices])
plt.xlabel('Gini Importance')
plt.title('Top 10 Features for Kupffer Cell Activation')
```

**Typical Results:**
```
TOP FEATURES FOR KUPFFER CELL ACTIVATION:
1. cpg_island_count (0.18)        ← CpG motifs trigger TLR9
2. ps_fraction (0.14)             ← Phosphorothioate backbone immune-stimulatory
3. dose_mg_kg (0.11)              ← Higher dose = more activation
4. gc_content (0.09)              ← GC-rich sequences more immunogenic
5. chemical_mod_2ome (0.08)       ← 2'-OMe reduces activation (negative weight)
6. formulation_lnp (0.07)         ← LNP delivery enhances uptake
7. target_gene_expression (0.06)  ← Highly expressed targets amplify response
8. length (0.05)                  ← Longer oligos = more TLR binding sites
9. route_iv (0.04)                ← IV vs SC affects liver exposure
10. conjugate_galnac (0.04)       ← GalNAc targets hepatocytes
```

**Interpretation:**
- **Positive importance:** Feature increases toxic probability
- **Example:** `cpg_count` ↑ → KC activation ↑
- **Modification impact:** `2'-OMe` modification ↓ → Safer designs

---

## 🧪 Data Flow Summary

```
RAW LITERATURE DATA (Papers, Patents, FDA Filings)
         ↓
CURATION PIPELINE
  ├─ Extract oligonucleotide sequences
  ├─ Annotate chemical modifications
  ├─ Label toxicity outcomes (13 endpoints)
  └─ Aggregate multi-study data
         ↓
FEATURE ENGINEERING
  ├─ Chemistry: GC%, modifications, backbone
  ├─ Motifs: CpG, homopolymers, complexity
  ├─ Structure: MFE, hairpins, Tm
  ├─ Context: Target gene, tissue, PK
  └─ Omics: Expression, pathways, metabolites
         ↓
TRAINING DATASET (1,614 × 54)
         ↓
LABEL MATRIX (1,614 × 19)
         ↓
MODEL TRAINING (19 Random Forests)
  ├─ Split: 80% train, 20% test
  ├─ Fit: Learn feature→toxicity patterns
  ├─ Evaluate: AUROC, precision, recall
  └─ Save: 19 .pkl model files
         ↓
INFERENCE (Streamlit App)
  ├─ User inputs new sequence
  ├─ Extract 54 features
  ├─ Load 19 models
  ├─ Predict 19 probabilities
  └─ Generate risk profile
         ↓
OUTPUT: Toxicity Risk Profile + Safety Recommendations
```

---

## 🔍 Key Takeaways

### What's Working:
- ✅ **Pipeline is operational** (data → features → models → predictions)
- ✅ **13 endpoints trained** (comprehensive toxicity coverage)
- ✅ **Feature engineering complete** (54 mechanistically relevant features)
- ✅ **Inference deployed** (Streamlit app functional)

### What Needs Improvement:
- ⚠️ **AUROC scores are low** (most <0.50)
  - Fix: Add more data (target: 5K+ oligos)
  - Fix: Better features (protein binding, cellular uptake, metabolism)
  - Fix: Advanced algorithms (XGBoost, neural networks, ensembles)
- ⚠️ **Extreme class imbalance** (some endpoints <1% positive)
  - Fix: SMOTE oversampling
  - Fix: Focal loss (from object detection)
  - Fix: Cost-sensitive learning
- ⚠️ **Small dataset** (1,614 total)
  - Fix: Transfer learning from public databases (ChEMBL, PubChem)
  - Fix: Data augmentation (synthetic oligo generation)
  - Fix: Multi-task learning (share representations across endpoints)

### Immediate Next Steps:
1. **Audit feature importance** → Identify which features are useless
2. **Implement SMOTE** → Balance minority classes
3. **Try XGBoost** → Often outperforms Random Forest
4. **Add mechanistic features** → Protein binding, off-target scores, PK parameters
5. **Ensemble models** → Combine RF + XGBoost + Neural Net

---

## 📚 References

- **Random Forest Paper:** Breiman, L. (2001). Machine Learning, 45(1), 5-32.
- **Class Imbalance:** He, H., & Garcia, E. A. (2009). IEEE Trans. Knowl. Data Eng., 21(9), 1263-1284.
- **Oligonucleotide Toxicity:** Crooke, S. T. et al. (2020). Nucleic Acid Ther., 30(1), 1-15.
- **Feature Engineering:** Kuhn, M., & Johnson, K. (2019). Feature Engineering and Selection. CRC Press.

---

**Last Updated:** June 30, 2026  
**Maintainer:** Beshoy Armanios  
**Version:** 1.0.0
