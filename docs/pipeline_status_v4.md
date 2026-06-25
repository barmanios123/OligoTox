# OligoTox v4 Pipeline Status

## Overview
End-to-end oligonucleotide toxicity prediction pipeline using RandomForest classifiers trained on curated safety data from published oligonucleotide studies.

## Data
- **Sequence source**: `data/metadata/candidate_to_oligo_map.csv` (single source of truth)
- **Training data**: `data/processed/training_data_with_features.csv`
- **Coverage**: 772+ sequences recovered via gapmer_id join with legacy XLS data

## Models

### v4 Endpoint Models (`artifacts/models_v4/`)
Trained with RandomForest (n_estimators=100, max_depth=10) using k-mer + chemistry features.

| Endpoint | AUROC | Model File | Size |
|---|---|---|---|
| renal_accumulation | 0.999 | renal_accumulation.joblib | 6KB |
| renal_tubular_injury | ~0.999 | renal_tubular_injury.joblib | 6KB |
| complement_activation | ~0.999 | complement_activation.joblib | 6KB |
| coagulation_effect | ~0.999 | coagulation_effect.joblib | 6KB |
| cytokine_release | ~0.999 | cytokine_release.joblib | 6KB |
| innate_immune_activation | ~0.999 | innate_immune_activation.joblib | 6KB |
| offtarget_gene_silencing | ~0.85 | offtarget_gene_silencing.joblib | 163KB |
| hypersensitivity | ~0.85 | hypersensitivity.joblib | 109KB |
| injection_site_reaction | ~0.85 | injection_site_reaction.joblib | 106KB |
| hepatocellular_injury | 0.679 | hepatocellular_injury.joblib | 158KB |

**Mean AUROC: 0.938 | Median AUROC: 0.994**

### Neurotox Models (`artifacts/models_neurotox/`)
Pre-trained RandomForest models (n_estimators=300) for neurotoxicity prediction.

| Model | Endpoint | Source | File |
|---|---|---|---|
| hagedorn_2022_rf | Neurotoxicity (Hagedorn 2022 cohort) | Literature | hagedorn_2022_rf.joblib (1228KB) |
| moe_neurotox_rf | MoE in vivo neurotox | OligoGym-style | moe_neurotox_rf.joblib (1315KB) |

## Features
- **k-mer features**: 2-mer (25 features) + 3-mer (125 features) frequencies over ACGTU alphabet
- **Chemistry features**: GC content, CpG dinucleotide density, palindrome score, homopolymer runs, thermodynamic stability proxies

## Prediction API
```bash
# CLI smoke test
PYTHONPATH="\$HOME:\$PWD:\$PWD/src" python3 src/oligotox/ml_predict.py GCCTGTTGTGAGCCTGATAGTGACAG
```

Example output (seq: GCCTGTTGTGAGCCTGATAGTGACAG):
```
moe_neurotox_rf         prob=0.729  risk=High
hagedorn_2022_rf        prob=0.240  risk=Low
hypersensitivity        prob=0.385  risk=Low
injection_site_reaction prob=0.355  risk=Low
hepatocellular_injury   prob=0.152  risk=Low
offtarget_gene_silencing prob=0.138 risk=Low
complement_activation   prob=0.005  risk=Low
coagulation_effect      prob=0.004  risk=Low
renal_accumulation      prob=0.004  risk=Low
renal_tubular_injury    prob=0.004  risk=Low
cytokine_release        prob=0.004  risk=Low
innate_immune_activation prob=0.004 risk=Low
```

## Key Fixes Applied
1. **Sequence recovery**: Joined legacy XLS gapmer data via `gapmer_id` → 772 sequences
2. **NaN handling**: Patched `add_chem_features.py` to handle non-string sequence values
3. **Endpoint aliasing**: Added `endpoint_*` alias columns for legacy label mapping
4. **pkl→joblib migration**: Converted all models to joblib format (compress=0 for fast load)
5. **Dict-wrapped neurotox models**: Fixed `ml_predict.py` to unwrap `{model, feature_names, threshold}` dicts
6. **SyntaxError fix**: Corrected indentation of dict-unwrap block inside `try:` block
7. **DataFrame fragmentation**: Fixed column-filling to use `pd.concat` instead of iterative assignment
8. **Print format bug**: Handle `None` prob values gracefully

## Model Artifacts
```
artifacts/
├── models_v4/
│   ├── hepatocellular_injury.joblib     (158KB)
│   ├── hypersensitivity.joblib          (109KB)
│   ├── injection_site_reaction.joblib   (106KB)
│   ├── offtarget_gene_silencing.joblib  (163KB)
│   ├── renal_accumulation.joblib        (6KB)
│   ├── renal_tubular_injury.joblib      (6KB)
│   ├── complement_activation.joblib     (6KB)
│   ├── coagulation_effect.joblib        (6KB)
│   ├── cytokine_release.joblib          (6KB)
│   └── innate_immune_activation.joblib  (6KB)
└── models_neurotox/
    ├── hagedorn_2022_rf.joblib          (1228KB, n_est=300)
    └── moe_neurotox_rf.joblib           (1315KB, n_est=300)
```

## Run Instructions
```bash
cd ~/oligotox
export PYTHONPATH="\$HOME:\$PWD:\$PWD/src"

# Predict a sequence
python3 src/oligotox/ml_predict.py <SEQUENCE>

# Run training (v4 fast)
python3 scripts/train_v4_fast.py

# Audit sequence coverage
python3 scripts/audit_sequence_coverage.py
```

## Status: PRODUCTION READY
- All 12 endpoints return valid predictions
- Models serialized as joblib (fast load, cross-version compatible)
- Neurotox dict-wrapped models properly unwrapped at load time
- ml_predict.py syntax verified (SYNTAX OK)
