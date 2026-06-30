# OligoTox Data Directory

## Structure

```
data/
├── raw/              # Original curated datasets (not in repo - too large)
├── interim/          # Intermediate processing steps
└── processed/        # Final training-ready datasets
    ├── training_dataset.csv      # 1,614 × 54 features
    └── label_matrix.csv          # 1,614 × 19 toxicity endpoints
```

## Data Description

### Training Dataset (`training_dataset.csv`)
- **1,614 oligonucleotides** with **54 engineered features**
- Features span chemistry, structure, delivery, pharmacokinetics, and omics
- See `docs/TRAINING_PIPELINE_ARCHITECTURE.md` for full feature list

### Label Matrix (`label_matrix.csv`)
- **19 binary toxicity endpoints** (1=toxic, 0=safe)
- Endpoints cover hepatotoxicity, nephrotoxicity, immunogenicity, and more
- Data curated from 200+ publications, patents, and FDA filings

## Data Availability

**Note:** Raw data files are not included in this repository due to size constraints and potential proprietary information. The processed datasets used for training are available upon request for academic/research purposes.

For data inquiries: [Contact Information]
