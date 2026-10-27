# OligoTox Data Directory

## What is in the repo

Four files, 1.2 MB total, are committed so every number published in the root
README can be reproduced from a clone:

```
data/processed/
├── feature_matrix.csv      1,614 rows x 55 engineered features
├── label_matrix.csv        13 binary toxicity endpoints
├── training_dataset.csv    sequences, source dataset and curated labels
└── dataset_metadata.json   column manifest and training hash
```

Everything else under `data/` is a local intermediate or a superseded variant
and is gitignored: `data/raw/` (160 MB) and the rest of `data/processed/`
(445 MB) are not distributed.

## Dataset

- **1,614 rows covering 74 unique sequences.** Many rows are replicate
  measurements of the same oligo, which is why all cross-validation is grouped
  by normalized sequence. Treat 74, not 1,614, as the sample size.
- **55 features**, all derivable from a design at scoring time. The curated
  `label_*` / `endpoint_*` annotations and the post-dose assay readouts
  (`cas_avg`, `kd_avg`) were removed as label leakage; see the root README.
- **13 endpoints.** Only hepatocellular injury survives the baseline
  comparison; 9 endpoints have 7 or fewer positive sequences.

### Sources

`source` in `training_dataset.csv` names where each row came from. Two
systematic sources supply 1,572 rows (`oligotox_main` 800,
`papargyri_2020` 772). The remaining 42 rows are individually curated from
papers, PMC articles and FDA labels, and 36 of those 42 carry a positive
label, because they were curated on the strength of a reported finding. That
selection bias is measured in the root README's Baselines section, not
assumed away.

Full citations for the upstream datasets are in
`docs/provenance/DATA_AND_MODEL_PROVENANCE.md`.

## Reproducing

```bash
python -m oligotox.training.train       # ~10 s, retrains all 13 models
python scripts/evaluate_baselines.py    # ~95 s, baselines and bootstrap CIs
python scripts/derive_evidence_tiers.py # rewrites the evidence tier config
```

Verified in a clean clone: the metrics, the baseline CSVs, the tier config and
all 13 `.joblib` models come out identical to the committed ones.
