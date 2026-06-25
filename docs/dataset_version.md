# Dataset versioning

## Current curated pack
The current curated dataset release is `curated_oligo_tox_master_v2.csv`.

Associated v2 assets:
- `data/processed/curated_oligo_tox_master_v2.csv`
- `data/processed/endpoint_summary_v2.csv`
- `data/catalog/molecule_catalog_v2.json`
- `data/sequences/oligo_sequences_v2.fasta`
- `docs/labeling_scheme_v2.csv`
- `docs/curation_priority_v2.csv`

## How it flows into training
The v2 curated table is used upstream by `scripts/curated_csv_to_training.py` and the central dataset build pipeline.

The trainer itself does not fit directly on the raw curated CSV. Instead:
1. curated data is normalized into processed training inputs,
2. `pipelines.build_training_dataset` generates `feature_matrix.csv`, `label_matrix.csv`, and metadata,
3. `src/oligotox/training/train.py` trains models from those processed matrices.

## Why training rows exceed curated rows
The curated v2 pack contains 18 compounds, but the final training matrices can contain many more rows because the build pipeline expands and combines processed sources before model fitting.
