# v0.3 package integration

This repo integrated the downloaded v0.3 package from `downloads/` into the canonical project structure.

## Canonical destinations

- Curated dataset: `data/curated/curated_oligo_tox_master_v3.csv`
- FDA review PDFs: `data/raw/fda_reviews/`
- Extractors and training helpers: `scripts/`
- Labeling docs and model notes: `docs/`
- Gap notebook: `notebooks/`
- Endpoint summary report: `reports/endpoint_summary_v3.csv`

## Recommended next actions

1. Compare v3 schema with the current master dataset before overwriting any existing canonical master file.
2. Run `scripts/papargyri_extractor.py` once Table S1 Excel is downloaded locally.
3. Update `scripts/fda_nda_platelet_extractor.py` with real platelet nadirs from the NDA PDFs.
4. Decide whether `scripts/train_v3.py` becomes the canonical trainer or is merged into `train.py` / `pipelines/`.
5. Regenerate model artifacts only after the unified curated dataset is finalized.
