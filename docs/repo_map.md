# Repository map

## Canonical entrypoints
- `oligotox_app.py` — Streamlit application
- `train.py` — training entrypoint
- `src/oligotox/inference/api.py` — inference API
- `src/oligotox/training/train.py` — training implementation

## Core directories
- `src/` — package namespace and production code
- `data/` — curated, processed, raw, catalog, sequence, and versioned datasets
- `docs/` — dataset notes, labeling scheme, curation priorities, architecture notes
- `scripts/` — operational scripts for ingestion, curation, normalization, and maintenance
- `tests/` — automated tests
- `artifacts/` — trained models, metrics, configs, and generated outputs
- `archive/` — deprecated, historical, or staging content retained for provenance

## Scientific engine layers
- `rules/` — mechanistic toxicity rules
- `translation/` — translational bridge and species mapping
- `integrations/` — seed match, TLR, GSEA, phylogenetic, and context integrations
- `omics/` — omics ingestion, QC, reduction, mapping, and mechanistic joins
- `features/` — feature engineering and matrix building
- `pipelines/` — dataset and training orchestration

## Presentation note
For external presentation, treat OligoTox as a single platform with one app, one training path, one dataset lineage, and one active release version.
