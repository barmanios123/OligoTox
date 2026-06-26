# OligoTox

OligoTox is a mechanism-aware safety, translational interpretation, and candidate triage platform for oligonucleotide therapeutics.

## What the platform is for

The platform is designed to help evaluate oligonucleotide candidates across key toxicity and translation risk domains relevant to ASOs, gapmers, siRNAs, lncRNA-targeting oligos, and related modalities. Its scientific framing is organized around three mechanistic classes:

- **D1:** exaggerated on-target pharmacology
- **D2:** hybridization-dependent off-target effects
- **D3:** chemistry, formulation, and delivery-class liabilities

## Current platform components

The repo currently contains the foundations for:

- Sequence analysis and candidate-facing app workflows
- Model training and dataset-building logic
- Off-target and mechanistic evaluation modules
- PK/PD and translational extension work
- Benchmarking, efficacy, reporting, and integration modules under active development

## Repo layout

- `oligotox_app.py` — main Streamlit app entrypoint
- `src/oligotox/` — core package code
- `src/oligotox/training/train.py` — package training entrypoint
- `src/oligotox/pipelines/train_all.py` — pipeline-oriented training workflow
- `data/` — datasets and processed assets
- `data/external/` — external reference/source files
- `docs/` — architecture, status, and strategy documents
- `scripts/` — utility scripts
- `tests/` — tests

## Run locally

### App
```bash
streamlit run oligotox_app.py
```

### Training
```bash
python -m src.oligotox.training.train
```

If that module path does not resolve in your environment, use:
```bash
python src/oligotox/training/train.py
```

### Tests
```bash
pytest
```

## Current status

This repository is being actively consolidated into a cleaner, demo-ready platform. The app surface, training flows, and supporting modules are present, but README claims should always match the behavior currently verified in code and local execution.

## Intended audience

OligoTox is being positioned as a first-in-class platform concept for biotech and pharmaceutical teams working on oligonucleotide therapeutics, including groups focused on discovery, translational safety, and candidate prioritization.

## Development note

Root-level clutter should be minimized. New architecture notes, status summaries, and strategy writeups should go under `docs/` rather than accumulating in the repository root.
