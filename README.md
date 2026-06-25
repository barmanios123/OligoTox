# OligoTox v1.0

Mechanism-aware toxicity prediction and safety decision support for RNA and oligonucleotide therapeutics.

## Overview

OligoTox is a research platform for profiling candidate oligonucleotide safety risk using curated data, model-based endpoint prediction, mechanistic rule layers, off-target analysis, assay planning, and translational support views.

The current repository reflects the cleaned **OligoTox v1.0** codebase, with production-oriented code centered on the active app, package modules under `src/`, and model/data artifacts under `artifacts/` and `data/`.

## Current status

Current active capabilities include:

- Streamlit application for candidate scoring and decision support
- Model-backed inference through the unified scoring path
- Mechanistic evidence layers, including seed-match and innate immune context
- Candidate comparison and assay planning support
- Translational bridge views
- PDF safety brief generation when local PDF dependencies are available 

## Canonical entrypoints

Primary entrypoints:

- `oligotox_app.py` — Streamlit application
- `api.py` — top-level inference API shim
- `train.py` — top-level training entrypoint
- `src/oligotox/inference/api.py` — package inference implementation
- `src/oligotox/training/train.py` — package training implementation 

## Repository layout

Primary active areas:

- `src/` — production package code
- `data/` — curated, processed, raw, catalog, and sequence datasets
- `artifacts/` — trained models, configs, metrics, and evaluation outputs
- `integrations/` — mechanistic and external integration logic
- `rules/` — mechanistic and decision-rule engines
- `translation/` — translational support logic
- `assays/` — assay planning support
- `docs/` — dataset notes, architecture notes, and reference documentation
- `tests/` — automated tests
- `archive/` — historical, transitional, and non-canonical materials retained for provenance 

## Run the app

Start the app from the repository root:

```bash
cd ~/oligotox
source .venv/bin/activate
PYTHONPATH=$HOME:$PWD streamlit run oligotox_app.py
```

The app supports:

- Candidate profiling
- Mechanistic evidence review
- Off-target evaluation
- Candidate comparison
- Assay planning support
- Translational bridge views
- PDF safety brief generation when local dependencies are available 

## Training

Run training from the repository root:

```bash
cd ~/oligotox
source .venv/bin/activate
PYTHONPATH=$HOME:$PWD python train.py
```

`train.py` delegates to the package training implementation under `src/oligotox/training/train.py`. 

## Data and artifacts

Important data and artifact areas include:

- `data/processed/`
- `data/catalog/`
- `data/sequences/`
- `artifacts/models/`
- `artifacts/models_v3/`
- `artifacts/config/`
- `artifacts/eval/` 

See `docs/dataset_version.md` for dataset lineage and version details. 

## Recommended workflow

Typical workflow:

1. Launch the app and validate candidate scoring behavior.
2. Review mechanistic evidence, Off-Target outputs, and translational support views.
3. Retrain models with `train.py` when training inputs or labeling assets change.
4. Re-check app behavior after retraining. 

Example:

```bash
cd ~/oligotox
source .venv/bin/activate
PYTHONPATH=$HOME:$PWD python train.py
PYTHONPATH=$HOME:$PWD streamlit run oligotox_app.py
```

## Notes on legacy material

Historical prototypes, archived workspace modules, and transitional files are retained under `archive/` for provenance. These are not the primary runtime path for the current release. 

## Documentation

Useful supporting documentation includes:

- `docs/dataset_version.md`
- `docs/repo_map.md`

---

## License & IP

> **Proprietary software** — All rights reserved.

This repository is released under a custom proprietary license. No part of
this codebase, trained models, curated datasets, or derived works may be
reproduced, distributed, or used commercially without explicit written
permission from the author.

See [`LICENSE`](LICENSE) for the full terms.

For partnership, pilot access, or licensing inquiries:
- See [`docs/IP_AND_LICENSING.md`](docs/IP_AND_LICENSING.md) for IP strategy
  and commercial deal structures.
- Contact: **beshoyarmanios@gmail.com**
- `docs/architecture.md` (if present)
- `docs/labeling_scheme_v2.csv`
- `docs/curation_priority_v2.csv` 