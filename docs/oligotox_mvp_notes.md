# OligoTox Notes

## Purpose

OligoTox is a Streamlit-based decision-support interface for RNA and oligonucleotide therapeutic safety triage.

## Current screens

1. Candidate Profiler
2. Candidate Comparator
3. Assay Planner
4. Off-Target Evaluator
5. Translational Bridge

## Backed by current repository artifacts

### Real repository-backed components

- Packaged Python project under `src/oligotox/`
- Editable install via `pyproject.toml`
- Training entrypoint at `src/oligotox/training/train.py`
- Inference entrypoint at `src/oligotox/inference/api.py`
- Trained artifacts written under `artifacts/`
- Processed matrices written under `data/processed/`
- Streamlit application app entrypoint at `oligotox_app.py`
- Assay planner rule table at `app/assay_rules.csv`

### platform scaffolded components

- Candidate scoring may fall back to heuristic logic when a callable prediction API is not exposed cleanly from `src/oligotox/inference/api.py`
- Off-target evaluator currently uses placeholder seed-match rows
- Translational bridge currently uses placeholder rule summaries
- Safety brief export is lightweight HTML rather than final PDF rendering

## Intended next upgrades

1. Replace fallback scoring with direct model artifact inference
2. Replace placeholder off-target logic with transcriptome or seed-search backend
3. Replace translational bridge heuristics with explicit species and biomarker rules
4. Add PDF export for candidate safety brief
5. Add procurement-ready assay workbook output with richer formatting

## Recommended collaborator framing

The platform should be presented as a working product shell with real repository structure and training outputs, plus several decision-support layers that are intentionally scaffolded for rapid iteration.
