# OligoTox

OligoTox is a mechanism-aware safety, translational interpretation, and candidate triage platform for oligonucleotide therapeutics. It combines multi-endpoint toxicity scoring, efficacy model integration, PK/PD context, and regulatory-style reporting in a Streamlit-based operator workflow.

## Overview

OligoTox is designed to evaluate oligonucleotide candidates across the major preclinical liability classes relevant to translational decision-making. Its scientific framing is organized around three broad mechanistic domains:

- **D1:** exaggerated on-target pharmacology
- **D2:** hybridization-dependent off-target effects
- **D3:** chemistry, formulation, conjugate, and delivery-related liabilities

The platform is intended to support candidate review for modalities such as ASOs, gapmers, siRNAs, and related oligonucleotide formats.

## Current capabilities

The current platform supports the following verified capabilities:

- Streamlit app launch through `oligotox_app.py`
- Multi-endpoint toxicity inference through the OligoTox backend
- Modality- and chemistry-aware candidate input
- Endpoint-level output including composite score, dominant endpoint, and per-endpoint detail
- Efficacy model integration through the package registry
- PK/PD context integration in the operator console workflow
- Regulatory-style safety reporting
- JSON-style result export
- Local validation through `validate_platform.py`
- Project health checks through `scripts/project_health.sh`

## Repository structure

Key areas in the repository include:

- `oligotox_app.py` — primary Streamlit operator console entrypoint
- `src/oligotox/` — core package code for inference, training, efficacy, PK/PD, and reporting
- `scripts/` — launchers, health checks, and utility scripts
- `docs/` — architecture, notes, integrations, and supporting documentation
- `configs/` — configuration and integration files
- `benchmarks/` — benchmarking utilities
- `data/` — local data assets used in the workflow
- `reports/` — generated and development-stage reports

## Quickstart

### 1. Clone the repository

```bash
git clone https://github.com/barmanios123/OligoTox.git
cd OligoTox
```

### 2. Create and activate an environment

Using `venv`:

```bash
python -m venv .venv
source .venv/bin/activate
```

Or using Conda:

```bash
conda create -n oligotox python=3.11 -y
conda activate oligotox
```

### 3. Install dependencies

```bash
pip install -U pip
pip install -e .
pip install streamlit pandas numpy scikit-learn joblib
```

### 4. Validate the platform

```bash
python validate_platform.py
bash scripts/project_health.sh
```

### 5. Launch the app

```bash
streamlit run oligotox_app.py --server.port 8502
```

Open:

```text
http://localhost:8502
```

## Typical usage flow

1. Launch the app locally.
2. Enter a candidate oligonucleotide sequence.
3. Select modality and chemistry.
4. Optionally provide dose and PK/PD context.
5. Run the analysis.
6. Review toxicity, efficacy, PK/PD, and regulatory-style outputs.
7. Save or export the resulting run outputs.

## Training

If retraining is needed and the local environment includes the training module:

```bash
python -m oligotox.training.train
```

After retraining, rerun validation and health checks:

```bash
python validate_platform.py
bash scripts/project_health.sh
```

## Documentation

For detailed setup and operating instructions, see:

- `docs/RUNNING_OLIGOTOX.md`

## Status

The repository currently contains a working Streamlit app entrypoint, validated backend components, health-check scripts, and a cleaned public-facing repository structure. The platform is suitable for internal review, demo workflows, and iterative client-facing refinement.

## License

**Proprietary and Confidential**

This software is proprietary and confidential. It is provided for evaluation purposes only to authorized clients and collaborators.

- No commercial use without license
- No modification or derivative works
- No distribution or sublicensing
- Evaluation and collaboration discussions only

For licensing inquiries, contact the repository owner. See the `LICENSE` file for full terms.