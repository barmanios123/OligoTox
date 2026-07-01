# OligoTox

OligoTox is a mechanism-aware safety, translational interpretation, and candidate triage platform for oligonucleotide therapeutics. It combines multi-endpoint toxicity scoring, efficacy model integration, PK/PD context, and regulatory-style reporting in a Streamlit-based operator workflow.

## Why this repository exists

Oligonucleotide programs can fail for very different reasons: exaggerated on-target pharmacology, hybridization-dependent off-target effects, and liabilities introduced by chemistry, formulation, conjugates, or delivery strategy. OligoTox is organized around those risk classes to support earlier, more structured candidate review.

The repository is designed to present a clean package + app + docs story:

- `src/oligotox/` contains the core package code
- `oligotox_app.py` is the primary Streamlit operator console entrypoint
- `validate_platform.py` is the primary end-to-end validation command
- `docs/` contains setup, architecture, and supporting documentation
- `scripts/` contains health checks and launch helpers

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

### Core runtime

- `oligotox_app.py` — primary Streamlit operator console entrypoint
- `validate_platform.py` — primary end-to-end validation entrypoint
- `src/oligotox/` — core package code for inference, training, efficacy, PK/PD, and reporting

### Project support

- `scripts/` — launchers, health checks, and utilities
- `docs/` — architecture, setup, notes, and integration guidance
- `configs/` — configuration and integration files
- `benchmarks/` — benchmarking utilities
- `reports/` — generated and development-stage reports

### Data and artifacts

- `data/` — local data assets used by the workflow
- `artifacts/` — local runtime artifacts such as model outputs, configs, and run products

Some local-only, generated, archival, or environment-specific content may exist in the working tree without being intended as part of the clean public repo surface. The public repository should be interpreted primarily through the package code, app entrypoint, validation entrypoint, scripts, and docs.

## One obvious app entrypoint

Launch the app with:

```bash
bash scripts/run_app.sh
```

This wrapper runs the canonical Streamlit entrypoint:

```bash
streamlit run oligotox_app.py --server.port 8502
```

Then open:

```text
http://localhost:8502
```

## One obvious validation command

Validate the platform with:

```bash
python validate_platform.py
```

For a broader operational check, run:

```bash
bash scripts/project_health.sh
```

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
bash scripts/run_app.sh
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

Start here:

- `docs/RUNNING_OLIGOTOX.md` — step-by-step local setup and execution guide
- `docs/REPOSITORY_CONVENTIONS.md` — what is source code vs local/generated material

## Status

The repository contains a working Streamlit app entrypoint, validated backend components, health-check scripts, and a cleaner public-facing structure. The platform is suitable for internal review, demo workflows, and iterative client-facing refinement.

## License

**Proprietary and Confidential**

This software is proprietary and confidential. It is provided for evaluation purposes only to authorized clients and collaborators.

- No commercial use without license
- No modification or derivative works
- No distribution or sublicensing
- Evaluation and collaboration discussions only

For licensing inquiries, contact the repository owner. See the `LICENSE` file for full terms.
