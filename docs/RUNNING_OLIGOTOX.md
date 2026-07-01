# Running OligoTox Locally

This guide explains how to clone, set up, validate, and run the OligoTox platform from the GitHub repository.

## Prerequisites

Before starting, make sure the local machine has:

- Python 3.10 or newer
- Git installed
- A terminal available
- Internet access for cloning the repository and installing dependencies

Recommended:

- A virtual environment created with `venv` or Conda
- Streamlit installed in the same Python environment used for the project
- Basic familiarity with terminal commands

## 1. Clone the repository

```bash
git clone https://github.com/barmanios123/OligoTox.git
cd OligoTox
```

## 2. Create and activate a Python environment

### Option A: venv

```bash
python -m venv .venv
source .venv/bin/activate
```

### Option B: Conda

```bash
conda create -n oligotox python=3.11 -y
conda activate oligotox
```

Use one environment only.

## 3. Install dependencies

```bash
pip install -U pip
pip install -e .
pip install streamlit pandas numpy scikit-learn joblib
```

## 4. Confirm the repository structure

Key paths:

- `oligotox_app.py`
- `validate_platform.py`
- `src/`
- `scripts/`
- `docs/`
- `data/`
- `artifacts/`

Quick check:

```bash
ls
```

## 5. Validate the platform

Run:

```bash
python validate_platform.py
```

A successful run should show passing checks for imports, data/model artifacts, efficacy, PK/PD, reporting, toxicity, and app presence.

## 6. Run the project health check

Run:

```bash
bash scripts/project_health.sh
```

This verifies repository layout, imports, training inputs, model artifacts, and a basic inference smoke test.

## 7. Launch the app

Preferred command:

```bash
bash scripts/run_app.sh
```

Equivalent direct command:

```bash
streamlit run oligotox_app.py --server.port 8502
```

Then open:

```text
http://localhost:8502
```

If port `8502` is already in use, stop the previous process or use another port.

## 8. Use the app

1. Enter an oligonucleotide sequence.
2. Choose the modality.
3. Choose the chemistry.
4. Enter dose information if PK/PD context is needed.
5. Run the analysis.
6. Review toxicity, efficacy, PK/PD, regulatory report, saved runs, and build info.

## 9. Retrain models (optional)

If retraining is needed:

```bash
python -m oligotox.training.train
```

Then rerun:

```bash
python validate_platform.py
bash scripts/project_health.sh
```

## 10. Common issues

### Port already in use

```bash
lsof -tiTCP:8502 -sTCP:LISTEN | xargs kill -9
```

### Wrong environment

```bash
which python
python -V
```

### Missing package

```bash
pip install <package_name>
```

### Missing artifacts

If validation reports missing models or configs, ensure the required `artifacts/` and `data/` content is present locally, or regenerate it through training if appropriate.

## 11. Recommended startup sequence

```bash
cd OligoTox
source .venv/bin/activate
python validate_platform.py
bash scripts/project_health.sh
bash scripts/run_app.sh
```

Then open:

```text
http://localhost:8502
```
