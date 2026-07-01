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

Open a terminal and run:

```bash
git clone https://github.com/barmanios123/OligoTox.git
cd OligoTox
```

This downloads the repository and moves into the project folder.

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

Use one environment only. Do not mix multiple Python environments while setting up the platform.

## 3. Install dependencies

From the repository root, run:

```bash
pip install -U pip
pip install -e .
```

If additional packages are needed, install them in the same environment:

```bash
pip install streamlit pandas numpy scikit-learn joblib
```

If the repository later adds more formal dependency definitions, prefer those over manual installs.

## 4. Confirm the repository structure

From the repository root, confirm that the key files and folders exist:

- `oligotox_app.py`
- `validate_platform.py`
- `pyproject.toml`
- `src/`
- `scripts/`
- `docs/`
- `data/`
- `artifacts/`

Quick check:

```bash
ls
```

## 5. Run the platform validation

Before launching the app, validate that the backend is working:

```bash
python validate_platform.py
```

A successful run should show passing checks for:

- module imports
- data loading and model artifacts
- efficacy models
- PK/PD simulation
- FDA/report generation
- toxicity prediction
- Streamlit app presence

If validation fails, resolve the reported environment, artifact, or dependency issue before continuing.

## 6. Run the project health check

Run:

```bash
bash scripts/project_health.sh
```

This checks repository layout, Python imports, training inputs, model artifacts, and a basic inference smoke test.

If this passes, the local platform is in good shape for app launch.

## 7. Launch the Streamlit app

Run:

```bash
streamlit run oligotox_app.py --server.port 8502
```

Then open:

```text
http://localhost:8502
```

If port `8502` is already in use, either stop the previous process or use a different port:

```bash
streamlit run oligotox_app.py --server.port 8503
```

## 8. Use the app

Once the app opens:

1. Enter an oligonucleotide sequence.
2. Choose the modality.
3. Choose the chemistry.
4. Enter dose information if PK/PD context is needed.
5. Run the analysis.
6. Review:
   - candidate input summary
   - toxicity results
   - efficacy results
   - PK/PD output
   - regulatory report
   - saved runs
   - build information

The app should save run outputs locally and display endpoint-level results.

## 9. Run training (optional)

If retraining is needed and the local environment includes the training module, run:

```bash
python -m oligotox.training.train
```

After retraining, rerun:

```bash
python validate_platform.py
bash scripts/project_health.sh
```

This confirms the retrained state is still working.

## 10. Run tests (optional)

If tests are configured in the environment, run:

```bash
pytest
```

Depending on the local repo state, this may be limited to currently configured tests.

## 11. Common issues

### Port already in use

If Streamlit reports that the port is already in use:

```bash
lsof -tiTCP:8502 -sTCP:LISTEN | xargs kill -9
```

Then restart the app.

### Missing package

If an import fails, install the missing package into the active environment:

```bash
pip install <package_name>
```

### Validation fails because of missing artifacts

If validation reports missing models or config files, make sure the local repository includes the expected `artifacts/` and `data/` content. If necessary, rerun training to regenerate local outputs.

### Wrong environment

If behavior is inconsistent, confirm which Python is active:

```bash
which python
python -V
```

Make sure Streamlit, validation, training, and package imports all use the same environment.

## 12. Recommended startup sequence

For a clean local startup, use:

```bash
cd OligoTox
source .venv/bin/activate
python validate_platform.py
bash scripts/project_health.sh
streamlit run oligotox_app.py --server.port 8502
```

Then open:

```text
http://localhost:8502
```

## 13. Files most users should know

- `oligotox_app.py` — main app entrypoint
- `validate_platform.py` — end-to-end validation
- `scripts/project_health.sh` — platform health check
- `src/oligotox/` — core implementation
- `docs/` — documentation
- `README.md` — overview and quickstart

## 14. Demo checklist

Before showing the platform to a client or collaborator:

1. Activate the correct environment.
2. Run platform validation.
3. Run the health check.
4. Launch the Streamlit app.
5. Test one known example sequence.
6. Keep the working browser tab open before the meeting.

This helps reduce environment drift and last-minute demo issues.