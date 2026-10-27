# OligoTox

[![CI](https://github.com/barmanios123/OligoTox/actions/workflows/ci.yml/badge.svg)](https://github.com/barmanios123/OligoTox/actions/workflows/ci.yml)

OligoTox scores antisense oligos and siRNAs for 13 preclinical toxicity endpoints from design-time features only: sequence, chemistry and modality. It wraps the models in a Streamlit console that adds a GalNAc-siRNA PK/PD simulation and a regulatory-style safety report.

The main result is in the evaluation. The first version scored PR-AUC 0.93 on hepatocellular injury. Most of that came from label leakage and from splitting replicate rows of the same oligo across train and test. After fixing both, the honest number is 0.29 against a 0.13 base rate. This README reports the fixed numbers.

## Results

![Per-endpoint PR-AUC before and after the leakage fix](reports/results/pr_auc_by_endpoint.png)

Each bar is an out-of-fold PR-AUC from a calibrated logistic regression. Blue is the number to trust. In that setup, all rows of a sequence (RNA and DNA spellings count as the same sequence) stay on the same side of every split. That matches the real use case: scoring an oligo the model has never seen.

| Endpoint | Positive rows / sequences | Base rate | Before fix | Row-level CV | **Sequence-grouped CV** |
|---|---|---|---|---|---|
| kupffer cell activation | 34 / 27 | 0.021 | 0.93 | 0.95 | **0.89** |
| injection site reaction | 12 / 12 | 0.007 | 0.42 | 0.41 | **0.40** |
| thrombocytopenia | 4 / 4 | 0.002 | 0.15 | 0.39 | **0.36** |
| hepatocellular injury | 203 / 19 | 0.126 | 0.93 | 0.57 | **0.29** |
| coagulation effect | 5 / 5 | 0.003 | 0.32 | 0.35 | **0.28** |
| hepatic vacuolation | 12 / 6 | 0.007 | 0.57 | 0.57 | **0.23** |
| renal accumulation | 7 / 7 | 0.004 | 0.40 | 0.30 | **0.16** |
| offtarget gene silencing | 68 / 12 | 0.042 | 0.64 | 0.35 | **0.09** |
| renal tubular injury | 7 / 7 | 0.004 | 0.20 | 0.37 | **0.05** |
| innate immune activation | 7 / 7 | 0.004 | 0.09 | 0.34 | **0.02** |
| cytokine release | 7 / 7 | 0.004 | 0.09 | 0.34 | **0.02** |
| complement activation | 6 / 6 | 0.004 | 0.00 | 0.00 | **0.00** |

CNS toxicity has only 2 positives and is not evaluated. 1,614 rows cover 74 unique sequences.

### What the fix changed

- **Removed outcome columns from the features.** 23 `label_*` / `endpoint_*` columns were curated annotations of the outcome; `label_kupffer_activation` correlated 0.98 with its own target. Two assay readouts (`cas_avg`, `kd_avg`) were also removed, because they are measured after dosing and don't exist for a new design. 111 constant columns were dropped. 55 features remain.
- **Grouped CV by sequence.** Many rows are replicate measurements of the same oligo, so row-level folds leak. Folds now come from `StratifiedGroupKFold` on the normalized sequence.
- **Reported with the base rate.** For endpoints with fewer than 10 positive sequences, PR-AUC is noisy. Read the gap to the base rate, not the absolute value.

### How to read it

- **Kupffer cell activation** is the one endpoint with a strong, sequence-held-out signal. Its top predictors are design features (low GC, PS fraction, length). Those positives also come from a small number of curated sources, so part of the signal may be source rather than biology.
- **Hepatocellular injury** keeps a real but modest signal (2.3× the base rate).
- **Endpoints with 4–7 positive sequences** are not reliable estimates. Treat them as placeholders until more labeled data arrives.

## Run it

```bash
git clone https://github.com/barmanios123/OligoTox.git
cd OligoTox
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt && pip install -e . --no-deps
streamlit run oligotox_app.py
```

The trained models and their config ship in `artifacts/`, so the app scores candidates from a fresh clone without the training data.

### Tests

```bash
pip install pytest
pytest -q
```

The suite checks that leaky columns can't re-enter training, that published metrics are sequence-grouped, that app-time features match the trained columns, and that the app completes a full run headlessly. CI runs it on every push.

### Retrain and regenerate the figure

Retraining needs the processed training data in `data/processed/`, which is not in the repo.

```bash
python -m oligotox.training.train
cp artifacts/metrics/training_summary.json reports/results/training_summary_after.json
python scripts/make_results_figure.py
```

## What's in the app

| Tab | What it does | Status |
|---|---|---|
| Toxicity | 13 per-endpoint probabilities, a severity-weighted composite and a risk band | Trained models, validated above |
| PK/PD | Liver AUC, Cmax and half-life from a GalNAc-siRNA compartment simulator | Mechanistic simulator, not fitted to this data |
| Regulatory report | OSWG-style JSON safety summary built from the other tabs | Report formatting |
| Efficacy | OligoFormer, ENsiRNA, siDPT and other model slots | **Untrained stubs.** The app labels them; their outputs are not predictions |

## Limitations

- 74 unique sequences is small. Nine of the 13 endpoints have 7 or fewer positive sequences.
- Some training features (`ps_x_len`, `lna_x_gc`) were zeroed by upstream processing in about half the source rows. The app computes them consistently, so scores for real designs may shift slightly from CV behaviour.
- The structure features use a placeholder MFE, and the sequence embedding is a hash, not a learned model.
- The efficacy tab is not implemented beyond model interfaces.

## Repository layout

```
oligotox_app.py          Streamlit console
src/oligotox/            package: features, training, scoring, PK/PD, reporting
artifacts/models|config  shipped trained models and scoring config
reports/results/         metrics and the results figure
tests/                   pytest suite (run in CI)
benchmarks/oligogym/     leakage-aware benchmark scripts
docs/                    architecture and provenance notes
```

## License

Proprietary; see `LICENSE`. Available for evaluation and collaboration discussions.
