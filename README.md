# OligoTox

[![CI](https://github.com/barmanios123/OligoTox/actions/workflows/ci.yml/badge.svg)](https://github.com/barmanios123/OligoTox/actions/workflows/ci.yml)

OligoTox scores antisense oligos and siRNAs for 13 preclinical toxicity endpoints from design-time features only: sequence, chemistry and modality. It wraps the models in a Streamlit console that adds a GalNAc-siRNA PK/PD simulation and a regulatory-style safety report.

The main result is in the evaluation. The first version scored PR-AUC 0.93 on hepatocellular injury. Most of that came from label leakage and from splitting replicate rows of the same oligo across train and test. After fixing both, the honest number is 0.27 against a 0.13 base rate. Baselines then showed that the one apparently strong endpoint is mostly dataset recognition, not biology. Of the 13 endpoints, 1 survives that scrutiny, and the app presents only that one as a prediction — the rest run but are labelled unsupported, with the reason attached.

## Results

![Per-endpoint PR-AUC before and after the leakage fix](reports/results/pr_auc_by_endpoint.png)

Each bar is an out-of-fold PR-AUC from a calibrated logistic regression. Blue is the number to trust. In that setup, all rows of a sequence (RNA and DNA spellings count as the same sequence) stay on the same side of every split. That matches the real use case: scoring an oligo the model has never seen.

| Endpoint | Positive rows / sequences | Base rate | Before fix | Row-level CV | **Sequence-grouped CV** |
|---|---|---|---|---|---|
| kupffer cell activation | 34 / 27 | 0.021 | 0.93 | 0.75 | **0.48** |
| hepatocellular injury | 203 / 19 | 0.126 | 0.93 | 0.57 | **0.27** |
| injection site reaction | 12 / 12 | 0.007 | 0.42 | 0.36 | **0.26** |
| offtarget gene silencing | 68 / 12 | 0.042 | 0.64 | 0.32 | **0.19** |
| hepatic vacuolation | 12 / 6 | 0.007 | 0.57 | 0.58 | **0.16** |
| complement activation | 6 / 6 | 0.004 | 0.00 | 0.24 | **0.09** |
| coagulation effect | 5 / 5 | 0.003 | 0.32 | 0.37 | **0.06** |
| renal tubular injury | 7 / 7 | 0.004 | 0.20 | 0.01 | **0.06** |
| innate immune activation | 7 / 7 | 0.004 | 0.09 | 0.12 | **0.06** |
| cytokine release | 7 / 7 | 0.004 | 0.09 | 0.12 | **0.06** |
| thrombocytopenia | 4 / 4 | 0.002 | 0.15 | 0.25 | **0.03** |
| renal accumulation | 7 / 7 | 0.004 | 0.40 | 0.05 | **0.03** |

CNS toxicity has only 2 positives and is not evaluated. 1,614 rows cover 74 unique sequences.

### What the fix changed

- **Removed outcome columns from the features.** 23 `label_*` / `endpoint_*` columns were curated annotations of the outcome; `label_kupffer_activation` correlated 0.98 with its own target. Two assay readouts (`cas_avg`, `kd_avg`) were also removed, because they are measured after dosing and don't exist for a new design. 111 constant columns were dropped. 55 features remain.
- **Grouped CV by sequence.** Many rows are replicate measurements of the same oligo, so row-level folds leak. Folds now come from `StratifiedGroupKFold` on the normalized sequence.
- **Recomputed features from sequence.** `tlr7_score`, `ps_x_len`, `lna_x_gc`, `seed_gc` and `cpg_density` had been filled inconsistently across the source datasets, and were zero in about half the rows. They are now derived from the sequence for every training row, the same way the app derives them at scoring time. The shipped models were retrained on the recomputed features, so these numbers differ slightly from the first post-leakage run.
- **Reported with the base rate.** For endpoints with fewer than 10 positive sequences, PR-AUC is noisy. Read the gap to the base rate, not the absolute value.

## Baselines

![Model vs baselines with 95% bootstrap CIs](reports/results/baselines.png)

A PR-AUC above the base rate is not evidence that the model learned design chemistry. It could be learning which dataset a row came from, or one crude number like GC content. `scripts/evaluate_baselines.py` scores four baselines against the shipped model on the same sequence-grouped folds, with 95% bootstrap CIs resampled over sequences. The **data source only** baseline sees nothing but which dataset a row came from.

| Endpoint | Base rate | Logistic (55 features) | Gradient boosting | GC only | Length only | Data source only |
|---|---|---|---|---|---|---|
| kupffer cell activation | 0.021 | 0.48 [0.32–0.70] | 0.42 [0.29–0.60] | 0.07 | 0.14 | **0.69 [0.53–0.87]** |
| hepatocellular injury | 0.126 | 0.27 [0.19–0.52] | 0.37 [0.21–0.56] | 0.15 | 0.11 | 0.20 [0.10–0.35] |
| injection site reaction | 0.007 | 0.26 [0.06–0.49] | 0.17 | 0.03 | 0.20 | 0.01 |
| offtarget gene silencing | 0.042 | 0.19 [0.05–0.59] | 0.20 | 0.06 | 0.04 | 0.06 |
| hepatic vacuolation | 0.007 | 0.16 [0.04–0.34] | 0.24 | 0.01 | 0.02 | 0.20 [0.06–0.42] |

All 12 evaluated endpoints are in [`reports/results/baselines.csv`](reports/results/baselines.csv); the within-source comparison is in [`reports/results/within_source.csv`](reports/results/within_source.csv).

### How to read it

- **Why a source indicator predicts anything at all.** 1,572 of the 1,614 rows come from two systematic sources (`oligotox_main`, `papargyri_2020`); the remaining 42 are single rows hand-curated from individual papers and FDA labels. Those 42 exist *because* someone published a toxicity finding about them, so 36 of 42 carry a positive label, against 214 of 1,572 in the systematic sources. "Which source is this row from" is therefore close to the label itself, by construction. That is the selection bias the source-only baseline measures.
- **Kupffer cell activation is a dataset artifact, not a finding.** It has the highest grouped PR-AUC (0.48), but the data-source-only baseline beats it (0.69) and the two CIs overlap. All 34 positive rows (27 sequences) sit in that 42-row curated bucket and none in the two systematic sources, and its base rate there is 0.81. Scored inside that source alone the model gets 0.79, below its own base rate, and GC content alone scores higher (0.87). The model is recognising which dataset a row came from. Do not read this endpoint as a toxicity signal.
- **Hepatocellular injury is the most defensible endpoint.** The 55-feature model (0.27) and gradient boosting (0.37) both clear the data-source baseline (0.20) and the GC-only and length-only baselines, and the signal survives within a single source: inside `papargyri_2020` alone the model gets 0.34 against a 0.25 base rate. It is a modest signal, roughly 2.2× the base rate, on 19 positive sequences.
- **Injection site reaction and off-target gene silencing** clear their source-only baselines but have wide CIs reaching near zero. Treat them as weak and unconfirmed.
- **Endpoints with 4–7 positive sequences** are not reliable estimates, and several are beaten by the length-only baseline. Treat them as placeholders until more labeled data arrives.
- **Gradient boosting matches or beats the shipped logistic model** on several endpoints, but every CI overlaps. With 74 sequences there is no evidence to prefer either, so the shipped model stays the simpler, calibrated one.

## What the app will and won't tell you

The evaluation above is enforced in the app, not just described here. `scripts/derive_evidence_tiers.py` turns `baselines.csv` into `artifacts/config/endpoint_evidence.json`, which the scoring pipeline and the report read at runtime. Each endpoint lands in one tier:

| Tier | Endpoints | How the app presents it |
|---|---|---|
| **Validated** | hepatocellular injury | A prediction. Enters the composite score. |
| **Weak signal** | injection site reaction, offtarget gene silencing | Model output shown with its CI and the baseline it fails to clear. Not a prediction. |
| **Source-confounded** | kupffer cell activation | Model output shown, labelled as dataset recognition. Not a prediction. |
| **Insufficient data** | the remaining 9 | Model output shown with its positive-sequence count. Not a prediction. |

The practical consequences:

- **The composite score averages one endpoint, not 13,** and the tab says so ("Composite covers 1 of 13 endpoints"). Previously it was a severity-weighted mean over all 13, so eleven unsupported models diluted the one that works. On a test candidate the old composite named *renal tubular injury* — 7 positive sequences, beaten by a length-only baseline — as the dominant endpoint.
- **Unsupported endpoints never return `predicted_positive`.** They return `null` with an `evidence_tier` and a plain-English `evidence_note` saying which baseline beat them.
- **The regulatory report raises no recommendation from an unsupported endpoint,** whatever it scores. It used to emit "Immunogenicity assessment required per OSWG guidelines" from a complement-activation model fit on 6 positive sequences. A report with no validated endpoint is classified `UNCLASSIFIED_NO_VALIDATED_ENDPOINT` rather than `LOW_RISK`.

Nothing is hidden: every model output is still shown and still in the downloadable JSON. The tier decides whether it is labelled a prediction.

## Run it

```bash
git clone https://github.com/barmanios123/OligoTox.git
cd OligoTox
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt && pip install -e . --no-deps
streamlit run oligotox_app.py
```

The trained models and their config ship in `artifacts/`, so the app scores candidates from a fresh clone without the training data. The app, the tests and CI need no deep learning stack; `torch` is an optional extra (`pip install -e ".[efficacy]"`) used only by unshipped experiments.

### Tests

```bash
pip install pytest
pytest -q
```

The suite checks that leaky columns can't re-enter training, that published metrics are sequence-grouped, that app-time features match the trained columns, that the baseline results cover every evaluated endpoint with valid CIs, and that the app completes a full run headlessly. It also re-derives every evidence tier from `baselines.csv` and asserts the shipped config agrees, so the tiers can't drift from the evaluation; that the predictions table lists only validated endpoints; and that an unsupported endpoint scoring 1.0 raises no regulatory recommendation. CI runs it on every push.

## Reproducing the numbers

The data every published number is computed from is in the repo, so the evaluation is checkable rather than claimed. Four files, 1.2 MB total: `feature_matrix.csv` (1,614 x 55), `label_matrix.csv`, `training_dataset.csv` (sequences, source and curated labels) and `dataset_metadata.json`. Every source in them is public literature, a PMID/PMC paper or an FDA label; the other contents of `data/processed/` are intermediates and stay out.

From a fresh clone, training takes about 10 seconds and the baselines about 95 seconds:

```bash
python -m oligotox.training.train       # retrains all 13 endpoint models
python scripts/evaluate_baselines.py    # recomputes baselines and bootstrap CIs
```

Verified in a clean clone of this commit: all 13 per-endpoint PR-AUCs reproduce to the published figure exactly, `baselines.csv` and `within_source.csv` reproduce with zero difference across all 60 and 36 rows, `endpoint_evidence.json` comes out byte-identical, and so do all 13 shipped `.joblib` models. The models in `artifacts/` are the ones this data and this script produce.

### Retrain and regenerate the figures

```bash
python -m oligotox.training.train
cp artifacts/metrics/training_summary.json reports/results/training_summary_after.json
python scripts/make_results_figure.py
python scripts/evaluate_baselines.py     # writes baselines.csv and within_source.csv
python scripts/make_baselines_figure.py  # writes baselines.png
python scripts/derive_evidence_tiers.py  # writes artifacts/config/endpoint_evidence.json
```

The last step is what the app reads, so rerun it after any change to the baselines, or the tier test will fail.

## What's in the app

| Tab | What it does | Status |
|---|---|---|
| Candidate Input | Sequence, chemistry, modality and dose entry, with validation | — |
| Toxicity Results | Predictions for the validated endpoints, a composite over those only, and the remaining model outputs behind an expander with the reason each is excluded | Gated by evidence tier |
| PK/PD | Liver AUC, Cmax and half-life from a GalNAc-siRNA compartment simulator | Mechanistic simulator, not fitted to this data |
| Regulatory Report | OSWG-style JSON safety summary, with an evidence block and recommendations raised only from validated endpoints | Report formatting |
| Saved Runs | Runs kept during the session | — |
| Build Info | Model versions, feature count and config provenance | — |

An earlier Efficacy tab held untrained model slots and produced no predictions. It was removed rather than shipped behind a disclaimer.

## Limitations

- 74 unique sequences is small. Nine of the 13 endpoints have 7 or fewer positive sequences.
- **Only 1 of the 13 endpoints is presented as a prediction.** Hepatocellular injury is the one whose signal survives both the data-source baseline and a within-source check; Kupffer cell activation, the highest raw score, does not. The other 12 models ship and run, but the app labels them as unsupported.
- The one validated endpoint rests on 19 positive sequences, and the within-source confirmation comes from a single dataset. A second independent source is what would settle it.
- Several endpoints draw all their positives from one source dataset, so source and biology cannot be separated for them with this data.
- The structure features use a placeholder MFE, and the sequence embedding is a hash, not a learned model.
- The PK/PD simulator is mechanistic and has not been fitted or validated against measured exposure data.

## Repository layout

```
oligotox_app.py          Streamlit console
src/oligotox/            package: features, training, scoring, PK/PD, reporting
artifacts/models|config  shipped trained models, scoring config, evidence tiers
reports/results/         metrics, baselines and the figures
scripts/                 results figure, baseline evaluation, tier derivation
tests/                   pytest suite (run in CI)
benchmarks/oligogym/     leakage-aware benchmark scripts
docs/                    architecture and provenance notes
```

## License

Proprietary; see `LICENSE`. Available for evaluation and collaboration discussions.
