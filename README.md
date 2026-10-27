# OligoTox

[![CI](https://github.com/barmanios123/OligoTox/actions/workflows/ci.yml/badge.svg)](https://github.com/barmanios123/OligoTox/actions/workflows/ci.yml)

OligoTox scores antisense oligos and siRNAs for 13 preclinical toxicity endpoints from design-time features only: sequence, chemistry and modality. It wraps the models in a Streamlit console that adds a GalNAc-siRNA PK/PD simulation and a regulatory-style safety report.

The main result is in the evaluation. The first version scored PR-AUC 0.93 on hepatocellular injury. Most of that came from label leakage and from splitting replicate rows of the same oligo across train and test. After fixing both, the honest number is 0.27 against a 0.13 base rate. Baselines then showed that the one apparently strong endpoint is mostly dataset recognition, not biology. Of the 13 endpoints, 2 now survive that scrutiny, one of them validated on an independent external assay. The app presents only those two as predictions; the rest run but are labelled unsupported, with the reason attached.

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
| **Validated** | hepatocellular injury, CNS neurotoxicity | A prediction. Enters the composite score. |
| **Weak signal** | injection site reaction, offtarget gene silencing | Model output shown with its CI and the baseline it fails to clear. Not a prediction. |
| **Source-confounded** | kupffer cell activation | Model output shown, labelled as dataset recognition. Not a prediction. |
| **Insufficient data** | the remaining 8 | Model output shown with its positive-sequence count. Not a prediction. |

The practical consequences:

- **The composite score averages one endpoint, not 13,** and the tab says so ("Composite covers 1 of 13 endpoints"). Previously it was a severity-weighted mean over all 13, so eleven unsupported models diluted the one that works. On a test candidate the old composite named *renal tubular injury* — 7 positive sequences, beaten by a length-only baseline — as the dominant endpoint.
- **Unsupported endpoints never return `predicted_positive`.** They return `null` with an `evidence_tier` and a plain-English `evidence_note` saying which baseline beat them.
- **The regulatory report raises no recommendation from an unsupported endpoint,** whatever it scores. It used to emit "Immunogenicity assessment required per OSWG guidelines" from a complement-activation model fit on 6 positive sequences. A report with no validated endpoint is classified `UNCLASSIFIED_NO_VALIDATED_ENDPOINT` rather than `LOW_RISK`.

Nothing is hidden: every model output is still shown and still in the downloadable JSON. The tier decides whether it is labelled a prediction.

## CNS neurotoxicity, validated on an independent assay

Every other number here is internal: one curation, cross-validated. CNS neurotoxicity is different, and it is the strongest result in the project.

Two public datasets measure it on sequences this repo never trained on, and they share no sequences with each other:

| Dataset | Sequences | Assay | Readout |
|---|---|---|---|
| `MOE_Neurotox_2024` | 1,685 | in vivo | FOB score at 3h |
| `Hagedorn_2022` | 1,653 | in vitro | calcium oscillation, LNA gapmers |

Training on one and testing on the other is a real external test, across different chemistry, different assays and different in vivo/in vitro settings.

| Direction | All sequence features | G-count rule | GC only |
|---|---|---|---|
| MOE to Hagedorn | 0.828 [0.808-0.847] | **0.839** | 0.717 |
| Hagedorn to MOE | 0.783 [0.761-0.803] | **0.797** | 0.628 |

ROC-AUC with 95% bootstrap CIs. Both directions clear chance comfortably, so the signal is real and it transfers.

**The honest caveat: the model does not beat the published rule.** Hagedorn et al. report that G count and the distance from the 3' end to the nearest G drive acute neurotoxicity. Those four features alone score the same or slightly better than the full feature set in both directions. So this endpoint is best read as a reimplementation of their rule, validated on an independent assay, rather than a new finding. The evidence config says so in the text the app displays.

The shipped `cns_toxicity` model is trained on all 3,338 sequences and uses sequence-derived features only, since the external data carries no chemistry annotation, so it has its own column list in `artifacts/config/endpoint_features.json`. It replaces a model fit on **2 positive sequences** that training had marked `fallback_low_count_fit`.

### A label direction worth flagging

Hagedorn's readout is calcium oscillation, and toxicity shows up as a **reduction** in it, so a low score means toxic. `data/oligogym/oligogym_tox_combined_fixed.csv` has this inverted: it labels high scores as toxic and ends up calling 83% of the compounds neurotoxic. The threshold used here is stated in `scripts/evaluate_cns_external.py` and swept over five values, and the result holds across all of them (MOE 0.81 at every cut-off, Hagedorn 0.90 at the strictest).

### External data ingestion

`scripts/normalize_external_tox.py` is the single entry point for the datasets in `data/oligogym/`. It exists because the pre-existing combined files could not be checked: one has a `y_binary` column holding rounded raw readouts, and the other is binary but has Hagedorn's direction inverted. Neither records a threshold or a direction.

Every readout's direction, threshold and citation is declared in one table in that script, and a dataset missing from it is an error rather than a silent default. It also separates the three potency datasets (McQuisten, Ichihara x2, 6,718 sequences) from the three toxicity ones, because a target-inhibition readout sitting in a file named `tox_combined` is how a potency number ends up trained as a safety label.

```bash
python scripts/normalize_external_tox.py
```

Writes `external_tox_normalized.csv`, `external_efficacy_normalized.csv` and a manifest recording every labelling decision. Tests assert the direction of each source and that no efficacy readout carries a toxic label.

## Efficacy, back with evidence behind it

The Efficacy tab was removed in an earlier pass because it held untrained model slots whose outputs were not predictions. It is back, trained on 6,718 potency sequences that were sitting unused in `data/oligogym`.

Potency has its own shortcut. It varies far more between target genes than within one, so a model that has seen other oligos against the same gene can score well while knowing nothing about design. Two guards: folds are grouped by target gene, and a target-gene-mean baseline that ignores the sequence entirely is scored alongside.

| Model | Data | Gene-grouped CV | Cross-study | GC only | Gene mean | Tier |
|---|---|---|---|---|---|---|
| siRNA potency | 2,838 duplexes, 41 genes | **0.53** [0.50-0.56] | **0.44** / **0.48** | 0.30 | -0.11 | validated |
| ASO activity | 3,867 PS-ASOs, 85 genes | **0.25** [0.22-0.28] | none available | 0.14 | -0.15 | weak |

Spearman between measured potency and out-of-fold prediction. The siRNA number transfers between two independent studies (Ichihara 1 and 2) that share no sequences, in both directions. The gene-mean baseline carries none of it, and 17-29% of the variance is between genes, which is why the grouping matters.

The tab only scores the modality a model was trained for: ask it about a PMO and it says no model is validated for that, rather than returning a number. Scores are raw readouts in the units of the training data, labelled as such, never probabilities.

**Caveat on the siRNA model:** the training data is 2007 chemistry with no modification annotation, so it scores sequence design only and knows nothing about the chemistry of a modern GalNAc-siRNA.

### A duplex parsing bug this turned up

Ichihara's entries are siRNA duplexes written as one string, `sense.antisense`. Read naively they look like 41-nt "sequences" whose length and GC content are computed over both strands and the separator. The normalizer now splits them and uses the 21-nt antisense guide, which is what drives potency.

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
| Efficacy | Potency for the candidate's modality, with its units, evidence and caveat | siRNA validated, ASO weak |
| Saved Runs | Runs kept during the session | — |
| Build Info | Model versions, feature count and config provenance | — |

The Efficacy tab previously held untrained model slots and was removed rather than shipped behind a disclaimer. It returned once the models were trained and validated against baselines.

## Limitations

- 74 unique sequences is small. Nine of the 13 endpoints have 7 or fewer positive sequences.
- **Only 2 of the 13 endpoints are presented as predictions.** Hepatocellular injury survives both the data-source baseline and a within-source check; Kupffer cell activation, the highest raw score, does not. CNS neurotoxicity is validated on external data but reproduces a published rule rather than improving on it. The other 11 models ship and run, but the app labels them unsupported.
- Hepatocellular injury rests on 19 positive sequences with within-source confirmation from a single dataset. A second independent source is what would settle it, and no external hepatotoxicity data is on hand.
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
