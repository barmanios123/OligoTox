# OligoTox

[![CI](https://github.com/barmanios123/OligoTox/actions/workflows/ci.yml/badge.svg)](https://github.com/barmanios123/OligoTox/actions/workflows/ci.yml)

OligoTox scores antisense oligos and siRNAs for 13 preclinical toxicity endpoints from design-time features only: sequence, chemistry and modality. It wraps the models in a Streamlit console that adds a GalNAc-siRNA PK/PD simulation and a regulatory-style safety report.

The main result is in the evaluation. The first version scored PR-AUC 0.93 on hepatocellular injury. Most of that came from label leakage and from splitting replicate rows of the same oligo across train and test. After fixing both, the honest number is 0.27 against a 0.13 base rate. Baselines then showed that the one apparently strong endpoint is mostly dataset recognition, not biology. Of the 13 endpoints, 3 now survive that scrutiny, one of them validated on an independent external assay. The app presents only those two as predictions; the rest run but are labelled unsupported, with the reason attached.

## Results

![Per-endpoint PR-AUC before and after the leakage fix](reports/results/pr_auc_by_endpoint.png)

Each bar is an out-of-fold PR-AUC from a calibrated logistic regression. Blue is the number to trust. In that setup, all rows of a sequence (RNA and DNA spellings count as the same sequence) stay on the same side of every split. That matches the real use case: scoring an oligo the model has never seen.

| Endpoint | Positive rows / sequences | Base rate | Before fix | Row-level CV | **Sequence-grouped CV** |
|---|---|---|---|---|---|
| kupffer cell activation | 34 / 27 | 0.021 | 0.93 | 0.73 | **0.48** |
| hepatocellular injury | 203 / 19 | 0.126 | 0.93 | 0.55 | **0.29** |
| injection site reaction | 12 / 12 | 0.007 | 0.42 | 0.36 | **0.27** |
| offtarget gene silencing | 68 / 12 | 0.042 | 0.64 | 0.31 | **0.23** |
| hepatic vacuolation | 12 / 6 | 0.007 | 0.57 | 0.44 | **0.16** |
| complement activation | 6 / 6 | 0.004 | 0.00 | 0.22 | **0.07** |
| coagulation effect | 5 / 5 | 0.003 | 0.32 | 0.36 | **0.07** |
| renal tubular injury | 7 / 7 | 0.004 | 0.20 | 0.01 | **0.06** |
| cytokine release | 7 / 7 | 0.004 | 0.09 | 0.12 | **0.06** |
| innate immune activation | 7 / 7 | 0.004 | 0.09 | 0.12 | **0.06** |
| thrombocytopenia | 4 / 4 | 0.003 | 0.15 | 0.25 | **0.02** |
| renal accumulation | 7 / 7 | 0.004 | 0.40 | 0.05 | **0.02** |

CNS toxicity has only 2 positives in this data and is not evaluated here; it is trained and validated separately on external assays, below. 1,614 rows cover 74 unique sequences.

### What the fix changed

- **Removed outcome columns from the features.** 23 `label_*` / `endpoint_*` columns were curated annotations of the outcome; `label_kupffer_activation` correlated 0.98 with its own target. Two assay readouts (`cas_avg`, `kd_avg`) were also removed, because they are measured after dosing and don't exist for a new design. 111 constant columns were dropped. 39 features remain.
- **Grouped CV by sequence.** Many rows are replicate measurements of the same oligo, so row-level folds leak. Folds now come from `StratifiedGroupKFold` on the normalized sequence.
- **Recomputed features from sequence.** `tlr7_score`, `ps_x_len`, `lna_x_gc`, `seed_gc` and `cpg_density` had been filled inconsistently across the source datasets, and were zero in about half the rows. They are now derived from the sequence for every training row, the same way the app derives them at scoring time. The shipped models were retrained on the recomputed features, so these numbers differ slightly from the first post-leakage run.
- **Reported with the base rate.** For endpoints with fewer than 10 positive sequences, PR-AUC is noisy. Read the gap to the base rate, not the absolute value.

### The hash embedding, removed

16 of the original 55 features were `embed_stub_hash_rna_00..15`, a sha256 of the sequence. The strategy roadmap described the platform as having three RNA foundation models integrated; what actually reached the models was this hash.

`scripts/ablate_hash_embedding.py` scored every evaluated endpoint with and without the block on the same sequence-grouped folds. Removing it improved **9 of 12** endpoints, mean +0.009 PR-AUC, and off-target gene silencing by **+0.073**. Training now drops the block, leaving 39 features.

The consequence is not cosmetic: off-target gene silencing moved from `weak` to `validated` once the hash features stopped masking it. Hepatocellular injury went from 0.27 to 0.29.

## Baselines

![Model vs baselines with 95% bootstrap CIs](reports/results/baselines.png)

A PR-AUC above the base rate is not evidence that the model learned design chemistry. It could be learning which dataset a row came from, or one crude number like GC content. `scripts/evaluate_baselines.py` scores four baselines against the shipped model on the same sequence-grouped folds, with 95% bootstrap CIs resampled over sequences. The **data source only** baseline sees nothing but which dataset a row came from.

| Endpoint | Base rate | Logistic (39 features) | Gradient boosting | GC only | Length only | Data source only |
|---|---|---|---|---|---|---|
| kupffer cell activation | 0.021 | 0.48 [0.32-0.70] | 0.41 [0.29-0.60] | 0.07 | 0.14 | **0.69 [0.53-0.87]** |
| hepatocellular injury | 0.126 | 0.29 [0.21-0.56] | 0.45 [0.26-0.63] | 0.15 | 0.11 | 0.20 [0.10-0.35] |
| injection site reaction | 0.007 | 0.27 [0.07-0.49] | 0.11 | 0.03 | 0.20 | 0.01 |
| offtarget gene silencing | 0.042 | 0.23 [0.07-0.65] | 0.45 | 0.06 | 0.04 | 0.06 |
| hepatic vacuolation | 0.007 | 0.16 [0.04-0.34] | 0.15 | 0.01 | 0.02 | 0.20 [0.06-0.42] |

All 12 evaluated endpoints are in [`reports/results/baselines.csv`](reports/results/baselines.csv); the within-source comparison is in [`reports/results/within_source.csv`](reports/results/within_source.csv).

### How to read it

- **Why a source indicator predicts anything at all.** 1,572 of the 1,614 rows come from two systematic sources (`oligotox_main`, `papargyri_2020`); the remaining 42 are single rows hand-curated from individual papers and FDA labels. Those 42 exist *because* someone published a toxicity finding about them, so 36 of 42 carry a positive label, against 214 of 1,572 in the systematic sources. "Which source is this row from" is therefore close to the label itself, by construction. That is the selection bias the source-only baseline measures.
- **Kupffer cell activation is a dataset artifact, not a finding.** It has the highest grouped PR-AUC (0.48), but the data-source-only baseline beats it (0.69) and the two CIs overlap. All 34 positive rows (27 sequences) sit in that 42-row curated bucket and none in the two systematic sources, and its base rate there is 0.81. Scored inside that source alone the model gets 0.79, below its own base rate, and GC content alone scores higher (0.87). The model is recognising which dataset a row came from. Do not read this endpoint as a toxicity signal.
- **Hepatocellular injury is the most defensible endpoint.** The 39-feature model (0.27) and gradient boosting (0.37) both clear the data-source baseline (0.20) and the GC-only and length-only baselines, and the signal survives within a single source: inside `papargyri_2020` alone the model gets 0.34 against a 0.25 base rate. It is a modest signal, roughly 2.2× the base rate, on 19 positive sequences.
- **Off-target gene silencing** clears its source-only baseline (0.23 against 0.06) and both simple baselines, on 12 positive sequences. It only emerged once the hash features were removed, and its CI is wide, so treat it as the weakest of the three validated endpoints.
- **Injection site reaction** clears its source-only baseline but not the length-only baseline, so it stays weak and unconfirmed.
- **Endpoints with 4–7 positive sequences** are not reliable estimates, and several are beaten by the length-only baseline. Treat them as placeholders until more labeled data arrives.
- **Gradient boosting matches or beats the shipped logistic model** on several endpoints, but every CI overlaps. With 74 sequences there is no evidence to prefer either, so the shipped model stays the simpler, calibrated one.

## What the app will and won't tell you

The evaluation above is enforced in the app, not just described here. `scripts/derive_evidence_tiers.py` turns `baselines.csv` into `artifacts/config/endpoint_evidence.json`, which the scoring pipeline and the report read at runtime. Each endpoint lands in one tier:

| Tier | Endpoints | How the app presents it |
|---|---|---|
| **Validated** | hepatocellular injury, off-target gene silencing, CNS neurotoxicity | A prediction. Enters the composite score. |
| **Weak signal** | injection site reaction | Model output shown with its CI and the baseline it fails to clear. Not a prediction. |
| **Source-confounded** | kupffer cell activation | Model output shown, labelled as dataset recognition. Not a prediction. |
| **Insufficient data** | the remaining 8 | Model output shown with its positive-sequence count. Not a prediction. |

The practical consequences:

- **The composite score averages one endpoint, not 13,** and the tab says so ("Composite covers 1 of 13 endpoints"). Previously it was a severity-weighted mean over all 13, so eleven unsupported models diluted the one that works. On a test candidate the old composite named *renal tubular injury* — 7 positive sequences, beaten by a length-only baseline — as the dominant endpoint.
- **Unsupported endpoints never return `predicted_positive`.** They return `null` with an `evidence_tier` and a plain-English `evidence_note` saying which baseline beat them.
- **The report states what its own evidence supports.** Every JSON report carries a `conclusions` block naming the evidence standard, what the report supports (relative ranking on the endpoints that met it), and what it does not: an absolute probability of toxicity, a clinical outcome or dose, or a substitute for in vivo toxicology. A report with no validated endpoint supports nothing and says so.
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
| siRNA potency | 2,838 duplexes, 41 genes | **0.57** [0.54-0.59] | **0.46** / **0.49** | 0.30 | -0.11 | validated |
| ASO activity | 3,867 PS-ASOs, 85 genes | **0.26** [0.23-0.29] | none available | 0.14 | -0.15 | weak |

Spearman between measured potency and out-of-fold prediction, using position-specific features (nucleotide identity at positions 1, 10 and the 3' end, seed AU content and duplex-end stability) alongside composition. Adding them lifted every split, held-out studies included. The siRNA number transfers between two independent studies (Ichihara 1 and 2) that share no sequences, in both directions. The gene-mean baseline carries none of it, and 17-29% of the variance is between genes, which is why the grouping matters.

The tab only scores the modality a model was trained for: ask it about a PMO and it says no model is validated for that, rather than returning a number. Scores are raw readouts in the units of the training data, labelled as such, never probabilities.

**Caveat on the siRNA model:** the training data is 2007 chemistry with no modification annotation, so it scores sequence design only and knows nothing about the chemistry of a modern GalNAc-siRNA.

### A duplex parsing bug this turned up

Ichihara's entries are siRNA duplexes written as one string, `sense.antisense`. Read naively they look like 41-nt "sequences" whose length and GC content are computed over both strands and the separator. The normalizer now splits them and uses the 21-nt antisense guide, which is what drives potency.

### Privacy on the public deployment

Runs live in the visitor's browser session and are never written to the server. An earlier version saved every run to `artifacts/runs/` and the Saved runs tab listed every file it found there, which on a shared container means one visitor's sequences are readable by the next. Download a run to keep it.

## Calibration, reported separately from ranking

PR-AUC says whether a model ranks. It says nothing about whether 0.30 means a three-in-ten chance. For a safety tool those come apart, so `scripts/evaluate_calibration.py` reports both on the same out-of-fold predictions.

| Endpoint | Tier | Observed rate | Mean predicted | Brier skill | Calibration MAE |
|---|---|---|---|---|---|
| cns toxicity | validated | 0.412 | 0.404 | **+0.30** | 0.009 |
| offtarget gene silencing | validated | 0.042 | 0.045 | +0.03 | 0.037 |
| hepatocellular injury | validated | 0.126 | 0.111 | **-0.01** | 0.066 |
| injection site reaction | weak | 0.007 | 0.029 | -2.37 | 0.021 |

Brier skill is `1 - brier/brier_base_rate`; at or below zero means the probabilities are no better than quoting the base rate.

**The finding that matters: hepatocellular injury ranks but does not calibrate.** It beats every baseline on PR-AUC and its signal survives a within-source check, yet its Brier skill is slightly negative. Read its score as a ranking, not a percentage. CNS neurotoxicity is the one endpoint whose probabilities are trustworthy as probabilities, which is unsurprising given it has 3,338 training sequences against 74.

Most unsupported endpoints over-predict badly, injection site reaction by about four times its base rate, which is a second reason not to read them as predictions. Every endpoint's `probabilities_calibrated` flag and Brier skill is in `artifacts/config/endpoint_evidence.json`, and the app shows a "reads as" column saying *probability* or *ranking only*.

![Calibration of the validated endpoints](reports/results/calibration.png)

### Do richer representations help? No.

The strategy roadmap lists three RNA foundation models under "what we have". Nothing called them, and what reached the models was the hash stub. With 3,338 CNS sequences and a cross-assay split, the question is finally testable: `scripts/evaluate_sequence_representations.py`.

| Representation | Features | MOE to Hagedorn | Hagedorn to MOE |
|---|---|---|---|
| **Hand features** | 17 | **0.828** [0.808-0.847] | **0.783** [0.761-0.804] |
| 3-mer counts | 64 | 0.791 | 0.772 |
| 3-mer + 4-mer | 320 | 0.792 | 0.767 |
| k-mer + hand | 337 | 0.827 | 0.779 |
| Hash stub (control) | 16 | 0.506 | 0.507 |

Seventeen hand-designed features beat 320 k-mer features in both directions, and adding k-mers to them changes nothing. Sequence content past what the design features already capture does not carry CNS neurotoxicity here. The hash control lands on 0.506, exactly chance, which is both a sanity check on the setup and confirmation that removing it was right.

Pretrained RNA language model embeddings are not evaluated, because no weights are available in this environment; downloading them is a prerequisite, not a result. But the k-mer arms test the thing those embeddings are meant to supply, and it is not what is limiting this endpoint. The roadmap's claim that foundation models are integrated remains untrue, and on this evidence acting on it is not the highest-value next step.

### Hepatocellular injury does not transfer

The one public set pairing gapmer sequences with liver enzymes is the hepatotoxicity screen in [NAR 2022;50(13):7224](https://doi.org/10.1093/nar/gkac562), Supplementary Table S1: 149 LNA gapmers, C57BL/6 mice at 20 mg/kg IV, AST/ALT at 96 h. `scripts/evaluate_hepatotox_external.py` parses it and scores every sequence with the shipped model.

| | Value |
|---|---|
| Sequences | 149 |
| Hepatotoxic (ALT > 100 U/L) | 5 |
| ROC-AUC | **0.40** [0.13-0.70] |
| PR-AUC | 0.035, against a 0.034 base rate |
| Rank of the toxic sequences | 23, 76, 87, 120, 139 of 149 |

**The model ranks the known-hepatotoxic gapmers no better than chance, and three of the five sit in the bottom half.** With 5 positives the interval spans almost everything, so this does not prove the model is wrong. What it does mean is that there is no evidence the endpoint transfers off its own source, and that sits alongside its negative Brier skill.

Read together: hepatocellular injury ranks within `papargyri_2020`, does not produce calibrated probabilities, and shows no transfer to an independent screen of a different chemistry. It keeps its tier, because the tier reflects the internal evaluation and 5 positives should not overturn that either way, but the evidence string the app displays now carries this result.

This is also the clearest statement of what the project needs: not more features, but hepatotoxicity data at the scale the CNS endpoint has.

## The evaluation is the product

The models here are not better than PFRED, ASOptimizer or any proprietary platform, and with 74 unique sequences they will not be. What this repo has instead is a harness that says, per endpoint, whether a result survives scrutiny, and an app that refuses to present anything that does not.

| Check | Script | Catches |
|---|---|---|
| Leakage | `src/oligotox/training/train.py` | Curated annotations and post-dose readouts used as features |
| Grouped splits | same | Replicate rows of one sequence spanning a split |
| Source baseline | `scripts/evaluate_baselines.py` | A model reading which dataset a row came from |
| Ablation | `scripts/ablate_hash_embedding.py` | Feature blocks that carry nothing |
| Representation | `scripts/evaluate_sequence_representations.py` | Whether a richer encoder would help |
| Calibration | `scripts/evaluate_calibration.py` | Scores that rank but are not probabilities |
| External transfer | `scripts/evaluate_cns_external.py`, `evaluate_hepatotox_external.py` | Internal validity that does not generalise |
| Gating | `scripts/derive_evidence_tiers.py` | Any of the above reaching the user as a prediction |

The checks are modality-agnostic, and they live in their own installable package, [`packages/evidence-audit`](packages/evidence-audit), because they need nothing from this project and nothing from biology:

```bash
pip install -e packages/evidence-audit
evidence-audit their_data.csv --label hepatotox --group subject_id --source study
```

It needs a feature matrix, a label column, a column to group on and optionally a source column. Group keys are compared as strings and otherwise left alone; `--nucleotide` additionally folds U onto T so RNA and DNA spellings of one sequence are one group, which is the default on this project's own import path (`python -m oligotox.audit`, kept as a thin wrapper). It exits non-zero when an endpoint fails a check, so it can gate a pipeline. Run blind against this repo's own training data it reproduces, unprompted, the three conclusions that took this project months to reach:

```
               endpoint                verdict  groups  pos_groups  grouped  ungrouped  source_only  brier_skill
  hepatocellular_injury ranks but uncalibrated      74          19    0.274      0.550        0.201        0.002
kupffer_cell_activation      source confounded      74          27    0.480      0.776        0.827        0.343
       thrombocytopenia      insufficient data      74           4    0.100      0.178        0.003        0.031
```

These vary in the third decimal between platforms, which is why the calibration verdict requires a Brier skill of 0.02 rather than merely a positive one. Hepatocellular injury scored -0.009 on one machine and +0.002 on another, and under a sign test that inverted its published verdict; [`NEGATIVE_RESULTS` entry 14](docs/NEGATIVE_RESULTS.md) has the detail.

Each verdict comes with its reasoning, for example *"a model seeing only the source scores 0.83 against 0.48; the model is recognising the dataset"*. Each result lands in `artifacts/config/endpoint_evidence.json`, which the app and the OSWG report read at runtime, so a failed check becomes a label in the interface rather than a footnote.

A checker that fails everything is not a checker. `scripts/demo_audit.py` runs it unchanged over this project's own data and over the external CNS set:

| Dataset | Endpoint | Verdict |
|---|---|---|
| internal, curated | hepatocellular injury | ranks but uncalibrated |
| internal, curated | kupffer cell activation | source confounded |
| internal, curated | thrombocytopenia | insufficient data |
| **external, two published assays** | cns toxicity | **passes** |

Three failures for three different reasons, one pass, no tuning between runs. The verdicts track the data, not the tool's disposition.

`scripts/check_deployment.py` is the companion for the deployed app: it compares the live deployment's branch against this checkout's evidence config and exits non-zero on drift. It exists because a push once left the container serving a stale config, labelling a calibrated endpoint as ranking-only, and nothing caught it until someone opened the app by hand.

**[`docs/NEGATIVE_RESULTS.md`](docs/NEGATIVE_RESULTS.md)** records the ten failure modes this found, with numbers and the script for each. **[`docs/DATA_REQUESTS.md`](docs/DATA_REQUESTS.md)** lists what the platform is actually limited by, which is labelled data rather than method.

## What each failing endpoint would need

"Insufficient data" is true but useless on its own. `scripts/endpoint_requirements.py` works out the binding constraint per endpoint and what would clear it.

The headline is more encouraging than the verdicts suggest. **Eight endpoints are blocked only on count, and 31 additional positive sequences in total would move all eight past the minimum** — 3 to 6 each:

| Endpoint | Positive sequences | Needs |
|---|---|---|
| thrombocytopenia | 4 | 6 more |
| coagulation effect | 5 | 5 more |
| complement activation | 6 | 4 more |
| hepatic vacuolation | 6 | 4 more |
| cytokine release, innate immune activation, renal accumulation, renal tubular injury | 7 each | 3 more each |

Clearing the count is necessary, not sufficient — most also need those positives to come from a *second* source, since all of their current ones sit in a single dataset.

The other two failures are different in kind:

- **Kupffer cell activation** is not short of data: it has 27 positive sequences, all 27 in `curated_other`. It needs positives from a systematic screen rather than hand-curated literature examples, because curated rows exist *because* a finding was published about them.
- **Injection site reaction** has a real estimate that is too noisy. Roughly 109 positive sequences, against 12 now, would be needed for its confidence interval to clear the length-only baseline — an order-of-magnitude estimate that assumes the point estimate holds, which more data may disprove.

Shortfalls against the minimum are exact arithmetic. Sample sizes needed to narrow a confidence interval are flagged `estimate: true` in `reports/results/endpoint_requirements.json`.

### A 190,927-row corpus that made things worse

The literature monitor surfaced [ASO Atlas](https://github.com/barneyhill/aso_atlas) (Apache-2.0, 165,555 unique sequences, 343 genes, with chemistry). Replacing the 3,867-sequence ASO activity model with it should have been the easiest win here. It was not:

| Setup | Spearman |
|---|---|
| McQuisten 2007, 3,867 sequences (shipped) | **0.26** |
| ASO Atlas pooled, 161,832 pairs | 0.18 |
| ASO Atlas transferred to McQuisten | **-0.06** |
| **A431 cell line at a single dose**, 12,435 pairs | **0.29** |

The corpus pools percent inhibition across 81 cell lines and many dosages, extracted from patents spanning 2001-2025, so the readout is not commensurable. One cell line at one dose, on 8% of the data, beats the whole corpus. Adding cell line and dose as features made it *worse*, so the heterogeneity is not an offset a covariate can absorb.

Adopting it needed per-context models rather than a swap, so that is what it got.

### Per-assay potency models

`scripts/train_context_potency.py` trains one model per experimental context and gates each against its own baselines: it must beat a target-gene-mean baseline that ignores the sequence, beat GC alone, and clear a 0.15 floor. Folds are grouped by target gene inside each context.

| Context | Sequences | Genes | Spearman | GC only | Gene mean | |
|---|---|---|---|---|---|---|
| A431, dose 2000, gymnosis | 15,236 | 9 | **0.33** | -0.22 | -0.33 | ships |
| A431, dose 4000, gymnosis | 17,679 | 13 | **0.30** | 0.09 | -0.17 | ships |
| SH-SY5Y, dose 4000, electroporation | 13,115 | 6 | 0.24 | -0.10 | -0.37 | ships |
| HepG2, dose 1000, electroporation | 2,668 | 6 | 0.24 | -0.18 | -0.17 | ships |
| T-24, any dose | 3,664 | 67 | 0.16 | 0.08 | -0.14 | ships |
| HepG2, dose 4500, electroporation | 6,308 | 7 | 0.15 | -0.23 | -0.27 | below floor |
| A549, any dose | 4,131 | 59 | -0.01 | -0.26 | -0.25 | below floor |

**Five of seven ship, and the best reaches 0.33** against 0.26 for the single-study model and 0.18 for the same corpus pooled. Two are rejected, which is the point of having a gate.

Each prediction names its assay: the app shows "A431 knockdown 32%", not a potency score. The caveat on every one says it is not general potency, and a test asserts both.

## Design advice

Scoring a candidate tells you where it stands. The Design advice tab tells you what to change.

Every single-base substitution is applied, re-scored through the same gated models the rest of the app uses, and ranked by how much predicted risk it removes. For a G-rich 20-mer starting at 0.66 predicted CNS risk:

| Edit | Risk reduction | Knockdown change | Trade-off |
|---|---|---|---|
| G17C | 0.172 | -0.1% | none worth speaking of |
| G17A | **0.174** | -3.7% | costs potency |
| G11T | 0.148 | **+0.8%** | improves both |
| G4C | 0.135 | +0.1% | improves both |

Two design decisions behind it.

**Advice comes only from validated endpoints.** Ten of the thirteen toxicity models are not presented as predictions, and an edit recommended to move a number that means nothing is worse than no advice. A test asserts the advised set equals the validated set.

**Nothing is read off a coefficient.** The models are linear, so coefficients exist, but with correlated features they describe the fit rather than the effect of an edit. Every number here is the difference between two real predictions, and a test re-scores the edited sequence to confirm the claimed change is the one that actually occurs.

Trade-offs are shown rather than hidden: an edit that lowers risk at the cost of knockdown says so, and the few that improve both are called out. As a sanity check the advice targets G residues, which is what the published rule the CNS endpoint reproduces would predict.

## Structure, accessibility and off-target burden

Three mechanistic features, two of which were placeholders.

**Real folding energies.** `struct_mfe_proxy` was computed from GC content and `struct_viennarna_used` was 0.0 on every row, so the shipped structure block contained no thermodynamics at all. ViennaRNA is now a pinned runtime requirement, and the code raises rather than falling back silently, because a placeholder where real values are expected is drift that nothing would catch.

**Off-target complementarity.** The published mechanism for gapmer hepatotoxicity is RNase H1 cleaving unintended transcripts, so complementarity burden should predict it. `scripts/offtarget_transcriptome.py` matches three reverse-complement probes per oligo in one Aho-Corasick pass over the whole Ensembl cDNA set, 465,769 transcripts and 1,276 Mb. The previous implementation returned `0.01 + gc * 0.02` — GC content rescaled and labelled a hit count.

It did not work: off-target burden alone scores **0.125 against a 0.126 base rate**. A dimensionality control makes the result honest — six random columns cost 0.030 and the off-target block costs 0.052, so most of the drop is the penalty of adding any six features to 19 positive sequences.

**Target-site accessibility.** P(binding site unpaired) from ViennaRNA's local partition function. Measured against cDNA it looked worthwhile, 0.100 to 0.157. Measured against genomic spans, where 93% of oligos locate instead of 34%, the gain falls to **+0.014** and accessibility alone scores −0.023. The first number was inflated by a coverage bias: the oligos that happen to locate in a cDNA are exonic, which is a selected subset. See `NEGATIVE_RESULTS.md` entry 12.

### The finding that explains both

Only 34% of ASO Atlas oligos locate in a canonical cDNA. Orientation is not the problem: the reverse complement matches and the forward sequence does not. Fetching the genomic span for DGAT2 and relocating the same 400 ASOs:

| Target sequence | ASOs located |
|---|---|
| canonical cDNA, 2,407 nt | 96 of 400 (24%) |
| **genomic span, 54,354 nt** | **400 of 400 (100%)** |

Gapmers target pre-mRNA, and a cDNA has no introns. So the off-target burden above was measured against mature mRNA while the published mechanism implicates long pre-mRNAs — the wrong molecule, which is a different statement from the mechanism being wrong. `fetch_transcripts.py --genomic` fetches the right sequence and accessibility folds a bounded window around the site, so a megabase pre-mRNA stays tractable.

Fixing it raised accessibility coverage from 34% to 93% across 59 genes, and in doing so **cut the measured benefit of accessibility from +0.057 to +0.014**. A missing-data problem that correlates with the feature being tested inflates that feature's apparent value, so low coverage was a bias problem rather than only a power problem.

## Continuous evaluation

The platform is meant to absorb new data without anyone watching, which is dangerous in the obvious way: a loop that retrains and promotes on its own is how the original PR-AUC 0.93 happened. So the loop is built as a ratchet. Evidence can promote an endpoint and evidence can demote it, but the gates never loosen and nothing can pin a tier.

```bash
python scripts/fetch_datasets.py          # what is declared and usable
python scripts/run_pipeline.py            # normalize, train, evaluate, re-tier
python scripts/run_pipeline.py --check    # non-zero if any tier moved
```

One run does: normalize declared sources, retrain the internal endpoints, recompute baselines and calibration, re-derive every tier from that evidence, compare against the previous run, and append what changed to `reports/results/endpoint_history.jsonl`.

**A source is only ingested if `data/registry/sources.json` declares its readout direction and label threshold.** A downloader must not decide what counts as toxic, so a dataset without a declared direction is reported as blocked rather than guessed at. Four sources are usable today, 9,637 sequences; four are blocked, each with the reason.

**Nothing in the loop can promote an endpoint.** Tiers are recomputed by `derive_evidence_tiers.py` from the measurements, applying the same gates every time. Demotion is automatic and immediate: an endpoint that stops clearing its baselines loses its tier on the next run, with no grace period and no override. There is deliberately no `force_tier` and a test asserts nobody adds one.

### The inbound stream

A loop with no new data recomputes the same answer forever. `scripts/monitor_literature.py` queries Europe PMC weekly for papers that plausibly carry an ingestible dataset, and scores them on what makes a dataset usable rather than interesting: supplementary data, open access, panel or screen language rather than a single compound, and a named toxicity readout. Anything resembling a source already ingested is scored down, and every candidate is deduplicated against previous runs, so a quiet week reports nothing rather than the same papers again.

It reports leads, never sources. A paper becomes ingestible only once someone reads it, establishes the readout direction and threshold, and declares it in `data/registry/sources.json`. That judgement is exactly what the registry exists to keep out of automation, and getting it wrong is how the shipped combined file ended up calling 83% of Hagedorn's compounds neurotoxic.

### Running while you are offline

`.github/workflows/continuous-evaluation.yml` runs the chain weekly in GitHub's cloud, so it keeps going with no machine involved. The split matters:

| Outcome | What the workflow does |
|---|---|
| No tier changed | Commits the refreshed evidence directly |
| Any tier moved | Opens a **pull request** instead |

A tier change alters which endpoints read as predictions, which is a claim about what the platform can predict, so it gets a human. The PR body says what to check and names the distinction worth catching: a promotion driven by more data is the point, a promotion driven by a new feature that happens to fit 74 sequences is the failure this project exists to catch.

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

The data every published number is computed from is in the repo, so the evaluation is checkable rather than claimed. Four files, 1.2 MB total: `feature_matrix.csv` (1,614 x 39), `label_matrix.csv`, `training_dataset.csv` (sequences, source and curated labels) and `dataset_metadata.json`. Every source in them is public literature, a PMID/PMC paper or an FDA label; the other contents of `data/processed/` are intermediates and stay out, as does `data/raw/`, which holds source papers and third-party datasets this repository has no right to redistribute.

From a fresh clone, training takes about 10 seconds and the baselines about 95 seconds. The scripts need a few libraries the app does not, kept in a separate extra so the deployed app stays small:

```bash
pip install -e ".[scripts]"             # matplotlib, pypdf, xlrd, openpyxl
python -m oligotox.training.train       # retrains all 13 endpoint models
python scripts/evaluate_baselines.py    # recomputes baselines and bootstrap CIs
```

Verified in a clean clone of this commit: all 13 per-endpoint PR-AUCs reproduce to the published figure exactly, `baselines.csv` and `within_source.csv` reproduce with zero difference across all 60 and 36 rows, `endpoint_evidence.json` comes out byte-identical, and so do all 13 shipped `.joblib` models when retrained on the same platform. Across platforms the models differ in bytes but not in behaviour: a scheduled Linux run reproduced macOS-built models to within 3e-9 in predicted probability, which is BLAS float noise rather than a different fit. The models in `artifacts/` are the ones this data and this script produce.

### Retrain and regenerate the figures

```bash
pip install -e ".[scripts]"              # figures and parsers need this extra
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
| Efficacy | Potency for the candidate's modality, with units, evidence and caveat, plus per-assay knockdown models | siRNA validated, ASO per-assay |
| Design advice | Single-base edits ranked by predicted risk reduction, with their potency cost | Validated endpoints only |
| Saved Runs | Runs kept during the session | — |
| Build Info | Model versions, feature count and config provenance | — |

The Efficacy tab previously held untrained model slots and was removed rather than shipped behind a disclaimer. It returned once the models were trained and validated against baselines.

## Limitations

- 74 unique sequences is small. Nine of the 13 endpoints have 7 or fewer positive sequences.
- **Only 3 of the 13 endpoints are presented as predictions.** Hepatocellular injury survives both the data-source baseline and a within-source check; Kupffer cell activation, the highest raw score, does not. CNS neurotoxicity is validated on external data but reproduces a published rule rather than improving on it. The other 10 models ship and run, but the app labels them unsupported.
- Hepatocellular injury rests on 19 positive sequences with within-source confirmation from a single dataset. A second independent source is what would settle it, and no external hepatotoxicity data is on hand.
- Several endpoints draw all their positives from one source dataset, so source and biology cannot be separated for them with this data.
- The structure features use a placeholder MFE. The hash 'sequence embedding' was removed after an ablation showed it helped nothing; no learned sequence model has replaced it.
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

## Citing this work

`CITATION.cff` carries the metadata; GitHub's "Cite this repository" button renders it. If you use the evaluation findings rather than the code, the nine failure modes in [`docs/NEGATIVE_RESULTS.md`](docs/NEGATIVE_RESULTS.md) are the citable part, and each names the script that reproduces it.

Cite the underlying studies and OligoGym for the data, not this repository; see `data/oligogym/ATTRIBUTION.md`.

## License

MIT; see `LICENSE`.

The curated datasets under `data/` are derived from published literature, PMC articles and FDA labels, each credited in `docs/provenance/DATA_AND_MODEL_PROVENANCE.md`. The MIT grant covers this repository's own code and curation; the underlying publications keep their own terms, so cite the original sources rather than this repo when using the data.
