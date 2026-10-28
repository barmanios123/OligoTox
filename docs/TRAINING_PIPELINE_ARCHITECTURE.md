# OligoTox Training Pipeline

**Revised 2026-10-27.** This replaces a version that described a RandomForest
pipeline evaluated with a random 80/20 split. That split was the leakage this
project later found and fixed, so the old document described the bug as the
design. Nothing in it should be followed.

The root `README.md` is authoritative for results.
`docs/provenance/DATA_AND_MODEL_PROVENANCE.md` covers where the data came
from. This file covers how training runs.

---

## What is being trained

13 binary classifiers, one per toxicity endpoint, from design-time features
only: sequence, chemistry and modality. Nothing measured after dosing.

12 of them train on the internal curated dataset. **CNS neurotoxicity does
not**: it is trained separately on external data with its own feature list,
and the trainer skips it by name. See `scripts/train_cns_external.py`.

Entry point:

```bash
python -m oligotox.training.train
```

About 10 seconds on the committed data.

---

## The data

**`data/processed/training_dataset.csv`** — 1,614 rows covering **74 unique
sequences**. The row count is not the sample size. Most rows are replicate
measurements of the same oligo, which is why every fold is grouped by
sequence.

**`data/processed/label_matrix.csv`** — 13 binary endpoints. Nine have 7 or
fewer positive sequences. CNS toxicity has 2 internally, which is why it is
trained externally instead.

`source` names the originating dataset per row. It is used for the
source-confounding baseline, never as a feature.

---

## Features: 191 candidates, 39 survive

`build_training_dataset` assembles 191 candidate columns. Training filters
them in `_drop_leaky_features`, and the result is written to
`artifacts/metrics/dropped_features.json`:

| Dropped | Count | Why |
|---|---|---|
| Outcome annotations (`label_*`, `endpoint_*`) | 23 | Curated annotations of the outcome. `label_kupffer_activation` correlated 0.98 with its own target. This was the leakage behind the original 0.93. |
| Assay readouts (`cas_avg`, `kd_avg`) | 2 | Measured after dosing; they do not exist for a design being screened. |
| Hash embedding (`embed_stub_hash_rna_*`) | 16 | A sha256 of the sequence. Removing it improved 9 of 12 endpoints. |
| Constant columns | 111 | No variance, so no information. Includes all 66 joined omics columns, constant because the join is at candidate level. |

### The 39 that remain

| Block | Count | Contents |
|---|---|---|
| Design and chemistry | 9 | `ps_fraction`, `lna_fraction`, `seq_length`, `gc_content`, `seed_gc`, `cpg_density`, `tlr7_score`, `ps_x_len`, `lna_x_gc` |
| `aso_*` | 22 | Chemistry and modality one-hots (LNA, MoE, 2'-OMe, PS, GalNAc; gapmer, siRNA, PMO, SSO, mRNA) and derived flags |
| `struct_*` | 8 | Length, GC, uracil fraction, longest homopolymer, palindrome fraction, dinucleotide repeats, hairpin score, `struct_mfe_proxy` |

Two caveats a reader should not have to find in the code:
`struct_mfe_proxy` is a placeholder rather than a folding energy, since
ViennaRNA is requested at build time but never used; and no omics, GSEA,
pathway or foundation-model feature reaches any shipped model.

### Features are recomputed, not trusted

`tlr7_score`, `ps_x_len`, `lna_x_gc`, `seed_gc` and `cpg_density` were filled
inconsistently by the upstream sources and were zero in about half the rows.
Training recomputes all five from the sequence using the same function the
app calls at scoring time, so training and scoring features cannot drift.

---

## Model

```python
from sklearn.calibration import CalibratedClassifierCV
from sklearn.linear_model import LogisticRegression

base = LogisticRegression(
    max_iter=2000,
    class_weight="balanced",   # upweight the rare positive class
    solver="liblinear",
)
model = CalibratedClassifierCV(base, method="sigmoid", cv=2)
```

`_make_model()` in `src/oligotox/training/train.py`. Each fitted model is
about 6.8 KB.

**Why linear.** With 74 unique sequences, a calibrated linear model is what
the data supports. Gradient boosting was evaluated as a baseline and matches
or beats it on several endpoints, but every confidence interval overlaps, so
there is no evidence to prefer either.

---

## Cross-validation

```python
from sklearn.model_selection import StratifiedGroupKFold, cross_val_predict

# RNA and DNA spellings of one sequence are the same group.
groups = raw["sequence"].str.upper().str.replace("U", "T", regex=False)

cv = StratifiedGroupKFold(n_splits=n_splits)   # 4 or 5, by positive count
prob = cross_val_predict(_make_model(), X, y, cv=cv, groups=groups,
                         method="predict_proba")[:, 1]
```

**There is no held-out test set.** With 4 to 27 positive sequences per
endpoint, a single split would be noise. Every published number is
out-of-fold under sequence-grouped CV, reported with its base rate and a 95%
bootstrap confidence interval resampled over sequences.

A DataSAIL split is also computed, to
`artifacts/metrics/datasail_split_diagnostic.json`. It is a diagnostic; no
published number comes from it.

---

## What a run produces

```
artifacts/models/<endpoint>.joblib        12 models; cns_toxicity is skipped
artifacts/config/endpoint_thresholds.json per-endpoint F1-selected threshold
artifacts/config/composite_config.json    severity weights, risk bands,
                                          trained feature columns
artifacts/metrics/training_summary.json   PR-AUC, Brier, positives, n_splits
artifacts/metrics/build_info.json         training CSV, row count, columns
artifacts/metrics/dropped_features.json   what was filtered and why
```

Training does **not** decide what the app presents. That comes from
`artifacts/config/endpoint_evidence.json`, written by
`scripts/derive_evidence_tiers.py` from the baseline, calibration and
external-transfer results. Rerun it after any change to the baselines, or the
tier test fails.

---

## Inference

```
oligotox_app.py
  └─ features.candidate.candidate_features(sequence, modality, chemistry, ...)
  └─ pipelines.score_candidate.score_candidate()
       ├─ loads endpoint_evidence.json, endpoint_features.json, thresholds
       ├─ per endpoint: builds a row from that endpoint's column list
       │    (cns_toxicity has its own sequence-only list; the rest share 39)
       ├─ attaches tier, evidence note and calibration status
       └─ composite over validated endpoints only
```

An endpoint that is not `validated` returns `predicted_positive: None` rather
than a call, and an endpoint whose Brier skill is at or below zero is
labelled "ranking only" in the interface.

---

## Reproducing

```bash
pip install -e ".[scripts]"
python -m oligotox.training.train
python scripts/evaluate_baselines.py        # ~95 s
python scripts/derive_evidence_tiers.py
```

Verified in a clean clone: all 13 per-endpoint PR-AUCs reproduce exactly,
`baselines.csv` and `within_source.csv` with zero difference, and all 13
`.joblib` models byte-identical.

---

## What this pipeline does not do

- It does not produce a usable model for 10 of the 13 endpoints. The
  evaluation says so and the app refuses to present them.
- It does not fit the PK/PD simulator, which has never been validated
  against measured exposure data.
- It does not use any learned sequence representation. Hand-designed
  features beat k-mer representations on an external split; see
  `scripts/evaluate_sequence_representations.py`.

`docs/NEGATIVE_RESULTS.md` records what went wrong and how each was caught.
