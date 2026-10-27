# OligoTox: Data & Model Provenance

**Revised**: 2026-10-27 (supersedes the June 26, 2026 version in full)
**Maintainer**: Beshoy Armanios, PhD (Pharmacology/Toxicology)

The June 2026 version of this file was written before label leakage was found
in the training data. It advertised 99.9% accuracy, RandomForest models, 1,576
sequences and 49 features. All of that was wrong or has since changed, and the
numbers in it came from leaked outcome columns. It has been rewritten against
the code and artifacts as they stand. The root `README.md` remains
authoritative for evaluation results; this file covers where the data and
models came from.

---

## 📖 PRIMARY DATA SOURCES & CITATIONS

### Published Datasets

#### 1. **scASO (Mortberg et al. 2023)**
**Full Citation**: Mortberg MA & Gentile JE (equal contribution), Nadaf NM, Vanderburg C, Simmons S, Dubinsky D, Slamin A, Maldonado S, Petersen CL, Jones N, Kordasiewicz HB, Zhao HT, Vallabh SM, Minikel EV. **A single-cell map of antisense oligonucleotide activity in the brain.** _Nucleic Acids Res._ 2023 May 16:gkad371.

- **DOI**: [10.1093/nar/gkad371](https://doi.org/10.1093/nar/gkad371)
- **PMID**: [37188501](https://pubmed.ncbi.nlm.nih.gov/37188501/)
- **Repository**: `data/raw/github_repos/scaso/`
- **Data Type**: Single-cell RNA-seq, ASO knockdown efficiency
- **Species**: Mouse (brain tissue - multiple cell types)
- **Key Contribution**: Cell-type-specific ASO activity maps, mechanistic insights
- **Integration Point**: `src/oligotox/omics/` (single-cell DEG analysis)

#### 2. **openASO (McQuisten & Peek 2007)**
**Full Citation**: McQuisten KA, Peek AS. "Identification of sequence motifs significantly associated with antisense activity." _BMC Bioinformatics._ 2007;8:184.

- **DOI**: 10.1186/1471-2105-8-184  
- **PMID**: [17555590](https://www.ncbi.nlm.nih.gov/pubmed/17555590)
- **Repository**: `data/raw/github_repos/openASO/`
- **Data Type**: Sequence-activity relationships, motif analysis
- **Key Contribution**: Predictive motifs for ASO potency
- **Integration Point**: `src/oligotox/features/aso_rich_features.py`

#### 3. **OligoGym Benchmark (Martinelli et al. 2023)**
**Description**: "A curated siRNA potency dataset that contain only sugar and base information."

- **Conference**: NeurIPS 2025 Track (Roche collaboration)
- **Repository**: `data/raw/oligogym/oligogym_bulk/Martinelli_2023_1_processed.json`
- **Data Type**: siRNA potency with chemical modification annotations
- **Key Contribution**: Gold-standard benchmark for chemistry-activity models
- **Integration Point**: `src/oligotox/benchmarks/oligogym.py`
- **Cache Directory**: `./data/benchmarks/`

#### 4. **ChEMBL/SIDER Approved Oligonucleotide Drugs**
**Source**: FDA-approved therapeutics database

- **Module**: `src/oligotox/integrations/chembl_benchmark.py`
- **Drug Count**: 13 approved oligonucleotide drugs
- **Data Type**: Real-world toxicity profiles, BLACK BOX warnings, REMS programs
- **Key Examples**:
  - **Mipomersen** (Kynamro): PS-ASO, BLACK BOX for thrombocytopenia + glomerulonephritis
  - **Eteplirsen** (Exondys 51): PMO backbone, lower toxicity profile benchmark
  - **Fomivirsen** (Vitravene): First approved ASO (discontinued)
- **Use Case**: External validation against clinical safety data

---

## TRAINING DATA

**File**: `data/processed/training_dataset.csv` (committed, 335 KB)

- **1,614 rows covering 74 unique sequences.** Many rows are replicate
  measurements of the same oligo. 74, not 1,614, is the effective sample size,
  and all cross-validation is grouped on the normalized sequence (RNA and DNA
  spellings of one sequence count as the same group).
- **13 binary toxicity endpoints** in `label_matrix.csv`. Nine have 7 or fewer
  positive sequences; CNS toxicity has 2 and is not evaluated.
- **`source`** names the originating dataset per row.

### Source composition and its selection bias

| Bucket | Rows | Rows positive on >=1 endpoint |
|---|---|---|
| `oligotox_main` | 800 | 214 of 1,572 across both systematic sources |
| `papargyri_2020` | 772 | (same 214) |
| individually curated (35 distinct sources) | 42 | 36 |

The 42 individually curated rows come from single papers, PMC articles and FDA
labels, and they exist *because* a toxicity finding was published about them.
So 86% of them carry a positive label, against 14% of the systematic rows. The
consequence is that "which source is this row from" is close to the label
itself, by construction. This is measured, not assumed: the root README's
Baselines section scores a data-source-only baseline against the real model,
and for Kupffer cell activation the source baseline wins (0.69 vs 0.48). All
34 Kupffer positive rows sit in the 42-row curated bucket and none in the
systematic sources.

---

## FEATURES

`build_training_dataset` assembles 191 candidate feature columns; training
filters them to **55**. The filter is in `artifacts/metrics/dropped_features.json`:

| Dropped | Count | Why |
|---|---|---|
| Outcome annotations (`label_*`, `endpoint_*`) | 23 | Curated annotations of the outcome. `label_kupffer_activation` correlated 0.98 with its own target. This was the leakage. |
| Assay readouts (`cas_avg`, `kd_avg`) | 2 | Measured after dosing; they do not exist for a new design. |
| Constant columns | 111 | No variance across rows, so no information. |

### The 55 surviving features

| Block | Count | Contents |
|---|---|---|
| Design and chemistry | 9 | `ps_fraction`, `lna_fraction`, `seq_length`, `gc_content`, `seed_gc`, `cpg_density`, `tlr7_score`, `ps_x_len`, `lna_x_gc` |
| `aso_*` | 22 | Chemistry and modality one-hots (LNA, MoE, 2'-OMe, PS, GalNAc; gapmer, siRNA, PMO, SSO, mRNA) and derived flags |
| `struct_*` | 8 | Length, GC, uracil fraction, longest homopolymer, palindrome fraction, dinucleotide repeats, hairpin score, `struct_mfe_proxy` |
| `embed_stub_hash_rna_00..15` | 16 | A hash of the sequence, not a learned embedding |

Two caveats a reader should not have to discover from the code:

- **`struct_mfe_proxy` is a placeholder, not a folding energy.** ViennaRNA is
  requested at build time (`viennarna_requested: true`) but never used
  (`struct_viennarna_used` is 0.0 for every row), so the structure block
  carries no thermodynamic information.
- **16 of the 55 features (29%) are a hash stub.** The column names say so.
  There is no trained sequence model behind them.

### Omics features add nothing to the shipped models

`build_training_dataset` loads 4 omics modalities and joins 66 omics feature
columns (`omics_features_added: 66`, `candidate_level_omics_joined: true`).
Because the join is at candidate level, every one of those columns is constant
across rows, and all 66 are dropped by the constant filter. **No omics, GSEA,
pathway, sentinel-marker or Geneformer feature reaches any shipped model.**
The modules under `src/oligotox/omics/` and `src/oligotox/integrations/` run
and are reachable from the API, but they do not contribute to the 55 features
the models are trained on.

---

## MODELS

**Location**: `artifacts/models/<endpoint>.joblib` (13 files, ~6.8 KB each)

- **Algorithm**: `CalibratedClassifierCV(LogisticRegression(max_iter=2000,
  class_weight="balanced", solver="liblinear"), method="sigmoid", cv=2)`
- **Defined in**: `_make_model()` in `src/oligotox/training/train.py`
- **Not RandomForest.** The earlier RandomForest models, their
  `endpoint_*.joblib` filenames and their reported accuracies are gone.

### Why accuracy is not reported

Base rates run from 0.126 down to 0.0012. A model predicting "no toxicity"
for every oligo scores above 99% accuracy on most of these endpoints, which is
how the old "99.9% average accuracy" claim arose. Evaluation uses PR-AUC
against the base rate, with 95% bootstrap CIs over sequences, plus baselines.
Current results, including the finding that only 1 of 13 endpoints is
defensible, are in the root `README.md`.

### Evidence gating

`artifacts/config/endpoint_evidence.json`, generated by
`scripts/derive_evidence_tiers.py` from the baseline results, assigns each
endpoint a tier. The app and the OSWG report read it at runtime: only
`validated` endpoints are presented as predictions, enter the composite score,
or can raise a regulatory recommendation.

---

## PIPELINE

**Training entry point**: `python -m oligotox.training.train`

```
1. build_training_dataset.run()        -> training_dataset.csv, feature_matrix.csv,
                                          label_matrix.csv, dataset_metadata.json
     features/aso_rich_features.py        chemistry and modality
     features/structure_features.py       sequence structure (placeholder MFE)
     features/embedding_features.py       hash stub
     omics/ + integrations/               joined, then dropped as constant
2. Drop outcome annotations, assay readouts and constant columns (55 remain)
3. Recompute tlr7_score, ps_x_len, lna_x_gc, seed_gc, cpg_density from the
   sequence, so training features match what the app derives at scoring time
4. Per endpoint: StratifiedGroupKFold on normalized sequence, out-of-fold
   predict_proba, PR-AUC and Brier
5. pipelines/thresholds.batch_select_thresholds() -> endpoint_thresholds.json
6. Save models, composite_config.json, training_summary.json, build_info.json,
   dropped_features.json
7. DataSAIL split kept as a diagnostic only (datasail_split_diagnostic.json)
```

**Cross-validation**: `StratifiedGroupKFold` grouped on the normalized
sequence. DataSAIL (`integrations_v03/datasail_splits.py`) is run for
comparison and is not the basis of any published number.

---

## REPRODUCIBILITY

Four files are committed so the published numbers can be checked:
`feature_matrix.csv`, `label_matrix.csv`, `training_dataset.csv` and
`dataset_metadata.json` (1.2 MB total). The rest of `data/processed/` (445 MB)
and `data/raw/` (160 MB) are intermediates and stay out of the repo.

```bash
python -m oligotox.training.train       # ~10 s
python scripts/evaluate_baselines.py    # ~95 s
python scripts/derive_evidence_tiers.py
```

Verified in a clean clone at commit `8d31e2c`: all 13 per-endpoint PR-AUCs
reproduce exactly, `baselines.csv` and `within_source.csv` reproduce with zero
difference over all 60 and 36 rows, `endpoint_evidence.json` is byte-identical,
and all 13 `.joblib` models are byte-identical. The shipped models are the
ones this data and this script produce.

### What is not versioned

`artifacts/models/` holds only the 13 current `.joblib` files. The timestamped
run directories, `model_manifest.json`, `models_index.json` and
`artifacts/models_v3/` described in the June version do not exist. Model
provenance is the git history plus `artifacts/metrics/build_info.json`, which
records the training CSV path, row count and the feature and label columns.

---

## MODULES THAT EXIST BUT DO NOT AFFECT PREDICTIONS

Present, importable, and outside the shipped prediction path. Listed so a
reader does not mistake the repository's surface area for validated
capability:

| Module | Status |
|---|---|
| `integrations/geneformer_context.py` | Curated proxy, no Geneformer model |
| `integrations/seed_match.py` | `score_seed_matches()` takes an optional transcriptome sample; there is no bundled 3'UTR database |
| `integrations/gsea_pathway.py` | Runs on a DEG table; its outputs are not training features |
| `integrations/candidate_mechanism.py` | Assembles the above for display |
| `integrations/chembl_benchmark.py` | 13 approved-drug comparison entries |
| `benchmarks/oligogym.py` | OligoGym loader; see `benchmarks/oligogym/README.md` |
| `models/comparator.py` | Ranks candidates by composite score |
| `inference/api.py` | `run_oligotox()`; the Streamlit app calls `score_candidate()` directly and does not use it |
| `efficacy/` | Model interfaces only. The Efficacy tab was removed rather than shipped untrained |

---

## CONTACT

- **Maintainer**: Beshoy Armanios, PhD
- **Repository**: https://github.com/barmanios123/OligoTox
- **License**: MIT; see `LICENSE`. The datasets above are derived from the
  published sources cited in this document and keep their own terms, so cite
  the original publications rather than this repository.
