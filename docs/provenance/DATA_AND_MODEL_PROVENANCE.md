# 📊 OligoTox: Data & Model Provenance Documentation

**Generated**: June 26, 2026
**Version**: v1.0 Production Release
**Maintainer**: Beshoy Armanios, PhD (Pharmacology/Toxicology)

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

## 🧬 TRAINED MODEL REGISTRY

### Production Models (v4 - Current)
**Storage Location**: `artifacts/models/`

**Training Date**: June 26, 2026
**Training Script**: `artifacts/models/latest_train_baseline/` pretrained bundle (root) + `src/oligotox/training/train.py`
**Algorithm**: RandomForestClassifier (scikit-learn)
**Hyperparameters**:
- `n_estimators=100`
- `max_depth=10`
- `random_state=42`
- `n_jobs=-1`

#### Trained Endpoint Models:
```
artifacts/models/
├── endpoint_cns_neurotoxicity.joblib        117 KB  (acc=1.000)
├── endpoint_coagulation_effect.joblib        148 KB  (acc=0.999)
├── endpoint_complement_activation.joblib     165 KB  (acc=0.999)
├── endpoint_cytokine_release.joblib          159 KB  (acc=0.999)
├── endpoint_hepatocellular_injury.joblib     336 KB  (acc=0.999)
├── endpoint_injection_site_reaction.joblib   191 KB  (acc=0.998)
├── endpoint_innate_immune_activation.joblib  159 KB  (acc=0.999)
├── endpoint_renal_accumulation.joblib        171 KB  (acc=0.999)
├── endpoint_renal_tubular_injury.joblib      181 KB  (acc=0.999)
└── endpoint_thrombocytopenia.joblib          132 KB  (acc=1.000)

TOTAL: 10 models, ~1.7 MB, 99.9% average accuracy
```

#### Additional Variant Models:
- `hepatic_vacuolation.joblib`
- `kupffer_cell_activation.joblib`
- `offtarget_gene_silencing.joblib`

**Total Artifacts Registered**: 13 models (detected by UI Model Registry)

---

## 🗃️ TRAINING DATASET CONSTRUCTION

### Dataset Metadata
- **File**: `data/processed/training_dataset.csv`
- **Sequences**: 1,576 unique oligonucleotides
- **Features**: 49 engineered features
- **Endpoints**: 13 binary toxicity endpoints
- **Metadata**: `data/processed/dataset_metadata.json`

### Feature Categories (49 total)

#### 1. Structural Features (15 features)
**Module**: `src/oligotox/features/structure_features.py`
- GC content, melting temperature (Tm)
- Secondary structure propensity
- Hairpin stability, self-dimer formation
- Length, molecular weight

#### 2. Chemical Modification Features (8 features)
**Module**: `src/oligotox/features/aso_rich_features.py`
- **PS fraction**: Phosphorothioate backbone percentage
- **LNA fraction**: Locked nucleic acid content
- **2'-OMe fraction**: 2'-O-methyl sugar modifications
- **CAS average**: Chemical Analog Score (similarity to known toxic structures)
- **Kd average**: Average binding affinity (dissociation constant)

#### 3. Mechanistic/Pathway Features (18 features)
**Module**: `src/oligotox/integrations/`
- **TLR activation scores**: TLR7/8/9 immune stimulation (`tlr_immune.py`)
- **Seed match scores**: miRNA-like off-target potential (`seed_match.py`)
- **GSEA pathway enrichment**: Enriched biological pathways (`gsea_pathway.py`)
- **Sentinel markers**: Tissue-specific toxicity biomarkers (`omics/feature_dictionary.py`)

#### 4. Embedding Features (8 features)
**Module**: `src/oligotox/features/embedding_features.py`
- k-mer frequency vectors
- Chemistry-aware sequence encodings

### Endpoint Label Curation
**Mapping**: `CURATED_ENDPOINT_MAP` (line 297 in `build_training_dataset.py`)

```python
CURATED_ENDPOINT_MAP = {
    'endpoint_hepatocellular_injury': 'curated_hepatocellular_injury',
    'endpoint_renal_tubular_injury': 'curated_renal_tubular_injury',
    'endpoint_thrombocytopenia': 'curated_thrombocytopenia',
    'endpoint_coagulation_effect': 'curated_coagulation',
    'endpoint_complement_activation': 'curated_complement',
    'endpoint_innate_immune_activation': 'curated_innate_immune',
    'endpoint_cytokine_release': 'curated_cytokine',
    'endpoint_injection_site_reaction': 'curated_injection_site',
    'endpoint_cns_neurotoxicity': 'curated_cns_tox',
    'endpoint_renal_accumulation': 'curated_renal_accum',
}
```

**Normalization**: `_normalize_label_series()` converts categorical/numeric to binary (0/1)

---

## 🔗 DATA PIPELINE WORKFLOW

### Stage 1: Raw Data Ingestion
```
data/raw/
├── github_repos/
│   ├── scaso/          # Mortberg 2023 - scRNA-seq + ASO
│   └── openASO/        # McQuisten 2007 - Motifs
├── oligogym/
│   └── oligogym_bulk/  # Martinelli 2023 - siRNA benchmark
└── external_repos/     # Additional curated sources
```

### Stage 2: Omics Integration
**Module**: `src/oligotox/omics/`

**Registry**: `OMICS_REGISTRY` (defined in `__init__.py`)
- Centralized catalog of all omics data sources
- Metadata: modality, tissue, species, publication

**Loaders**: `loaders.py`
- `load_omics_table(name)` - Load specific omics dataset
- `join_omics_tables(tables)` - Merge multiple modalities
- Reducers: `passthrough`, `variance_topk`

**Quality Control**: `deg_qc.py`
- Expression threshold validation
- Replicate consistency checks
- P-value/FDR filtering

**Gene Symbol Normalization**: `symbol_cleaner.py`
- Uses `GENE_ALIAS_MAP` from `reference/gene_aliases.py`
- Resolves ambiguous gene names

**Mechanistic Readout Conversion**: `deg_to_readout.py`
- Converts DEG tables → `OmicsReadout` objects
- Maps to pathway features (`PATHWAY_FEATURES`, `SENTINEL_MARKERS`)

**Supported Omics Modalities** (`omics/schemas.py`):
1. RNA-seq (bulk, single-cell)
2. Proteomics (mass spec, antibody arrays)
3. Metabolomics (LC-MS, NMR)
4. miRNA profiling

### Stage 3: Feature Engineering
**Main Builder**: `src/oligotox/pipelines/build_training_dataset.py`

**Workflow** (run function, line 159):
1. **Load raw CSV sources** from `data/raw/`
2. **Join omics tables**: `join_omics_tables(OMICS_REGISTRY)`
3. **Build rich features**: `build_aso_rich_features()`
   - Structural: `build_structure_features()`
   - Embeddings: `build_embedding_features()`
   - Chemical mods from input columns
4. **Safe left join on candidate**: `safe_left_join_on_candidate()`
5. **Normalize endpoint labels**: Apply `CURATED_ENDPOINT_MAP`
6. **Export**:
   - CSV: `data/processed/training_dataset.csv`
   - Metadata: `data/processed/dataset_metadata.json`

### Stage 4: Model Training
**Script**: `src/oligotox/training/train.py`

**Full Training Pipeline** (train_all function, line 185):
```python
1. build_info = build_training_dataset(
       "data/processed/training_dataset.csv",
       ...,
       "data/processed/dataset_metadata.json"
   )
2. Load metadata: dataset_metadata.json
3. Build v03 feature DataFrame: build_v03_feature_df()
4. DataSAIL cross-validation splits: datasail_splits.py
5. Train RandomForest per endpoint:
   - Fit on training fold
   - Evaluate on validation fold
   - Threshold selection: batch_select_thresholds()
6. Save artifacts:
   - Models: artifacts/models/*.joblib
   - Manifest: model_manifest.json
   - Metadata: dataset_metadata.json (updated with training_dataset_hash)
7. Update Registry: ModelRecord (src/oligotox/registry.py)
```

**Cross-Validation Strategy**: `src/oligotox/integrations_v03/datasail_splits.py`
- Sequence-aware splitting (avoids data leakage)
- Chemistry-aware grouping
- Maintains endpoint prevalence balance

---

## 🧪 MECHANISTIC EVIDENCE INTEGRATION

### Geneformer Cell-Type Context
**Module**: `src/oligotox/integrations/geneformer_context.py`
- **Purpose**: Cell-type-specific transcriptome embedding
- **Current Status**: Curated proxy (placeholder)
- **Future (v2.0)**: HuggingFace Geneformer model integration
- **Output**: Context vector for tissue-specific predictions

### TLR/Innate Immune Activation
**Module**: `src/oligotox/integrations/tlr_immune.py`
- **Function**: `score_tlr_innate(sequence)`
- **Mechanisms Detected**:
  - CpG motifs (TLR9 agonists)
  - dsRNA mimics (TLR3)
  - ssRNA motifs (TLR7/8)
- **Output**: Immune activation risk score (0-1)

### Seed Match Off-Target Scoring
**Module**: `src/oligotox/integrations/seed_match.py`
- **Function**: `score_seed_matches(sequence)`
- **Algorithm**: miRNA seed region matching (2-8mer)
- **Database**: Human 3'UTR transcriptome
- **Output**: Off-target gene silencing potential

### GSEA Pathway Enrichment
**Module**: `src/oligotox/integrations/gsea_pathway.py`
- **Function**: `score_pathways(deg_table)`
- **Gene Sets**: `src/oligotox/omics/gene_sets.py` (GENE_SETS dictionary)
- **Pathways**: Inflammatory, apoptosis, complement, coagulation
- **Output**: Pathway enrichment scores (NES, FDR)

### Candidate Mechanistic Evidence Builder
**Module**: `src/oligotox/integrations/candidate_mechanism.py`
- **Function**: `build_candidate_mechanistic_evidence()`
- **Integrates**:
  1. Seed match analysis
  2. TLR immune scoring
  3. Pathway enrichment
- **Output**: Unified mechanistic report for UI

---

## 🔬 BENCHMARKING & VALIDATION

### 1. OligoGym Benchmark Integration
**Module**: `src/oligotox/benchmarks/oligogym.py`

**Class**: `OligoGymDataset`
- **Cache**: `./data/benchmarks/`
- **Methods**:
  - `load_task(task_name)` - Load specific benchmark
  - `generate_synthetic()` - Synthetic test data
- **Status Function**: `benchmark_status()` in `evaluation/benchmark_registry.py`

### 2. ChEMBL Drug Comparison
**Module**: `src/oligotox/integrations/chembl_benchmark.py`

**Function**: `benchmark_against_drug(predicted, drug_name)`
- **Input**: OligoTox predictions + approved drug name
- **Output**: `BenchmarkComparison` dataclass
- **Metrics**: Concordance with known clinical toxicity

**Drug Database** (13 entries):
- Mipomersen, Fomivirsen, Eteplirsen, Nusinersen, Inotersen, Patisiran, etc.
- Each with known BLACK BOX warnings, REMS, clinical AEs

### 3. Model Comparator
**Module**: `src/oligotox/models/comparator.py`
- Compare OligoTox predictions against external models
- Track prediction concordance
- Identify divergent cases for further investigation

---

## 📝 REPRODUCIBILITY & VERSION CONTROL

### Dataset Versioning
1. **Hash Tracking**: `training_dataset_hash` field in `ModelRecord`
2. **Metadata Preservation**: `dataset_metadata.json` stores:
   - `training_csv` path
   - Feature list (column names)
   - Endpoint definitions
   - Data source provenance
   - Build timestamp
3. **Audit Trail**: `data/processed/training_dataset_rebuilt_audit.csv`

### Model Versioning
1. **Timestamped Runs**: 
   - Format: `train_model_v4_YYYYMMDD_HHMMSS_<hash>/`
   - Location: `artifacts/models/`
2. **Manifest per Run**: `model_manifest.json`
   - Training params
   - Endpoint list
   - Accuracy metrics
   - Feature importance
3. **Model Index**: `models_index.json` (aggregates all runs)
4. **Legacy Archives**:
   - v3 models: `artifacts/models_v3/` (pkl format)
   - Backups: `_local_only/backups/retrain_*/models.bak/`

### Code Provenance
**Import Graph** (key dependencies):
```
oligotox_app.py
  └─ oligotox.inference.api.run_oligotox()
      ├─ oligotox.pipelines.score_candidate.score_candidate()
      │   └─ oligotox.models.loader.load_dataset_metadata()
      └─ oligotox.omics (OMICS_REGISTRY, load_omics_table, join_omics_tables)

src/oligotox/training/train.py
  ├─ oligotox.pipelines.build_training_dataset.run()
  │   ├─ oligotox.features.aso_rich_features.build_aso_rich_features()
  │   ├─ oligotox.features.structure_features.build_structure_features()
  │   └─ oligotox.features.embedding_features.build_embedding_features()
  ├─ oligotox.features.v03_bridge.build_v03_feature_df()
  ├─ oligotox.integrations_v03.datasail_splits (cross-validation)
  └─ oligotox.pipelines.thresholds.batch_select_thresholds()
```

---

## 🔍 QUALITY ASSURANCE

### 1. Data QC Pipeline
**Omics QC** (`src/oligotox/omics/deg_qc.py`):
- Expression threshold: log2FC > 0.5, padj < 0.05
- Replicate correlation > 0.8
- Batch effect detection

**Sequence Validation** (`src/oligotox/validators.py`):
- `validate_sequence_input()` - ATGCU alphabet only
- Length constraints (15-30 nt)
- `validate_mapping_file()` - Format checks

### 2. Feature Engineering QC
- **Missing Value Handling**: Imputation or exclusion
- **Outlier Detection**: Z-score > 3 flagged
- **Feature Correlation**: Multicollinearity checks
- **Audit File**: `data/processed/training_dataset_rebuilt_audit.csv`

### 3. Model Evaluation
- **Training Accuracy**: 99.9% average
- **Cross-Validation**: DataSAIL k-fold (sequence-aware)
- **Threshold Optimization**: Per-endpoint F1-score maximization
- **Calibration**: Probability calibration curves

### 4. Test Suite
- **Unit Tests**: `pytest src/oligotox/` (quick tests pass 3/3 in 0.82s)
- **Integration Tests**: End-to-end prediction pipeline
- **Regression Tests**: Model performance benchmarks

---

## 📚 ADDITIONAL DOCUMENTATION

### Key Files
1. `PLATFORM_COMPLETE.md` - Production deployment status
3. `DATA_AND_MODEL_PROVENANCE.md` - This document
4. `README.md` - User-facing documentation

### Data Preservation
- **Archives**: `_local_only/archive/`
- **Backups**: `_local_only/backups/`

### Contact
- **Maintainer**: Beshoy Armanios, PhD
- **Email**: [Redacted for privacy]
- **GitHub**: [Project repository]
- **Status**: Actively maintained, production-ready v1.0

---

**Last Updated**: June 26, 2026, 12:00 PM EDT
**Platform Version**: v1.0 Production Release
**Models Trained**: 10 endpoint models (13 total artifacts)
**Training Accuracy**: 99.9% average across all endpoints
**Dataset Size**: 1,576 sequences, 49 features, 13 endpoints

