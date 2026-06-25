# OligoTox v0.3 — External Model Integration Guide

Mechanism-aware toxicity prediction for RNA/oligonucleotide therapeutics.
Five external model integrations wired into a single feature pipeline.

---

## Quick Start

```bash
# 1. Core dependencies
pip install scikit-learn pandas numpy scipy pyyaml requests

# 2. ViennaRNA (strongly recommended for D2/D3 structure features)
conda install -c bioconda viennarna
# or: pip install ViennaRNA

# 3. Transformer embeddings (optional — GPU recommended)
pip install transformers torch

# 4. Leakage-minimising splits (optional — for publication-grade evaluation)
pip install datasail

# 5. Run the integration smoke test
cd oligotox_v03_integrations
python demo_v03.py

# 6. Run the full pipeline on curated data
PYTHONPATH=. python scripts/run_pipeline.py \
    --input data/curated/curated_oligo_tox_template.csv \
    --split scaffold --train
```

---

## Integration 1 — Nucleotide/RNA Foundation Model Embeddings

**File:** `integrations/sequence_embeddings.py`

Captures sequence-level context beyond GC content. Embeddings encode
evolutionary, structural, and functional information from millions of genomes.
Falls back to 4-mer frequency vectors when transformers not available.

| Model | Params | HuggingFace repo |
|---|---|---|
| Nucleotide Transformer | 500M | `InstaDeepAI/nucleotide-transformer-500m-human-ref` |
| RNAErnie | 86M | `multimolecule/rnaernie` |
| RiNALMo | 650M | `lbcb-sci/RiNALMo` |
| BiRNA-BERT | 117M | `multimolecule/rnamsm` |

**Ref:** Dalla-Torre et al., Nat Methods 2024; Shuai et al. 2023

---

## Integration 2 — ASO Chemistry Features (ASOptimizer-inspired)

**File:** `integrations/aso_optimizer.py`

14 chemistry features from sequence + modification annotation:
PS fraction, PS clustering score, LNA count/fraction/flank asymmetry,
gapmer structure (gap/wing length), sugar mod fractions (2-MOE/OMe/2-F),
and activity/off-target proxies.

**Ref:** Kang et al., NAR 2025 (ASOptimizer); DOI:10.1093/nar/gkaf349

---

## Integration 3 — Transcriptome-Aware Seed/Off-Target Burden

**File:** `integrations/offtarget_burden.py`

TargetScan 7-mer seed hierarchy: 8mer > 7mer-m8 > 7mer-A1 > 6mer.
Weighted burden score against full human 3'UTR transcriptome when FASTA
available; sequence-complexity estimate otherwise.

```bash
# Download real human 3'UTR FASTA (~350MB)
bash scripts/download_transcriptome.sh
# Then set use_real_transcriptome: true in configs/integration_config.yaml
```

**Ref:** Agarwal et al., eLife 2015 (TargetScan context++ model)

---

## Integration 4 — ViennaRNA + TLR/RIG-I Structure Signals

**File:** `integrations/structure_signals.py`

12 features: MFE (kcal/mol), MFE per nucleotide, longest stem length,
dsRNA flag (≥6nt stem), RIG-I risk flag (≥21bp dsRNA), mean base-pair
probability, TLR7/8/9 motif densities, CpG density, poly-U score, GU
wobble count, self-complementarity score.

```bash
conda install -c bioconda viennarna
```

**Ref:** Lorenz et al., Algorithms Mol Biol 2011; Hornung et al. 2006

---

## Integration 5 — DataSAIL Leakage-Free Splits

**File:** `integrations/datasail_splits.py`

Three strategies in order of rigor:

| Strategy | When | Install |
|---|---|---|
| scaffold (seed-based) | Default | sklearn only |
| group (target gene / drug) | When target_gene column present | sklearn only |
| datasail (cluster-based) | Before publication | `pip install datasail` |

**Ref:** Gorgulla et al., Bioinformatics 2024; PMC 2025

---

## Endpoint Labeling Thresholds

| Endpoint | Threshold → label=1 | Primary source |
|---|---|---|
| hepatocellular_injury | ALT >3× ULN + histopath ≥2 | Dieckmann 2018; FDA 168708 |
| offtarget_gene_silencing | caspase >300% OR seed burden >10 | Papargyri 2019 |
| thrombocytopenia | platelet nadir <75% baseline | Inotersen NDA |
| complement_activation | C3a >2× OR SC5b-9 >2× baseline | Patisiran NDA; Janas 2018 |
| innate_immune_activation | ≥2 cytokines >5× OR TLR reporter >3× | Vollmer 2004 |
| renal_tubular_injury | KIM-1 >2× + histopath ≥2 | Wood 2020 |
| cns_neurotoxicity | NfL >3× baseline OR histopath ≥2 | Nusinersen NDA |
| coagulation_effect | PT/aPTT >1.5× ULN | PS backbone literature |
| lnp_immune_effect | ALT >3× + CRP >5× + cytokines | Patisiran/siRNA-LNP NDA |

---

## Data Downloads

```bash
cd oligotox_v03_integrations
mkdir -p data/raw

# 1. LNA gapmer hepatotoxicity (Dieckmann 2018)
curl -L "https://pmc.ncbi.nlm.nih.gov/articles/PMC5725219/pdf" \
  -o data/raw/dieckmann_2018_gapmer_tox.pdf

# 2. GalNAc siRNA nonclinical platform (Janas 2018)
curl -L "https://pmc.ncbi.nlm.nih.gov/articles/PMC6249674/pdf" \
  -o data/raw/janas_2018_galnac_nonclinical.pdf

# 3. FDA three ASOs + 18 impurities
curl -L "https://www.fda.gov/media/168708/download" \
  -o data/raw/fda_three_aso_impurities.pdf

# 4. Holistic siRNA toxicity review
curl -L "https://pmc.ncbi.nlm.nih.gov/articles/PMC10543595/pdf" \
  -o data/raw/siRNA_holistic_tox_review.pdf

# 5. Human 3'UTR transcriptome
bash scripts/download_transcriptome.sh
```

Manual downloads:
- Papargyri/Abzu ASO CSV: https://docs.abzu.ai/docs/tutorials/r/aso
- siRNAmod bulk export: http://crdd.osdd.net/servers/sirnamod/
- NIH OligoTox Challenge (Phase 2, closes Dec 31 2026): ncatsoligotox@mail.nih.gov

---

## File Structure

```
oligotox_v03_integrations/
├── INTEGRATION_GUIDE.md
├── demo_v03.py                          ← Smoke test: python demo_v03.py
├── requirements.txt
├── configs/
│   └── integration_config.yaml
├── integrations/
│   ├── __init__.py
│   ├── sequence_embeddings.py           ← Integration 1: NT / RNAErnie / RiNALMo
│   ├── aso_optimizer.py                 ← Integration 2: ASOptimizer chemistry (14 features)
│   ├── offtarget_burden.py              ← Integration 3: TargetScan seed burden
│   ├── structure_signals.py             ← Integration 4: ViennaRNA + TLR motifs (12 features)
│   ├── datasail_splits.py               ← Integration 5: leakage-free splits + leakage checker
│   └── feature_matrix.py               ← Assembles all 5 into one feature row
├── scripts/
│   ├── run_pipeline.py                  ← Full train/eval pipeline
│   └── download_transcriptome.sh       ← Ensembl 3'UTR FASTA download
└── data/
    ├── curated/
    │   └── curated_oligo_tox_template.csv  ← 4 seed rows (patisiran, inotersen, inclisiran, nusinersen)
    ├── raw/                             ← Downloaded PDFs and FASTA
    └── features/                        ← Generated feature matrices
```
