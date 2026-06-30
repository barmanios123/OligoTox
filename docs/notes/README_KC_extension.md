# Kupffer Cell Activation — 30+ Labeled Examples Extension

Generated: June 29, 2026
Target: Improve `kupffer_activation` model from AUROC 0.750 (17 examples) → target AUROC > 0.85

---

## Problem

The initial training run had only 17 Kupffer cell activation examples — insufficient for
a Random Forest with StratifiedKFold CV. This package adds 33 positive cases + 5 negative
controls, bringing the total to 50 labeled sequences for this endpoint.

---

## Files

| File | Purpose |
|------|---------|
| `kupffer_activation_30plus_examples.csv` | 38 new rows (33 positive, 5 negative) |
| `load_kc_examples.py` | Merge script; produces training_dataset_v3.csv |

---

## Evidence Breakdown

| Tier | n | Sources |
|------|---|---------|
| **T1** (primary in vivo) | 26 | Jackson 2021 (PMC7856291); Janas 2018 (PMC6249674 Alnylam); Henry 2002 (PMID11991849); Patisiran/Lumasiran/Mipomersen/Inotersen FDA labels; CALAA-01 Phase 1; PD-L1 LNP (PMID23423360) |
| **T2** (review/class effects) | 11 | Frazier 2015 ToxicolPathol; SOT 2010 oligonucleotide class effects; DOAJ PS-TLR9 activation; ScienceDirect 2025 LNP immune review |
| **T3** (mechanism inference) | 1 | MASLD-primed KC context (your HNF1AOS1 DEG data) |

---

## Key Mechanistic Insights Captured

### 1. LNP → KC-PAF axis (Jackson 2021 T1)
- Kupffer cells preferentially engulf siRNA nanocarriers; endosome escape triggers PAF release
- Mechanism is **endosome-escape dependent** — zwitterionic micelles without endosomolytic core do NOT activate KCs
- MC3 (Patisiran-matched) LNP confirmed lethal via KC-PAF in tumor-bearing mice
- **Feature implication**: `is_lnp=1` + `is_endosomolytic` = highest KC risk

### 2. GalNAc-siRNA KC drug accumulation (Janas 2018 T1)
- All 6 Alnylam ESC GalNAc-siRNAs show KC basophilic granules at ≥100x pharmacologic dose in monkeys
- **NOT adverse** as an isolated finding; becomes adverse when accompanied by ALT elevation or inflammatory infiltration
- **Feature implication**: `is_galnac=1` + `dose_multiplier>100x` = KC granule finding (non-adverse tier)

### 3. First-gen PS ASO + TLR9 (T2)
- PS backbone alone activates TLR9 cooperatively; CpG motifs amplify the effect
- `has_cpg` + `has_ps` = KC activation risk regardless of modality
- **Feature implication**: `has_cpg=1` + `has_ps=1` is the rule-based signature

### 4. Approved drug clinical KC activation (T1)
- Mipomersen (black box hepatotoxicity) — clinical KC + hepatocellular vacuolation at therapeutic dose
- Inotersen — thrombocytopenia via KC-complement-platelet mechanism
- CALAA-01 Phase 1 discontinued in part due to infusion reactions consistent with KC-PAF axis

---

## Integration Steps

```bash
# 1. Copy files
cp kupffer_activation_30plus_examples.csv data/processed/
cp load_kc_examples.py .

# 2. Merge into training_dataset_v3.csv
python load_kc_examples.py \
    --kc_csv   data/processed/kupffer_activation_30plus_examples.csv \
    --training data/processed/training_dataset_v2.csv \
    --out      data/processed/training_dataset_v3.csv

# 3. Retrain both hepatic endpoints
python train_hepatic_endpoints.py \
    --dataset data/processed/training_dataset_v3.csv \
    --outdir  models/
```

---

## Expected AUROC Improvement

| Model | Before | Expected After |
|-------|--------|---------------|
| `hepatic_vacuolation_rf` | 0.999 (stable, 602 examples) | Unchanged |
| `kupffer_activation_rf`  | 0.750 (17 examples) | **>0.85** (50+ examples) |

The improvement depends on class balance in your full training_dataset_v2.csv.
If the negative class is still >>5:1, add `class_weight='balanced'` in the RF
(already coded in train_hepatic_endpoints.py automatically).

---

## Sources

1. Jackson et al. (2021) *Biomaterials* 268:120528 — PMC7856291
   KC-PAF axis as dose-limiting toxicity of siRNA nanocarriers; generalizable to LNP class

2. Janas et al. (2018) *Toxicol Pathol* 46(7):735-745 — PMC6249674
   Nonclinical safety of 6 ESC GalNAc-siRNAs; KC basophilic granules at supratherapeutic doses

3. Henry et al. (2002) *PubMed* 11991849
   Targeting KCs with ASOs; KC obliteration prevents liver damage

4. Frazier (2015) *Toxicol Pathol* 43:78-89
   ASO therapy toxicity review; KC class effects; complement; PS backbone immunostimulation

5. SOT 2010 Oligonucleotide Therapeutics toxicological considerations
   KC vacuolization class effect; TLR9/CpG signaling

6. DOAJ PS-TLR9 cooperative activation (Characterization of cooperative PS-oligo activation of TLR9)
   PS-oligos activate TLR9 independent of CpG in cooperative manner
