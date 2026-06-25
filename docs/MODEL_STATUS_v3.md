# OligoTox — Model Status & Limitations (v0.3)

## Endpoint coverage after v0.3 curation

| Endpoint | N+ | Rate | Status | Notes |
|---|---|---|---|---|
| hepatocellular_injury | 17 | 45% | **STRONG** | D2 hybridization + LNA fraction signal |
| hypersensitivity | 16 | 42% | **STRONG** | PS backbone + route feature |
| injection_site_reaction | 22 | 58% | STRONG | Near-trivial (SC PS = 1) |
| offtarget_gene_silencing | 8 | 21% | **MODERATE** | GC, LNA fraction, seed-match |
| renal_tubular_injury | 9 | 24% | **MODERATE** | PS/PMO backbone |
| innate_immune_activation | 7 | 18% | **MODERATE** | LNP, CpG motifs |
| renal_accumulation | 7 | 18% | **MODERATE** | PS/PMO class flag |
| cytokine_release | 6 | 16% | **MODERATE** | LNP delivery |
| thrombocytopenia | 4 | 10% | **THIN** | PS fraction + cyno species |
| coagulation_effect | 4 | 10% | **THIN** | PS fraction |
| complement_activation | 4 | 10% | **THIN** | LNP IV + PS IV |
| cns_neurotoxicity | 3 | 8% | **SPARSE → rules only** | IT route, high LNA |
| reproductive_dev_tox | 1 | 3% | **SPARSE → rules only** | LNP, human-specific seq |

## Correct interpretation for computational audiences

This is a **first-generation, mechanism-aware model**. A correct characterisation:

> "OligoTox v0.3 is a research-grade directional risk-ranking system with a
> mix of moderately trained and weakly trained endpoints. It is sufficient for
> hypothesis generation and candidate comparison. Endpoints with fewer than
> 4 training positives are mechanistic-rule-derived flags, not calibrated ML
> predictions. Confidence intervals for low-prevalence endpoints are wide and
> should not be treated as independent of the domain expert."

## Limitations to disclose

1. **Dataset size**: 38 compounds total; LNA gapmer positives dominated by caspase in vitro data — in vivo endpoints underrepresented.
2. **Batch effects**: mix of in vitro (caspase), animal (rat/cyno), and clinical data; no correction applied in v0.3.
3. **Class imbalance**: 5 endpoints at ≤10% positive rate; oversampling applied but cannot substitute for real data.
4. **Feature coverage**: no SMILES/molecular graph yet; v0.4 target is ChemProp/MPNN encoding of PS pattern and LNA positions.
5. **No external validation**: all metrics are cross-validated on training set; held-out test set requires NIH OligoTox Challenge or partner data.

## Priority data additions (v0.4 targets)

| Action | Impact |
|---|---|
| Add Papargyri 2020 full 768-gapmer table | +50-100 hepatocellular/offtarget positives → STRONG tier |
| Extract inotersen/eplontersen NDA platelet tables | thrombocytopenia THIN → MODERATE |
| Add patisiran cyno complement data (NDA Table 9) | complement_activation THIN → MODERATE |
| Add nusinersen/tofersen CNS histopath tables | cns_neurotoxicity SPARSE → THIN |
| Register for NIH OligoTox Challenge Phase 2 | Full external validation dataset on completion (Dec 2026) |

## Architecture summary

```
Input: Sequence + Chemistry + Delivery (30 features)
  ↓
Mechanistic Rules Engine (D1/D2/D3 three-domain framework)
  ↓
Per-endpoint Classifier
  STRONG/MODERATE:  Calibrated RF or GBM (class-balanced + oversampled)
  THIN:             Calibrated Logistic (C=0.1, high regularization)
  SPARSE/RULES_ONLY: Mechanistic rules + domain priors only
  ↓
Output: Risk score + Domain attribution + Confidence flag + Assay recommendation
```
