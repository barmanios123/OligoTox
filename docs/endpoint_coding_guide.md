# OligoTox Endpoint Coding Guide  
*For curators encoding NDA / published study findings into binary labels*

---

## General Rule

| Value | Use when |
|-------|----------|
| **1** | Label lists finding as: clinically meaningful adverse reaction; nonclinical finding with clear safety relevance; listed under Warnings and Precautions |
| **0** | Background lab variation without interpretation; isolated finding not interpreted as tox; not mentioned |

---

## Per-Endpoint Coding Rules (sourced from FDA labels)

### `endpoint_hepatocellular_injury`
Set to **1** if: ALT/AST ≥ 3× ULN as adverse reaction; histologic hepatocellular injury; listed in hepatic Warnings.  
Set to **0** if: "no hepatotoxicity signal" (e.g., PMO class — LiverTox score E for eteplirsen, golodirsen, viltolarsen, casimersen).

### `endpoint_renal_tubular_injury`
Set to **1** if: "proximal tubular injury/degeneration/necrosis" in nonclinical section OR drives a renal monitoring Warning.  
Set to **0** if: no renal histopathology; normal urinalysis/cystatin C in all studies.  
**PMO class note**: all 4 approved PMO drugs = 1 (animal data; mandatory monitoring label).

### `endpoint_renal_accumulation`
Set to **1** if: explicit renal accumulation or long renal tissue half-life with associated risk language.  
All naked PMO IV drugs = 1 (renal excretion pathway documented).

### `endpoint_thrombocytopenia`
Set to **1** if: platelet nadir < 75% baseline OR thrombocytopenia listed as adverse reaction.  
PMO drugs = 0 (no PS backbone; platelet-PF4 interaction requires PS).  
PS-ASO drugs (e.g., inotersen, mipomersen) = evaluate carefully; inotersen = 1 (Black Box Warning).

### `endpoint_coagulation_effect`
Set to **1** if: aPTT/PT prolongation listed as safety signal; anticoagulation monitoring required.

### `endpoint_complement_activation`
Set to **1** if: C3a/C5a/SC5b-9 elevation OR Warnings mention complement activation.  
PMO drugs = 0.  
LNP-siRNA (e.g., patisiran) = 1 (IRR / complement infusion reactions).

### `endpoint_innate_immune_activation`
Set to **1** if: TLR7/8/9, RIG-I, MDA-5 pathway mentioned; cytokine elevations (IFN-α, IL-6) as safety signal.

### `endpoint_cytokine_release`
Set to **1** if: cytokine storm, flu-like syndrome, or specific cytokine elevations described as adverse reaction.

### `endpoint_hypersensitivity`
Set to **1** if: hypersensitivity, anaphylaxis, or angioedema described as Warnings, Contraindications, or ≥ grade 2 AE.  
Golodirsen = 1 (anaphylaxis contraindication). Casimersen = 1 (anaphylaxis + angioedema Warnings).  
Viltolarsen = 1 (urticaria, Warning). Eteplirsen = 1 (rash/pyrexia Warning).

### `endpoint_cns_neurotoxicity`
Set to **1** if: neurotoxicity, delayed CNS findings, NfL elevation, demyelination, or intrathecal safety concerns described.

### `endpoint_offtarget_gene_silencing`
Set to **1** if: off-target transcript knockdown; caspase > 300% threshold (Dieckmann/Papargyri standard); seed-match tox.

### `endpoint_reproductive_dev_tox`
Set to **1** if: DART findings; teratogenicity; fetal weight reduction; developmental toxicity.

### `endpoint_injection_site_reaction`
Set to **1** if: injection/infusion site reaction listed as ≥ 10% frequency AE.

---

## Severity Coding

| Level | Numeric | Use when |
|-------|---------|----------|
| `none` | 0 | No signal in that organ system |
| `low` | 1 | Mild/infrequent; managed by routine monitoring |
| `medium` | 2 | Clear signal; specific monitoring required; notable frequency |
| `high` | 3 | Boxed warning; serious clinical events; strong restrictions |

---

## PMO Drug Reference (ground truth, v1.0)

| Drug | Renal | Immune | BBW | Renal Monitor | Hypersens Warn |
|------|-------|--------|-----|---------------|----------------|
| eteplirsen (Exondys 51) | medium | low | No | Yes | Yes |
| golodirsen (Vyondys 53) | medium | medium | No | Yes | Yes (contraind) |
| viltolarsen (Viltepso) | medium | low | No | Yes | No |
| casimersen (Amondys 45) | medium | medium | No | Yes | Yes |

**Primary sources**: FDA labels 206488s035, 211970s002, 213026; LiverTox NBK555333, NBK588131; NDA PharmR 206488Orig1s000

---

## Next Steps for Curation

1. Open each source PDF listed in `primary_source` field of `pmo_label_training_rows.csv`  
2. Apply the rules above to fill in any remaining drugs (ASO PS gapmers, GalNAc siRNAs, LNP-siRNAs)  
3. Run: `python scripts/curation/build_pmo_label_dataset.py`  
4. Run: `python scripts/curation/merge_into_master.py`  
5. Run: `python scripts/curated_csv_to_training.py` → `python train.py`
