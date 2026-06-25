# Biology validity check

## Overview

- Curated v2 rows: 18
- Curated v2 columns: 43
- Training label rows: 3259
- Training label columns: 13

## Endpoint coverage vs training prevalence

| Endpoint | Curated positives | Curated rate | Coverage call | Training positives | Training rate | Verdict |
|---|---:|---:|---|---:|---:|---|
| cns_neurotoxicity | 1 | 0.056 | Sparse |  |  | missing_from_training_matrix |
| cns_toxicity |  |  |  | 1016 | 0.312 | directionally_consistent |
| coagulation_effect | 2 | 0.111 | Sparse | 4 | 0.001 | expected_sparse |
| complement_activation | 2 | 0.111 | Sparse | 6 | 0.002 | expected_sparse |
| cytokine_release | 4 | 0.222 | Thin | 8 | 0.002 | warning_thin_endpoint_very_sparse_in_training |
| hepatic_vacuolation |  |  |  | 6 | 0.002 | directionally_consistent |
| hepatocellular_injury | 8 | 0.444 | Good | 341 | 0.105 | directionally_consistent |
| hypersensitivity | 10 | 0.556 | Good |  |  | missing_from_training_matrix |
| injection_site_reaction | 12 | 0.667 | Good | 17 | 0.005 | directionally_consistent |
| innate_immune_activation | 4 | 0.222 | Thin | 7 | 0.002 | warning_thin_endpoint_very_sparse_in_training |
| kupffer_cell_activation |  |  |  | 9 | 0.003 | directionally_consistent |
| offtarget_gene_silencing | 4 | 0.222 | Thin | 339 | 0.104 | directionally_consistent |
| renal_accumulation | 7 | 0.389 | Thin | 3 | 0.001 | warning_thin_endpoint_very_sparse_in_training |
| renal_tubular_injury | 8 | 0.444 | Good | 3 | 0.001 | warning_good_endpoint_looks_sparse_in_training |
| reproductive_dev_tox | 1 | 0.056 | Sparse |  |  | missing_from_training_matrix |
| thrombocytopenia | 2 | 0.111 | Sparse | 4 | 0.001 | expected_sparse |

## PMO class sanity check

PMO rows found by name match: 4

| Endpoint | Mean label | Positive count | N rows |
|---|---:|---:|---:|
| renal_tubular_injury | 1.000 | 4 | 4 |
| renal_accumulation | 1.000 | 4 | 4 |
| hypersensitivity | 1.000 | 4 | 4 |
| thrombocytopenia | 0.000 | 0 | 4 |
| complement_activation | 0.000 | 0 | 4 |
| innate_immune_activation | 0.000 | 0 | 4 |
| hepatocellular_injury | 0.000 | 0 | 4 |

Expected pattern: renal_tubular_injury / renal_accumulation / hypersensitivity should be elevated relative to thrombocytopenia / complement_activation / innate_immune_activation for PMO-class molecules.

## Interpretation rules

- Good curated endpoints should usually not collapse to near-zero prevalence in the full training label matrix.
- Sparse curated endpoints are expected to remain sparse in training and should be interpreted cautiously.
- PMO-class rows should preserve the backbone logic: renal and hypersensitivity signals can persist while PS-driven platelet/complement endpoints remain low.
- If a biologically strong curated endpoint looks nearly absent in training, inspect the builder pipeline for dropped or renamed columns.

## Next actions

- Review any row marked `warning_good_endpoint_looks_sparse_in_training`.
- Review PMO rows if thrombocytopenia or complement_activation are unexpectedly high.
- Add more rows for sparse endpoints such as thrombocytopenia, complement_activation, innate_immune_activation, and reproductive/CNS toxicity.
- Use the Streamlit app for candidate-level sniff tests after this file-level consistency check.
