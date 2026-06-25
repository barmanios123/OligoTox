# OligoTox Training Dataset v2.0
*Curated June 2026 | 18 compounds | 5 modality classes | 13 endpoints*

## Files
- `processed/curated_oligo_tox_master_v2.csv` — 18-row master training table, 43 columns
- `processed/endpoint_summary_v2.csv` — per-endpoint coverage stats
- `processed/molecule_catalog_v2.json` — machine-readable molecule catalog with sequences
- `sequences/oligo_sequences_v2.fasta` — all sequences in FASTA format
- `docs/labeling_scheme_v2.csv` — 9-endpoint coding rules with thresholds and assays
- `docs/curation_priority_v2.csv` — 12-item next-curation queue
- `scripts/oligotox_download_v2.sh` — download script for PDFs and labels

## Modality coverage
| Modality | N | Key drugs |
|---|---|---|
| ASO_gapmer | 8 | inotersen, mipomersen, vupanorsen, LNA_HiTox_A/B, LNA_MidTox_C, LNA_LoTox_D/E |
| PMO | 4 | eteplirsen, golodirsen, viltolarsen, casimersen |
| siRNA | 4 | inclisiran, givosiran, lumasiran, patisiran |
| SSO | 1 | nusinersen |
| mRNA_LNP | 1 | BNT162b2 |

## Sparse endpoints (priority for curation v3)
- thrombocytopenia: 2/18 positive (11%) — add nusinersen detail, eplontersen
- complement_activation: 2/18 positive (11%) — add mRNA-1273, Alnylam ESC GalNAc
- cns_neurotoxicity: 1/18 positive (6%) — add intrathecal oligos, anti-tau ASO
- reproductive_dev_tox: 1/18 positive (6%) — add patisiran DART details, givosiran DART

## Sources
All data from: FDA NDA reviews, published drug labels, PMC11022857, Papargyri 2020,
Janas 2018, Dieckmann 2018, NEJM clinical trial publications.
See `source_doi` column for per-row primary reference.
