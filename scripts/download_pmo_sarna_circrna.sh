#!/usr/bin/env bash
set -euo pipefail

mkdir -p data/raw/tier1_labelable data/raw/tier2_sequence_chemistry data/raw/tier3_mechanistic_omics

# PMO approved drugs — FDA NDA safety summaries
curl -L "https://www.accessdata.fda.gov/drugsatfda_docs/label/2016/206488s000lbl.pdf" \
  -o data/raw/tier1_labelable/eteplirsen_nda_label.pdf

curl -L "https://www.accessdata.fda.gov/drugsatfda_docs/label/2019/211970s000lbl.pdf" \
  -o data/raw/tier1_labelable/golodirsen_nda_label.pdf

curl -L "https://www.accessdata.fda.gov/drugsatfda_docs/label/2020/213026s000lbl.pdf" \
  -o data/raw/tier1_labelable/viltolarsen_nda_label.pdf

curl -L "https://www.accessdata.fda.gov/drugsatfda_docs/label/2021/213569s000lbl.pdf" \
  -o data/raw/tier1_labelable/casimersen_nda_label.pdf

# saRNA safety — PMC8759012 (LNP-nCoVsaRNA Phase 1)
curl -L "https://pmc.ncbi.nlm.nih.gov/articles/PMC8759012/pdf" \
  -o data/raw/tier2_sequence_chemistry/sarba_lnp_phase1_safety_pmc8759012.pdf

# circRNA therapeutics review 2026
curl -L "https://pmc.ncbi.nlm.nih.gov/articles/PMC13095804/pdf" \
  -o data/raw/tier3_mechanistic_omics/circrna_therapeutics_review_2026.pdf
