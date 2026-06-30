#!/usr/bin/env bash
# ============================================================
# OligoTox Platform — Full Dataset Download Script
# Verified URLs as of June 2026
# Run from: ~/oligotox/
# Usage:  cd ~/oligotox && bash oligotox_dataset_download.sh
# ============================================================

set -euo pipefail
LOGFILE="data/download_log.txt"
FAILFILE="data/manual_download_required.txt"

mkdir -p data/raw/tier1_labelable \
         data/raw/tier2_sequence_chemistry \
         data/raw/tier3_mechanistic_omics \
         data/raw/github_repos \
         data/raw/regulatory

log()  { echo "[$(date '+%H:%M:%S')] $*" | tee -a "$LOGFILE"; }
fail() { echo "[MANUAL NEEDED] $*" | tee -a "$FAILFILE"; }

log "=== OligoTox Dataset Download Started ==="
log "Destination: $(pwd)/data/raw/"

# ============================================================
# TIER 1 — Directly labelable tox datasets
# These populate feature_matrix.csv + label_matrix.csv NOW
# ============================================================
log ""
log "--- TIER 1: Directly Labelable Tox Datasets ---"

# 1a. FDA: Toxicity of 3 ASO drugs + 18 impurities (primary human hepatocytes)
#    Endpoints: hepatocellular_injury, renal_tubular_injury
#    Source: FDA Science Forum / NCTR — direct PDF
log "Downloading: FDA 3-ASO + 18 impurities tox PDF..."
curl -L --retry 3 --max-time 60 \
  "https://www.fda.gov/media/168708/download" \
  -o "data/raw/tier1_labelable/fda_three_aso_impurities.pdf" \
  && log "  OK: fda_three_aso_impurities.pdf" \
  || fail "  fda_three_aso_impurities.pdf — try: https://www.fda.gov/media/168708/download"

# 1b. Janas et al. 2018 — GalNAc siRNA nonclinical safety (PMC full text)
#    Endpoints: complement_activation, cytokine_release, hepatocellular_injury, renal_tubular_injury
#    DOI: 10.1177/0192623318793297 | PMC: PMC6249674
log "Downloading: Janas 2018 GalNAc siRNA nonclinical safety (PMC PDF)..."
curl -L --retry 3 --max-time 60 \
  "https://pmc.ncbi.nlm.nih.gov/articles/PMC6249674/pdf/10.1177_0192623318793297.pdf" \
  -o "data/raw/tier1_labelable/janas_2018_galnac_nonclinical_safety.pdf" \
  && log "  OK: janas_2018_galnac_nonclinical_safety.pdf" \
  || { log "  Trying alternate PMC PDF URL..."; \
  curl -L --retry 3 --max-time 60 \
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC6249674/pdf" \
    -o "data/raw/tier1_labelable/janas_2018_galnac_nonclinical_safety.pdf" \
    && log "  OK (alt): janas_2018_galnac_nonclinical_safety.pdf" \
    || fail "  janas_2018_galnac_nonclinical_safety.pdf — visit: https://pmc.ncbi.nlm.nih.gov/articles/PMC6249674/"; }

# 1c. Dieckmann et al. 2018 — Hybridization-dependent gapmer tox (caspase/hepatotox)
#    PMC: PMC5725219 | Endpoints: offtarget_gene_silencing, hepatocellular_injury
log "Downloading: Dieckmann 2018 hybridization-dependent gapmer tox (PMC)..."
curl -L --retry 3 --max-time 60 \
  "https://pmc.ncbi.nlm.nih.gov/articles/PMC5725219/pdf/10.1016_j.omtn.2017.09.013.pdf" \
  -o "data/raw/tier1_labelable/dieckmann_2018_hybridization_gapmer_tox.pdf" \
  && log "  OK: dieckmann_2018_hybridization_gapmer_tox.pdf" \
  || { curl -L --retry 3 --max-time 60 \
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC5725219/pdf" \
    -o "data/raw/tier1_labelable/dieckmann_2018_hybridization_gapmer_tox.pdf" \
    && log "  OK (alt): dieckmann_2018_hybridization_gapmer_tox.pdf" \
    || fail "  dieckmann_2018_hybridization_gapmer_tox.pdf — visit: https://pmc.ncbi.nlm.nih.gov/articles/PMC5725219/"; }

# 1d. Burel et al. 2016 — RNaseH1-mediated promiscuous gapmer hepatotoxicity
#    PMC: PMC4797265 | Endpoints: hepatocellular_injury, offtarget_gene_silencing
log "Downloading: Burel 2016 RNaseH1 gapmer hepatotox (PMC)..."
curl -L --retry 3 --max-time 60 \
  "https://pmc.ncbi.nlm.nih.gov/articles/PMC4797265/pdf" \
  -o "data/raw/tier1_labelable/burel_2016_rnash1_gapmer_hepatotox.pdf" \
  && log "  OK: burel_2016_rnash1_gapmer_hepatotox.pdf" \
  || fail "  burel_2016_rnash1_gapmer_hepatotox.pdf — visit: https://pmc.ncbi.nlm.nih.gov/articles/PMC4797265/"

# 1e. scASO — single-cell ASO activity in mouse brain (Minikel/Mortberg 2023)
#    GitHub repo with analytical dataset CSV + R code
#    Endpoints: cns_delayed_neurotoxicity (D1: dose/cell-type knockdown profiles)
log "Cloning: scASO single-cell ASO brain dataset (ericminikel/scaso)..."
if [ ! -d "data/raw/github_repos/scaso" ]; then
  git clone --depth 1 https://github.com/ericminikel/scaso \
    "data/raw/github_repos/scaso" \
    && log "  OK: scaso cloned" \
    || fail "  scaso — visit: https://github.com/ericminikel/scaso"
else
  log "  SKIP: scaso already cloned"
fi

# 1f. NHP deep brain ASO dataset (Frei/Minikel 2025) — NHP intrathecal cell-type data
log "Cloning: NHP deep brain ASO dataset (nhp_deepbrain_aso)..."
if [ ! -d "data/raw/github_repos/nhp_deepbrain_aso" ]; then
  git clone --depth 1 https://github.com/ericminikel/nhp_deepbrain_aso \
    "data/raw/github_repos/nhp_deepbrain_aso" \
    && log "  OK: nhp_deepbrain_aso cloned" \
    || fail "  nhp_deepbrain_aso — visit: https://github.com/ericminikel/nhp_deepbrain_aso"
else
  log "  SKIP: nhp_deepbrain_aso already cloned"
fi

# 1g. UniTox — 2,418 FDA-approved drugs with GPT-4o toxicity ratings (zenodo)
#    Use to anchor hepatotoxicity, nephrotoxicity, hematotox calibration for approved oligos
#    Zenodo: https://zenodo.org/records/14042913
log "Downloading: UniTox 2024 dataset (zenodo CSV)..."
curl -L --retry 3 --max-time 120 \
  "https://zenodo.org/records/14042913/files/UniTox.csv" \
  -o "data/raw/tier1_labelable/unitox_2418_drugs.csv" \
  && log "  OK: unitox_2418_drugs.csv" \
  || fail "  unitox_2418_drugs.csv — visit: https://zenodo.org/records/14042913"

# 1h. DILIrank/DILIst — FDA NCTR liver tox benchmark (primary label reference)
#    1036 drugs, 4-class hepatotox labels — use for hepatocellular_injury calibration
log "Downloading: FDA LTKB DILIst benchmark dataset..."
curl -L --retry 3 --max-time 60 \
  "https://www.fda.gov/media/99399/download" \
  -o "data/raw/tier1_labelable/dilirank_fda_ltkb.csv" \
  && log "  OK: dilirank_fda_ltkb.csv" \
  || { log "  Trying DILIst alternate..."; \
  curl -L --retry 3 --max-time 60 \
    "https://www.fda.gov/media/134581/download" \
    -o "data/raw/tier1_labelable/dilirank_fda_ltkb.csv" \
    && log "  OK (alt): dilirank_fda_ltkb.csv" \
    || fail "  dilirank — visit: https://www.fda.gov/science-research/liver-toxicity-knowledge-base-ltkb/ltkb-benchmark-dataset  (download 'LTKB-BD Excel file')"; }

# 1i. SIDER 4.1 — side effects for 1,430 approved drugs (adverse event ground truth)
#    Use: populate benchmark.py with hematologic, renal, hepatic, immune AE labels
log "Downloading: SIDER 4.1 drug adverse effects (meddra_all_se.tsv.gz)..."
curl -L --retry 3 --max-time 120 \
  "http://sideeffects.embl.de/media/download/meddra_all_se.tsv.gz" \
  -o "data/raw/tier1_labelable/sider4_meddra_all_se.tsv.gz" \
  && log "  OK: sider4_meddra_all_se.tsv.gz" \
  || fail "  sider4_meddra_all_se.tsv.gz — visit: http://sideeffects.embl.de/download/"

log "Downloading: SIDER 4.1 drug names..."
curl -L --retry 3 --max-time 60 \
  "http://sideeffects.embl.de/media/download/drug_names.tsv" \
  -o "data/raw/tier1_labelable/sider4_drug_names.tsv" \
  && log "  OK: sider4_drug_names.tsv" \
  || fail "  sider4_drug_names.tsv — visit: http://sideeffects.embl.de/download/"

log "Downloading: SIDER 4.1 side effect frequencies..."
curl -L --retry 3 --max-time 120 \
  "http://sideeffects.embl.de/media/download/meddra_freq.tsv.gz" \
  -o "data/raw/tier1_labelable/sider4_meddra_freq.tsv.gz" \
  && log "  OK: sider4_meddra_freq.tsv.gz" \
  || fail "  sider4_meddra_freq.tsv.gz — visit: http://sideeffects.embl.de/download/"


# ============================================================
# TIER 2 — Sequence/chemistry datasets with partial endpoint coverage
# ============================================================
log ""
log "--- TIER 2: Sequence/Chemistry Datasets ---"

# 2a. siRNAmod — chemically modified siRNA database (all entries, CSV format)
#    Endpoints: innate_immune_activation, offtarget_gene_silencing, offtarget_mirna_loading
log "Downloading: siRNAmod database (CSV)..."
curl -L --retry 3 --max-time 60 \
  "http://crdd.osdd.net/servers/sirnamod/dd_excel.php" \
  -o "data/raw/tier2_sequence_chemistry/sirnamod_all_entries.xls" \
  && log "  OK: sirnamod_all_entries.xls" \
  || { curl -L --retry 3 --max-time 60 \
    "http://crdd.osdd.net/servers/sirnamod/sirnamod_download.csv" \
    -o "data/raw/tier2_sequence_chemistry/sirnamod_all_entries.csv" \
    && log "  OK (alt): sirnamod_all_entries.csv" \
    || fail "  sirnamod — visit: http://crdd.osdd.net/servers/sirnamod/dd.php — click 'Download data in CSV format'"; }

# 2b. ClinASO — computational-experimental gapmer ASO platform (GitHub)
#    Contains ASO efficacy + design rule data; use for D1 on-target feature training
log "Cloning: ClinASO gapmer platform (wusuw/ClinASO)..."
if [ ! -d "data/raw/github_repos/ClinASO" ]; then
  git clone --depth 1 https://github.com/wusuw/ClinASO \
    "data/raw/github_repos/ClinASO" \
    && log "  OK: ClinASO cloned" \
    || fail "  ClinASO — visit: https://github.com/wusuw/ClinASO"
else
  log "  SKIP: ClinASO already cloned"
fi

# 2c. openASO — regulatory region targeting for ASOs (GitHub)
log "Cloning: openASO (lackeylela/openASO)..."
if [ ! -d "data/raw/github_repos/openASO" ]; then
  git clone --depth 1 https://github.com/lackeylela/openASO \
    "data/raw/github_repos/openASO" \
    && log "  OK: openASO cloned" \
    || fail "  openASO — visit: https://github.com/lackeylela/openASO"
else
  log "  SKIP: openASO already cloned"
fi

# 2d. TOXsiRNA 2026 preprint + supplementary (biorxiv)
#    882 chemically modified siRNAs with toxicity scores 0-1
#    Endpoints: offtarget_gene_silencing, innate_immune_activation
log "Downloading: TOXsiRNA 2026 preprint PDF (biorxiv)..."
curl -L --retry 3 --max-time 60 \
  "https://www.biorxiv.org/content/10.64898/2026.02.12.705521v1.full.pdf" \
  -o "data/raw/tier2_sequence_chemistry/toxsirna_2026_preprint.pdf" \
  && log "  OK: toxsirna_2026_preprint.pdf" \
  || fail "  toxsirna_2026_preprint.pdf — visit: https://www.biorxiv.org/content/10.64898/2026.02.12.705521v1"

# 2d-ii. TOXsiRNA supplemental data (check biorxiv for .xlsx/.csv attachment)
log "Downloading: TOXsiRNA 2026 supplemental data (Excel)..."
curl -L --retry 3 --max-time 60 \
  "https://www.biorxiv.org/content/10.64898/2026.02.12.705521v1.supplemental" \
  -o "data/raw/tier2_sequence_chemistry/toxsirna_2026_supplemental.html" \
  && log "  OK: toxsirna_2026_supplemental.html (parse for attachment URLs)" \
  || fail "  toxsirna_supp — visit: https://www.biorxiv.org/content/10.64898/2026.02.12.705521v1 and click Supplementary Material"

# 2e. Papargyri et al. 2020 — 768 ASO caspase/tox dataset (PMC6965521)
#    Gold-standard D2 hybridization-dependent tox; caspase >300% = tox label
log "Downloading: Papargyri 2020 ASO tox design rules (PMC PDF)..."
curl -L --retry 3 --max-time 60 \
  "https://pmc.ncbi.nlm.nih.gov/articles/PMC6965521/pdf" \
  -o "data/raw/tier2_sequence_chemistry/papargyri_2020_aso_design_rules.pdf" \
  && log "  OK: papargyri_2020_aso_design_rules.pdf" \
  || fail "  papargyri_2020 — visit: https://pmc.ncbi.nlm.nih.gov/articles/PMC6965521/ and download PDF + supplementary"

# 2f. siRNA ML paper — "Using machine learning to predict chemically modified siRNA"
#    ScienceDirect 2024; check supplementary for CSV
log "Downloading: siRNAmod ML paper 2024 (ScienceDirect - check for supp CSV)..."
curl -L --retry 3 --max-time 60 \
  "https://www.sciencedirect.com/science/article/pii/S0888754324000363" \
  -o "data/raw/tier2_sequence_chemistry/sirnamod_ml_2024_landing.html" \
  && log "  OK: sirnamod_ml_2024_landing.html (parse for supplementary data link)" \
  || fail "  sirnamod_ml — visit: https://www.sciencedirect.com/science/article/pii/S0888754324000363 (download Supplementary data Excel)"

# 2g. Benchmarking LLMs for ASO efficacy (biorxiv 2026) — sequence/efficacy features
log "Downloading: LLM ASO efficacy benchmark 2026 (biorxiv PDF)..."
curl -L --retry 3 --max-time 60 \
  "https://www.biorxiv.org/content/10.64898/2026.02.17.706455v1.full.pdf" \
  -o "data/raw/tier2_sequence_chemistry/llm_aso_efficacy_benchmark_2026.pdf" \
  && log "  OK: llm_aso_efficacy_benchmark_2026.pdf" \
  || fail "  llm_aso_benchmark — visit: https://www.biorxiv.org/content/10.64898/2026.02.17.706455v1"

# 2h. Alnylam Capella GalNAc platform slides (2018 PDF — safety profile data)
log "Downloading: Alnylam Capella GalNAc platform slides PDF..."
curl -L --retry 3 --max-time 60 \
  "https://capella.alnylam.com/wp-content/uploads/2018/05/ABS_Capella-Slides_FINAL05092018.pdf" \
  -o "data/raw/tier2_sequence_chemistry/alnylam_capella_galnac_platform_slides_2018.pdf" \
  && log "  OK: alnylam_capella_galnac_platform_slides_2018.pdf" \
  || fail "  alnylam_capella — visit: https://capella.alnylam.com/ and search for platform safety data"

# 2i. Toxicity synergy paper 2025 — chemistry × sequence interplay
#    PubMed: 41078095
log "Downloading: ASO tox chemistry×sequence interplay 2025 (PubMed landing)..."
curl -L --retry 3 --max-time 60 \
  "https://pubmed.ncbi.nlm.nih.gov/41078095/" \
  -o "data/raw/tier2_sequence_chemistry/aso_tox_chemistry_sequence_2025_pubmed.html" \
  && log "  OK: aso_tox_chemistry_sequence_2025_pubmed.html" \
  || fail "  aso_tox_chem_seq — visit: https://pubmed.ncbi.nlm.nih.gov/41078095/"

# 2j. Hybridization-dependent off-target risk assessment (SAGE 2024)
log "Downloading: Hybridization off-target risk assessment 2024 (SAGE journal page)..."
curl -L --retry 3 --max-time 60 \
  "https://journals.sagepub.com/doi/10.1089/nat.2024.0072" \
  -o "data/raw/tier2_sequence_chemistry/hybridization_offtarget_risk_2024_sage.html" \
  && log "  OK: hybridization_offtarget_risk_2024_sage.html" \
  || fail "  hybridization_offtarget — visit: https://journals.sagepub.com/doi/10.1089/nat.2024.0072"


# ============================================================
# TIER 3 — Mechanistic/omics datasets (v0.4 calibration)
# ============================================================
log ""
log "--- TIER 3: Mechanistic/Omics Datasets ---"

# 3a. Holistic siRNA tox review 2023 (PMC10543595) — label logic + D1/D2/D3 framework
log "Downloading: Holistic siRNA tox review 2023 (PMC PDF)..."
curl -L --retry 3 --max-time 60 \
  "https://pmc.ncbi.nlm.nih.gov/articles/PMC10543595/pdf" \
  -o "data/raw/tier3_mechanistic_omics/holistic_sirna_tox_review_2023.pdf" \
  && log "  OK: holistic_sirna_tox_review_2023.pdf" \
  || fail "  holistic_sirna_review — visit: https://pmc.ncbi.nlm.nih.gov/articles/PMC10543595/"

# 3b. Toxicological landscape of short nucleic acid therapeutics (2024 review)
log "Downloading: Toxicological landscape nucleic acid therapeutics review (insights.bio)..."
curl -L --retry 3 --max-time 60 \
  "https://www.insights.bio/nucleic-acid-insights/journal/article/3629/the-toxicological-landscape-of-short-nucleic-acid-therapeutics-a-short-review-of-toxicity-mechanisms-clinical-translation-and-regulatory-frameworks" \
  -o "data/raw/tier3_mechanistic_omics/tox_landscape_nucleic_acid_review_2024.html" \
  && log "  OK: tox_landscape_nucleic_acid_review_2024.html" \
  || fail "  tox_landscape — visit: https://www.insights.bio/nucleic-acid-insights/journal/article/3629/"

# 3c. NCBI Bookshelf — Preclinical Safety Assessment of Therapeutic Oligonucleotides
log "Downloading: Preclinical safety assessment of oligos (NCBI Bookshelf)..."
curl -L --retry 3 --max-time 60 \
  "https://www.ncbi.nlm.nih.gov/books/NBK584239/" \
  -o "data/raw/tier3_mechanistic_omics/preclinical_safety_oligos_ncbi_bookshelf.html" \
  && log "  OK: preclinical_safety_oligos_ncbi_bookshelf.html" \
  || fail "  ncbi_bookshelf_oligos — visit: https://www.ncbi.nlm.nih.gov/books/NBK584239/"

# 3d. HepatoToxicity Portal (HTP, KOBIC 2025) — DILI database with GNN predictions
log "Downloading: HepatoToxicity Portal landing (HTP, kobic.re.kr)..."
curl -L --retry 3 --max-time 60 \
  "https://kobic.re.kr/htp/" \
  -o "data/raw/tier3_mechanistic_omics/htp_hepatotoxicity_portal_landing.html" \
  && log "  OK: htp_hepatotoxicity_portal_landing.html (navigate to download section)" \
  || fail "  htp_hepatotox — visit: https://kobic.re.kr/htp/ — look for Download/Export section"

# 3e. DILIsym / srijitseal DILI prediction dataset (GitHub) — DILI featurization
log "Cloning: srijitseal/DILI prediction dataset (GitHub)..."
if [ ! -d "data/raw/github_repos/DILI_srijitseal" ]; then
  git clone --depth 1 https://github.com/srijitseal/DILI \
    "data/raw/github_repos/DILI_srijitseal" \
    && log "  OK: DILI_srijitseal cloned" \
    || fail "  DILI_srijitseal — visit: https://github.com/srijitseal/DILI"
else
  log "  SKIP: DILI_srijitseal already cloned"
fi

# 3f. UniTox small molecule benchmark (zenodo) — SMILES-level tox ratings
log "Downloading: UniTox small-molecule benchmark CSV (zenodo)..."
curl -L --retry 3 --max-time 60 \
  "https://zenodo.org/records/14042913/files/UniTox-Small-Molecule-Benchmark.csv" \
  -o "data/raw/tier3_mechanistic_omics/unitox_small_molecule_benchmark.csv" \
  && log "  OK: unitox_small_molecule_benchmark.csv" \
  || fail "  unitox_smb — visit: https://zenodo.org/records/14042913"

# 3g. FDA regulatory guidance — 2024 nonclinical safety of antisense oligos
log "Downloading: FDA draft guidance oligonucleotide nonclinical safety 2024..."
curl -L --retry 3 --max-time 60 \
  "https://www.fda.gov/media/179461/download" \
  -o "data/raw/regulatory/fda_2024_draft_guidance_oligo_nonclinical_safety.pdf" \
  && log "  OK: fda_2024_draft_guidance_oligo_nonclinical_safety.pdf" \
  || { curl -L --retry 3 --max-time 60 \
    "https://www.fda.gov/regulatory-information/search-fda-guidance-documents/nonclinical-safety-assessment-investigational-oligonucleotide-therapeutics" \
    -o "data/raw/regulatory/fda_2024_oligo_safety_guidance_landing.html" \
    && log "  OK (landing): fda_2024_oligo_safety_guidance_landing.html" \
    || fail "  fda_guidance_oligos — visit: https://www.fda.gov and search 'oligonucleotide nonclinical safety guidance 2024'"; }

# 3h. PS backbone + complement mechanism key reference (Vollmer/Krieg)
log "Downloading: Vollmer 2004 PS complement activation (PubMed landing)..."
curl -L --retry 3 --max-time 60 \
  "https://pubmed.ncbi.nlm.nih.gov/15107527/" \
  -o "data/raw/tier3_mechanistic_omics/vollmer_2004_ps_complement_pubmed.html" \
  && log "  OK: vollmer_2004_ps_complement_pubmed.html" \
  || fail "  vollmer_2004 — visit: https://pubmed.ncbi.nlm.nih.gov/15107527/"

# 3i. ChEMBL oligonucleotide drug data (API — live download of nucleic acid drugs)
log "Downloading: ChEMBL oligonucleotide drugs (API JSON — nucleic acid biotherapeutics)..."
curl -L --retry 3 --max-time 90 \
  "https://www.ebi.ac.uk/chembl/api/data/molecule.json?molecule_type=Oligonucleotide&max_phase=4&limit=100" \
  -o "data/raw/regulatory/chembl_approved_oligonucleotides.json" \
  && log "  OK: chembl_approved_oligonucleotides.json" \
  || fail "  chembl_oligos — visit: https://www.ebi.ac.uk/chembl/explore/compounds/ and filter molecule_type=Oligonucleotide"

# Also pull the broader biotherapeutics (covers ASOs, siRNAs, mRNA)
curl -L --retry 3 --max-time 90 \
  "https://www.ebi.ac.uk/chembl/api/data/biotherapeutic.json?limit=500" \
  -o "data/raw/regulatory/chembl_biotherapeutics.json" \
  && log "  OK: chembl_biotherapeutics.json" \
  || fail "  chembl_biotherapeutics — API: https://www.ebi.ac.uk/chembl/api/data/biotherapeutic.json"


# ============================================================
# TIER 1b — GEO datasets (key transcriptomic studies)
# Use: fasterq-dump or GEOparse; listed as accessions for the Python getter below
# ============================================================
log ""
log "--- GEO Accession List (fetch with oligotox_geo_fetch.py) ---"
cat > "data/raw/geo_accessions_to_fetch.txt" << 'GEOEOF'
# GEO accessions for OligoTox training data
# Format: ACCESSION  DESCRIPTION  PRIMARY_ENDPOINT
#
# ASO / siRNA off-target transcriptomics
GSE27638   ASO_hepatic_fibrosis_EGFR   hepatocellular_injury
GSE47334   siRNA_off_target_seed_match_transcriptomics   offtarget_gene_silencing
GSE68580   siRNA_innate_immune_RISC_loading_HEK293   innate_immune_activation
#
# Complement / immune activation
GSE152641  innate_immune_two_gene_marker_panel   innate_immune_activation
#
# LNP / delivery system transcriptomics
GSE171110  LNP_lipid_nanoparticle_immune_response_mouse_liver   complement_activation
#
# SEARCH THESE MANUALLY at https://www.ncbi.nlm.nih.gov/gds
# Query: "antisense oligonucleotide toxicity transcriptomics"
# Query: "siRNA off-target hepatotoxicity gene expression"
# Query: "LNP innate immune activation liver RNA-seq"
GEOEOF
log "  OK: geo_accessions_to_fetch.txt written"

log ""
log "=== Download Phase Complete ==="
log "Log: $LOGFILE"
log "Manual items: $FAILFILE"

# Print summary counts
DOWNLOADED=$(find data/raw -name "*.pdf" -o -name "*.csv" -o -name "*.tsv.gz" -o -name "*.json" -o -name "*.xls" 2>/dev/null | wc -l)
CLONED=$(find data/raw/github_repos -maxdepth 1 -mindepth 1 -type d 2>/dev/null | wc -l)
log ""
log "Files downloaded: $DOWNLOADED"
log "Git repos cloned: $CLONED"
log ""
log "Next steps:"
log "  1. Run: python oligotox/scripts/geo_fetch.py  (fetch GEO accessions above)"
log "  2. Run: python oligotox/scripts/sider_parse.py  (parse SIDER TSVs for oligo drugs)"
log "  3. Run: python oligotox/scripts/chembl_oligo_parse.py  (extract approved oligo profiles)"
log "  4. Manually extract tables from PDFs → add rows to data/raw/curated_oligo_tox_template.csv"
log "  5. Run: python oligotox/scripts/curated_csv_to_training.py && python oligotox/train.py"
