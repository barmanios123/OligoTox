#!/usr/bin/env bash
# =============================================================================
# OligoTox Dataset Curation Download Script
# Run from: cd ~/Downloads/oligotox && bash ~/Downloads/oligotox_download.sh
# All 200-status verified URLs as of June 2026
# =============================================================================

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DATA="$SCRIPT_DIR/data/raw"
mkdir -p "$DATA"/{lna_gapmer,galnac_sirna,aso_tox,sirna_modified,reviews_labeling,regulatory,github_data}

UA="Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 Safari/537.36"
DL() { echo "  → Downloading: $2"; curl -sSL -A "$UA" --max-time 60 --retry 3 -o "$1" "$3" && echo "    ✓ $(du -sh "$1" | cut -f1)  $2" || echo "    ✗ FAILED (may be gated): $2"; }

# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "============================================================"
echo " BLOCK 1: LNA GAPMER — HYBRIDIZATION-DEPENDENT TOX"
echo " D2 endpoint: offtarget_gene_silencing, hepatocellular_injury"
echo "============================================================"
# Papargyri 2019 — 768 LNA gapmers, caspase + hepatotox readouts
# THE core dataset; caspase >300% = toxic (threshold from paper)
DL "$DATA/lna_gapmer/papargyri_2019_lna_768gapmers.pdf" \
   "Papargyri 2019: 768 LNA gapmers (caspase/hepatotox) [PMC6965521]" \
   "https://pmc.ncbi.nlm.nih.gov/articles/PMC6965521/pdf"

# Dieckmann 2018 — sensitive in vitro hybridization-dep tox assay
DL "$DATA/lna_gapmer/dieckmann_2018_hybridization_dep_tox.pdf" \
   "Dieckmann 2018: hybridization-dep hepatotox, LNA gapmers [PMC5725219]" \
   "https://pmc.ncbi.nlm.nih.gov/articles/PMC5725219/pdf"

# Burel 2015 — RNaseH1-mediated off-target hepatotox mechanism
DL "$DATA/lna_gapmer/burel_2015_rnaseH1_hepatotox.pdf" \
   "Burel 2015: RNaseH1 off-target hepatotox mechanism [Semantic Scholar]" \
   "https://pdfs.semanticscholar.org/d365/b53b86511b9d18d43a1958680ea5174e0f09.pdf"

# Nucleobase mods reducing hepatotox 2022
DL "$DATA/lna_gapmer/nucleobase_mods_reduced_hepatotox_2022.pdf" \
   "2022: Nucleobase modifications reduce LNA gapmer hepatotox [PMC9303313]" \
   "https://pmc.ncbi.nlm.nih.gov/articles/PMC9303313/pdf"

# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "============================================================"
echo " BLOCK 2: GalNAc siRNA NONCLINICAL SAFETY"
echo " D3 endpoints: complement, cytokine, renal, hepatic"
echo "============================================================"
# Janas 2018 — ESC GalNAc siRNA nonclinical platform tox (rat+monkey)
DL "$DATA/galnac_sirna/janas_2018_galnac_nonclinical_safety.pdf" \
   "Janas 2018: GalNAc siRNA nonclinical safety profile [PMC6249674]" \
   "https://pmc.ncbi.nlm.nih.gov/articles/PMC6249674/pdf"

# GalNAc follow-up study 2021
DL "$DATA/galnac_sirna/galnac_sirna_safety_followup_2021.pdf" \
   "GalNAc siRNA safety follow-up 2021 [PMC7710393]" \
   "https://pmc.ncbi.nlm.nih.gov/articles/PMC7710393/pdf"

# Alnylam Capella slide deck — platform nonclinical data
DL "$DATA/galnac_sirna/alnylam_capella_galnac_tox_slides.pdf" \
   "Alnylam Capella GalNAc siRNA tox slides [capella.alnylam.com]" \
   "https://capella.alnylam.com/wp-content/uploads/2018/05/ABS_Capella-Slides_FINAL05092018.pdf" || \
   echo "    ⚠ Capella slides may require browser download from: https://capella.alnylam.com/wp-content/uploads/2018/05/ABS_Capella-Slides_FINAL05092018.pdf"

# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "============================================================"
echo " BLOCK 3: FDA ASO TOXICITY DATA"
echo " All 3 endpoints + impurity tox data"
echo "============================================================"
# FDA mipomersen/inotersen/nusinersen + 18 impurities — primary human hepatocytes
DL "$DATA/aso_tox/fda_three_aso_18impurities_hepatocyte_tox.pdf" \
   "FDA: 3 ASO drugs + 18 impurities, human hepatocyte tox [media/168708]" \
   "https://www.fda.gov/media/168708/download"

# FDA 2024 draft guidance — nonclinical ASO safety framework
DL "$DATA/regulatory/fda_2024_guidance_oligo_nonclinical_safety.pdf" \
   "FDA 2024 Draft Guidance: Nonclinical Safety Assessment of Oligo Therapeutics" \
   "https://www.fda.gov/media/183540/download" || \
   echo "    ⚠ Try: https://www.fda.gov/regulatory-information/search-fda-guidance-documents/nonclinical-safety-assessment-oligonucleotide-based-therapeutics"

# SOT safety webinar — ASO platform class tox (Henry)
DL "$DATA/aso_tox/sot_henry_aso_platform_safety.pdf" \
   "Henry/SOT: ASO platform class safety data (2'-MOE) [toxicology.org]" \
   "https://www.toxicology.org/education/docs/CCTWebinarApr15_Henry.pdf"

# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "============================================================"
echo " BLOCK 4: siRNAmod + CHEMICALLY MODIFIED siRNA"
echo " D2 endpoint: offtarget silencing, innate immune"
echo "============================================================"
# siRNAmod ML paper supplemental (ScienceDirect - may be gated)
echo "  → Note: siRNAmod supplemental XLSX requires ScienceDirect login."
echo "    Manual download from: https://www.sciencedirect.com/science/article/pii/S0888754324000363"
echo "    Then save as: $DATA/sirna_modified/siRNAmod_ml_supplement.xlsx"
echo "  → siRNAmod DB bulk download: http://crdd.osdd.net/servers/sirnamod/dd.php"

# siRNA off-target/innate review 2022
DL "$DATA/sirna_modified/sirna_offtarget_innate_review_2022.pdf" \
   "siRNA off-target + innate immune review 2022 [PMC8604658]" \
   "https://pmc.ncbi.nlm.nih.gov/articles/PMC8604658/pdf"

# 2025 siRNA translational safety review
DL "$DATA/sirna_modified/sirna_translational_clinical_2025.pdf" \
   "siRNA translational + clinical safety 2025 [PMC12477600]" \
   "https://pmc.ncbi.nlm.nih.gov/articles/PMC12477600/pdf"

# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "============================================================"
echo " BLOCK 5: MECHANISTIC REVIEWS — LABEL DEFINITIONS"
echo " Use for mapping assay readouts → binary 0/1 endpoint labels"
echo "============================================================"
# Holistic siRNA tox mechanisms (2023) — maps to all 13 endpoints
DL "$DATA/reviews_labeling/siRNA_holistic_tox_mechanisms_2023.pdf" \
   "2023: Holistic siRNA tox mechanisms (all delivery/chem routes) [PMC10543595]" \
   "https://pmc.ncbi.nlm.nih.gov/articles/PMC10543595/pdf"

# 2023 oligo tox landscape review
DL "$DATA/reviews_labeling/oligo_tox_landscape_2023.pdf" \
   "2023: Oligonucleotide tox landscape review [PMC9908799]" \
   "https://pmc.ncbi.nlm.nih.gov/articles/PMC9908799/pdf"

# PS complement/platelet mechanism (foundational)
DL "$DATA/reviews_labeling/ps_complement_platelet_mechanism.pdf" \
   "PS-ASO complement + platelet mechanism (foundational) [PMC6491344]" \
   "https://pmc.ncbi.nlm.nih.gov/articles/PMC6491344/pdf"

# NCBI Bookshelf — Preclinical Safety Assessment chapter
echo "  → NCBI Bookshelf (HTML only): https://www.ncbi.nlm.nih.gov/books/NBK584239/"
echo "    Save page as PDF manually and place in: $DATA/reviews_labeling/"

# NIH OligoTox Challenge Phase 1 webinar (Jan 2026) — endpoint framework
DL "$DATA/reviews_labeling/nih_oligotox_challenge_2026_webinar.pdf" \
   "NIH NCATS OligoTox Challenge 2026 Phase 1 Webinar [ncats.nih.gov]" \
   "https://ncats.nih.gov/sites/default/files/2026-02/NCATS-OligoTox-Phase1-Webinar-2026-01-09-508.pdf"

# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "============================================================"
echo " BLOCK 6: GitHub REPOS — Code + Possible Companion Data"
echo "============================================================"
GH="$DATA/github_data"
echo "  → Cloning ClinASO (gapmer design platform, may include data)..."
git clone --quiet https://github.com/wusuw/ClinASO "$GH/ClinASO" 2>/dev/null || \
  echo "    ✓ ClinASO already cloned or skipped"

echo "  → Cloning openASO (ASO efficacy datasets)..."
git clone --quiet https://github.com/lackeylela/openASO "$GH/openASO" 2>/dev/null || \
  echo "    ✓ openASO already cloned or skipped"

# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "============================================================"
echo " BLOCK 7: ADDITIONAL APPROVED DRUG LABELS (NDA/BLA data)"
echo " Use with chembl_benchmark.py for ground-truth endpoint labels"
echo "============================================================"
echo "  → FDA Drugs@FDA links (manual browser download):"
echo "    Mipomersen (Kynamro): https://www.accessdata.fda.gov/scripts/cder/daf/index.cfm?event=BasicSearch.process&ApplNo=203568"
echo "    Inotersen (Tegsedi):  https://www.accessdata.fda.gov/scripts/cder/daf/index.cfm?event=BasicSearch.process&ApplNo=211945"
echo "    Nusinersen (Spinraza):https://www.accessdata.fda.gov/scripts/cder/daf/index.cfm?event=BasicSearch.process&ApplNo=209531"
echo "    Inclisiran (Leqvio):  https://www.accessdata.fda.gov/scripts/cder/daf/index.cfm?event=BasicSearch.process&ApplNo=214012"
echo "    Patisiran (Onpattro): https://www.accessdata.fda.gov/scripts/cder/daf/index.cfm?event=BasicSearch.process&ApplNo=210922"
echo "    Givosiran (Givlaari): https://www.accessdata.fda.gov/scripts/cder/daf/index.cfm?event=BasicSearch.process&ApplNo=212194"

# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "============================================================"
echo " BLOCK 8: NIH OligoTox Challenge (ACTIVE — Phase 2 until Dec 2026)"
echo "============================================================"
echo "  → Challenge page: https://www.nih.gov/challenges/oligonucleotide-toxicity-oligotox-open-data-challenge"
echo "  → Contact: ncatsoligotox@mail.nih.gov"
echo "  → Phase 2 runs May 1 – Dec 31, 2026. Winning datasets become public."
echo "  → You may want to SUBMIT to this challenge with OligoTox v0.3 data."

# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "============================================================"
echo " DOWNLOAD COMPLETE — Summary"
echo "============================================================"
echo ""
echo "  Files downloaded to: $DATA"
echo ""
find "$DATA" -type f | sort | while read f; do
  size=$(du -sh "$f" 2>/dev/null | cut -f1)
  echo "    $size  $f"
done

echo ""
echo "  MANUAL STEPS STILL NEEDED:"
echo "  1. siRNAmod XLSX:  http://crdd.osdd.net/servers/sirnamod/dd.php"
echo "  2. Abzu ASO CSV:   https://docs.abzu.ai/docs/tutorials/r/aso → follow data links"
echo "  3. FDA DAF labels: Links above for 6 approved drugs (PDF download)"
echo "  4. NCBI Bookshelf: https://www.ncbi.nlm.nih.gov/books/NBK584239/ → Print to PDF"
echo "  5. NIH Challenge:  Register at https://www.nih.gov/challenges/oligonucleotide-toxicity-oligotox-open-data-challenge"
