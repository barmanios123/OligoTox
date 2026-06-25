#!/usr/bin/env bash
# OligoTox v0.3 Dataset Downloader
# Usage: cd ~/oligotox && bash oligotox_datasets/scripts/download_datasets_v3.sh
set -e
ROOT="$HOME/oligotox/data/raw"
mkdir -p "$ROOT"
echo "=== OligoTox Dataset Downloader v0.3 ==="

dl() {
  local label="$1" url="$2" out="$3"
  [ -f "$ROOT/$out" ] && { echo "  [SKIP] $out"; return; }
  echo "  [GET] $label"
  curl -L --retry 3 --max-time 120 -A "Mozilla/5.0" -o "$ROOT/$out" "$url" 2>/dev/null \
    && echo "    ✓ $out" || echo "    ✗ FAILED — manual download required"
}

echo "--- BLOCK 1: LNA Gapmer Hepatotoxicity ---"
dl "Dieckmann 2018 MTNA" \
   "https://www.sciencedirect.com/science/article/pii/S2162253117302883/pdfft" \
   "dieckmann_2018_gapmer_tox.pdf"
dl "Burel 2015 RNaseH1 NAR PMC" \
   "https://pmc.ncbi.nlm.nih.gov/articles/PMC4344497/pdf/gkv117.pdf" \
   "burel_2015_rnaseH1_gapmer.pdf"
dl "Hagedorn 2018 nucleobase mods PMC" \
   "https://pmc.ncbi.nlm.nih.gov/articles/PMC6237733/pdf/gky752.pdf" \
   "hagedorn_2018_nucleobase.pdf"
echo "  MANUAL: Papargyri 2020 Table S1 (768 gapmers)"
echo "    URL: https://doi.org/10.1016/j.omtn.2019.12.011"
echo "    Save to: $ROOT/papargyri_2020_768_gapmers.csv"

echo "--- BLOCK 2: GalNAc siRNA Safety ---"
dl "Janas 2018 GalNAc nonclinical PMC" \
   "https://pmc.ncbi.nlm.nih.gov/articles/PMC6249674/pdf/10.1177_0192623318793117.pdf" \
   "janas_2018_galnac_nonclinical.pdf"
dl "Alnylam Capella tox slides 2018" \
   "https://capella.alnylam.com/wp-content/uploads/2018/05/ABS_Capella-Slides_FINAL05092018.pdf" \
   "alnylam_capella_galnac_tox.pdf"
dl "Prakash 2014 ESC GalNAc PMC" \
   "https://pmc.ncbi.nlm.nih.gov/articles/PMC4150748/pdf/gku484.pdf" \
   "prakash_2014_esc_galnac.pdf"

echo "--- BLOCK 3: FDA ASO Data ---"
dl "FDA 3-ASO + 18 impurities" \
   "https://www.fda.gov/media/168708/download" \
   "fda_three_aso_18_impurities.pdf"
echo "  MANUAL NDA labels:"
echo "    Inotersen:  https://www.accessdata.fda.gov/drugsatfda_docs/nda/2018/211172Orig1s000TOC.cfm"
echo "    Nusinersen: https://www.accessdata.fda.gov/drugsatfda_docs/nda/2016/209531Orig1s000TOC.cfm"
echo "    Patisiran:  https://www.accessdata.fda.gov/drugsatfda_docs/nda/2018/210922Orig1s000TOC.cfm"

echo "--- BLOCK 4: siRNA Off-Target & Immune ---"
dl "Holistic siRNA tox review 2023 PMC" \
   "https://pmc.ncbi.nlm.nih.gov/articles/PMC10543595/pdf" \
   "siRNA_tox_holistic_review_2023.pdf"
dl "siRNAmod ML paper 2024" \
   "https://www.sciencedirect.com/science/article/pii/S0888754324000363/pdfft" \
   "siRNAmod_ml_2024.pdf"
echo "  MANUAL: siRNAmod bulk export → http://crdd.osdd.net/servers/sirnamod/"

echo "--- BLOCK 5: CNS Safety ---"
dl "Holm 2022 CNS ASO biomarkers PMC" \
   "https://pmc.ncbi.nlm.nih.gov/articles/PMC9492605/pdf/s41392-022-01065-2.pdf" \
   "holm_2022_cns_aso_biomarkers.pdf"

echo "--- BLOCK 6: Mechanistic Reviews ---"
dl "Sewing 2017 PS-platelet PMC" \
   "https://pmc.ncbi.nlm.nih.gov/articles/PMC5363295/pdf/mt201617a.pdf" \
   "sewing_2017_ps_platelet.pdf"
dl "FDA 2023 oligo nonclinical guidance" \
   "https://www.fda.gov/media/167469/download" \
   "fda_2023_oligo_guidance.pdf"

echo "--- BLOCK 7: Terminated / Failed Drugs ---"
dl "Vupanorsen NEJM Evidence 2022" \
   "https://evidence.nejm.org/doi/pdf/10.1056/EVIDoa2200122" \
   "vupanorsen_2022.pdf"

echo "--- BLOCK 8: NIH OligoTox Challenge ---"
echo "  Phase 2 active through Dec 31, 2026"
echo "  Register: ncatsoligotox@mail.nih.gov"
echo "  Info: https://ncats.nih.gov/research/research-activities/oligotox"

echo ""
echo "=== Done. Files in $ROOT:"
ls -lh "$ROOT"/*.pdf 2>/dev/null | awk '{print "  "$5, $9}' || echo "  (none downloaded yet)"
