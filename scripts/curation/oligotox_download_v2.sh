#!/usr/bin/env bash
# oligotox_download_v2.sh - Run from ~/oligotox/
set -e; mkdir -p data/raw/pdfs data/raw/labels

echo "=== FDA ASO+impurities (opens in browser if blocked) ==="
curl -L "https://www.fda.gov/media/168708/download" -o data/raw/pdfs/fda_3aso_18impurities.pdf

echo "=== Janas 2018 GalNAc nonclinical safety (PMC) ==="
curl -L "https://pmc.ncbi.nlm.nih.gov/articles/PMC6249674/pdf" -o data/raw/pdfs/janas_2018_GalNAc.pdf

echo "=== siRNA holistic tox review 2023 (PMC) ==="
curl -L "https://pmc.ncbi.nlm.nih.gov/articles/PMC10543595/pdf" -o data/raw/pdfs/siRNA_holistic_tox_2023.pdf

echo "=== ASO adverse reactions 2024 (PMC) ==="
curl -L "https://pmc.ncbi.nlm.nih.gov/articles/PMC11022857/pdf" -o data/raw/pdfs/aso_adv_reactions_2024.pdf

echo "=== Burel 2015 RNaseH1 hepatotox (PMC) ==="
curl -L "https://pmc.ncbi.nlm.nih.gov/articles/PMC4005641/pdf" -o data/raw/pdfs/burel_2015_RNaseH1.pdf

echo "=== FDA Labels (browser auth required - manual download recommended) ==="
echo "Inotersen: https://www.accessdata.fda.gov/drugsatfda_docs/label/2019/211172s004lbl.pdf"
echo "Mipomersen: https://www.accessdata.fda.gov/drugsatfda_docs/label/2019/203568s011lbl.pdf"
echo "Eteplirsen: https://www.accessdata.fda.gov/drugsatfda_docs/label/2020/206488s016lbl.pdf"
echo "Patisiran: https://www.accessdata.fda.gov/drugsatfda_docs/label/2023/210922s011lbl.pdf"
echo "Inclisiran: https://www.accessdata.fda.gov/drugsatfda_docs/label/2022/214012s007lbl.pdf"
echo "Givosiran: https://www.accessdata.fda.gov/drugsatfda_docs/label/2021/212194s003lbl.pdf"

echo "=== MANUAL: siRNAmod export ==="
echo "Visit: http://crdd.osdd.net/servers/sirnamod/dd.php"
echo "Export all → save to data/raw/siRNAmod_export.xlsx"

echo "=== MANUAL: Papargyri aso.csv ==="
echo "Visit: https://docs.abzu.ai/docs/tutorials/r/aso"
echo "Follow data link → save to data/raw/papargyri_aso_768gapmers.csv"

echo "Download script complete."
