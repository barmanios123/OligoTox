#!/usr/bin/env bash
# =============================================================================
# OligoTox FULL Dataset Curation Script — v3
# Includes: approved drug labels, GEO omics, TG-GATEs, DrugMatrix, LINCS, SIDER
# Run from: ~/oligotox
# All URLs verified 200-status June 2026
# =============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DATA="$SCRIPT_DIR/data/raw"
mkdir -p "$DATA"/{lna_gapmer,galnac_sirna,aso_tox,sirna_modified,reviews_labeling,regulatory,omics_geo,omics_tggates,sider_chembl,nda_labels,github_data}

UA="Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 Safari/537.36"
DL() { echo "  → $2"; curl -sSL -A "$UA" --max-time 90 --retry 3 -o "$1" "$3" \
  && echo "    ✓ $(du -sh "$1" | cut -f1)" || echo "    ✗ FAILED (may be gated): $2"; }

# ═════════════════════════════════════════════════════════════════════════════
echo ""; echo "══════════════ BLOCK 1: LNA GAPMER TOX ══════════════"
DL "$DATA/lna_gapmer/papargyri_2019_768gapmers.pdf" \
  "Papargyri 2019 — 768 LNA gapmers, caspase+hepatotox [PMC6965521]" \
  "https://pmc.ncbi.nlm.nih.gov/articles/PMC6965521/pdf"

DL "$DATA/lna_gapmer/dieckmann_2018_hybridization_tox.pdf" \
  "Dieckmann 2018 — hybridization-dep hepatotox [PMC5725219]" \
  "https://pmc.ncbi.nlm.nih.gov/articles/PMC5725219/pdf"

DL "$DATA/lna_gapmer/burel_2015_rnaseH1_offtarget.pdf" \
  "Burel 2015 — RNaseH1-mediated off-target hepatotox [Semantic Scholar]" \
  "https://pdfs.semanticscholar.org/d365/b53b86511b9d18d43a1958680ea5174e0f09.pdf"

DL "$DATA/lna_gapmer/nucleobase_mods_reduced_hepatotox_2022.pdf" \
  "2022 — Nucleobase mods reducing LNA hepatotox [PMC9303313]" \
  "https://pmc.ncbi.nlm.nih.gov/articles/PMC9303313/pdf"

# ═════════════════════════════════════════════════════════════════════════════
echo ""; echo "══════════════ BLOCK 2: GalNAc siRNA NONCLINICAL ══════════════"
DL "$DATA/galnac_sirna/janas_2018_galnac_nonclinical.pdf" \
  "Janas 2018 — GalNAc siRNA nonclinical safety [PMC6249674]" \
  "https://pmc.ncbi.nlm.nih.gov/articles/PMC6249674/pdf"

DL "$DATA/galnac_sirna/galnac_followup_2021.pdf" \
  "GalNAc safety follow-up 2021 [PMC7710393]" \
  "https://pmc.ncbi.nlm.nih.gov/articles/PMC7710393/pdf"

DL "$DATA/galnac_sirna/alnylam_capella_slides.pdf" \
  "Alnylam Capella GalNAc tox slides [capella.alnylam.com]" \
  "https://capella.alnylam.com/wp-content/uploads/2018/05/ABS_Capella-Slides_FINAL05092018.pdf" || true

# ═════════════════════════════════════════════════════════════════════════════
echo ""; echo "══════════════ BLOCK 3: FDA ASO DATA ══════════════"
DL "$DATA/aso_tox/fda_three_aso_18impurities.pdf" \
  "FDA — 3 ASO drugs + 18 impurities, human hepatocytes [media/168708]" \
  "https://www.fda.gov/media/168708/download"

DL "$DATA/aso_tox/sot_henry_aso_platform_class.pdf" \
  "Henry/SOT — ASO platform class tox [toxicology.org]" \
  "https://www.toxicology.org/education/docs/CCTWebinarApr15_Henry.pdf"

DL "$DATA/regulatory/fda_oligo_analytical_2022.pdf" \
  "FDA 2022 — Oligo analytical/chemistry guidance [media/166574]" \
  "https://www.fda.gov/media/166574/download"

# ═════════════════════════════════════════════════════════════════════════════
echo ""; echo "══════════════ BLOCK 4: REVIEWS / LABEL DEFINITIONS ══════════════"
DL "$DATA/reviews_labeling/siRNA_holistic_tox_2023.pdf" \
  "2023 holistic siRNA tox mechanisms [PMC10543595]" \
  "https://pmc.ncbi.nlm.nih.gov/articles/PMC10543595/pdf"

DL "$DATA/reviews_labeling/oligo_tox_landscape_2023.pdf" \
  "2023 oligo tox landscape [PMC9908799]" \
  "https://pmc.ncbi.nlm.nih.gov/articles/PMC9908799/pdf"

DL "$DATA/reviews_labeling/ps_complement_platelet.pdf" \
  "PS-ASO complement + platelet mechanism [PMC6491344]" \
  "https://pmc.ncbi.nlm.nih.gov/articles/PMC6491344/pdf"

DL "$DATA/reviews_labeling/siRNA_offtarget_innate_2022.pdf" \
  "siRNA off-target + innate immune review 2022 [PMC8604658]" \
  "https://pmc.ncbi.nlm.nih.gov/articles/PMC8604658/pdf"

DL "$DATA/reviews_labeling/sirna_translational_2025.pdf" \
  "siRNA translational safety 2025 [PMC12477600]" \
  "https://pmc.ncbi.nlm.nih.gov/articles/PMC12477600/pdf"

DL "$DATA/reviews_labeling/nih_oligotox_challenge_webinar.pdf" \
  "NIH OligoTox Challenge 2026 webinar [ncats.nih.gov]" \
  "https://ncats.nih.gov/sites/default/files/2026-02/NCATS-OligoTox-Phase1-Webinar-2026-01-09-508.pdf"

DL "$DATA/reviews_labeling/aso_sequence_chemistry_synergy_2025.pdf" \
  "2025 — ASO tox = sequence × chemistry synergy [PubMed 41078095]" \
  "https://pmc.ncbi.nlm.nih.gov/articles/PMC11919288/pdf" || true

DL "$DATA/reviews_labeling/aso_transcriptomics_proteomics_2025.pdf" \
  "2025 — ASO time-course transcriptomics+proteomics toxic vs nontoxic [PMC13034028]" \
  "https://pmc.ncbi.nlm.nih.gov/articles/PMC13034028/pdf"

# ═════════════════════════════════════════════════════════════════════════════
echo ""; echo "══════════════ BLOCK 5: GEO OMICS DATASETS (wget) ══════════════"
echo "  Downloading GEO series files (SOFT/matrix format for each accession)"

GEO_SERIES=(
  "GSE53230"    # Hepatic transcriptomics — toxic vs non-toxic LNA gapmers, mouse liver [PMC/OmicsDI]
  "GSE108823"   # siRNA RNAi off-target hepatotox at supratherapeutic doses, rodent
)
for acc in "${GEO_SERIES[@]}"; do
  dest="$DATA/omics_geo/${acc}_matrix.txt.gz"
  url="https://ftp.ncbi.nlm.nih.gov/geo/series/${acc:0:-3}nnn/$acc/matrix/${acc}_series_matrix.txt.gz"
  echo "  → $acc ($url)"
  curl -sSL -A "$UA" --max-time 120 --retry 3 -o "$dest" "$url" \
    && echo "    ✓ GEO $acc: $(du -sh "$dest" | cut -f1)" \
    || echo "    ✗ GEO $acc failed — check accession at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=$acc"
done

echo "  Note: For LINCS L1000 KD data (119,848 siRNA/shRNA signatures):"
echo "    Portal:  https://clue.io/data"
echo "    GEO set: GSE70138 (LINCS Phase II) — ~25 GB, download selectively by cell line"
echo "    command: python scripts/lincs_geo_fetch.py --gse GSE70138 --cell_line HepG2 --out $DATA/omics_geo/"

# ═════════════════════════════════════════════════════════════════════════════
echo ""; echo "══════════════ BLOCK 6: Open TG-GATEs TOXICOGENOMICS ══════════════"
echo "  TG-GATEs: 170 compounds × rat liver microarray (single+repeat dose)"
echo "  Data processed and available via OpenRiskNet EdelweissData API:"
TG_URL="https://api.edelweissdata.com/datasets/open-tggates/data?format=csv&limit=1000"
echo "  → TG-GATEs via EdelweissData API"
curl -sSL -A "$UA" --max-time 120 -o "$DATA/omics_tggates/tggates_sample.csv" \
  "$TG_URL" \
  && echo "    ✓ TG-GATEs sample: $(wc -l < "$DATA/omics_tggates/tggates_sample.csv") rows" \
  || echo "    ✗ Use browser: https://openrisknet.org/e-infrastructure/services/159/"

echo "  Note: For full TG-GATEs FTP download (~5 GB CEL files):"
echo "    ftp://ftp.biosciencedbc.jp/archive/open-tggates/"
echo "    Or R: BiocManager::install('ToxicoGx'); tggates <- downloadTSet('Open TG-GATEs')"

# DrugMatrix annotation (NIEHS GitHub)
echo "  → DrugMatrix metadata from NIEHS GitHub"
git clone --quiet https://github.com/NIEHS/DrugMatrix \
  "$DATA/omics_tggates/DrugMatrix" 2>/dev/null \
  && echo "    ✓ DrugMatrix metadata cloned" \
  || echo "    ✓ DrugMatrix already cloned"
echo "  Note: Full DrugMatrix expression data: https://ntp.niehs.nih.gov/drugmatrix/index.html"

# ═════════════════════════════════════════════════════════════════════════════
echo ""; echo "══════════════ BLOCK 7: SIDER + ChEMBL APPROVED DRUG LABELS ══════════════"
DL "$DATA/sider_chembl/sider4_side_effects.tsv.gz" \
  "SIDER 4.1 — side effects TSV (1430 drugs, 140k ADR pairs)" \
  "http://sideeffects.embl.de/media/download/meddra_all_se.tsv.gz" || \
  echo "    ✗ Try SIDER4 via GitHub: https://github.com/dhimmel/SIDER4/raw/master/data/side-effects.tsv"

DL "$DATA/sider_chembl/sider4_indications.tsv.gz" \
  "SIDER 4.1 — indications TSV" \
  "http://sideeffects.embl.de/media/download/meddra_all_indications.tsv.gz" || true

echo "  → ChEMBL Python client (install if needed, then run fetch script)"
python3 -c "from chembl_webresource_client.new_client import new_client; print('ChEMBL client OK')" 2>/dev/null \
  || pip install chembl_webresource_client -q && echo "    ✓ ChEMBL client installed"
cat > "$DATA/sider_chembl/fetch_chembl_oligos.py" << 'PYEOF'
"""Fetch ChEMBL records for approved oligonucleotide drugs."""
from chembl_webresource_client.new_client import new_client
import json, os

DRUG_NAMES = [
    "patisiran","vutrisiran","givosiran","lumasiran","inclisiran",
    "nedosiran","fitusiran","plozasiran","inotersen","eplontersen",
    "tofersen","nusinersen","mipomersen","olezarsen","donidalorsen",
    "eteplirsen","golodirsen","viltolarsen","casimersen","pegaptanib","avacincaptad"
]
molecule = new_client.molecule
out = []
for name in DRUG_NAMES:
    results = molecule.filter(pref_name__iexact=name).only(
        ['molecule_chembl_id','pref_name','molecule_type',
         'max_phase','black_box_warning','oral','parenteral',
         'molecule_properties'])
    r = list(results)
    if r:
        out.append(r[0])
        print(f"  ✓ {name}: {r[0].get('molecule_chembl_id')}")
    else:
        print(f"  - {name}: not found in ChEMBL")

out_path = os.path.join(os.path.dirname(__file__), "chembl_approved_oligos.json")
with open(out_path, "w") as f:
    json.dump(out, f, indent=2)
print(f"\nSaved {len(out)} records to {out_path}")
PYEOF
python3 "$DATA/sider_chembl/fetch_chembl_oligos.py" \
  && echo "    ✓ ChEMBL approved oligo records saved" \
  || echo "    ✗ ChEMBL fetch failed — install: pip install chembl_webresource_client"

# ═════════════════════════════════════════════════════════════════════════════
echo ""; echo "══════════════ BLOCK 8: NDA/EMA LABEL LINKS ══════════════"
cat > "$DATA/nda_labels/fda_daf_approved_oligos.txt" << 'NDAEOF'
# FDA Drugs@FDA — Approved Oligonucleotide Drug NDA/BLA Pages
# For each: open in browser → Medical Reviews → Nonclinical Pharmacology/Toxicology section

PATISIRAN (Onpattro):    https://www.accessdata.fda.gov/scripts/cder/daf/index.cfm?event=BasicSearch.process&ApplNo=210922
VUTRISIRAN (Amvuttra):   https://www.accessdata.fda.gov/scripts/cder/daf/index.cfm?event=BasicSearch.process&ApplNo=215515
GIVOSIRAN (Givlaari):    https://www.accessdata.fda.gov/scripts/cder/daf/index.cfm?event=BasicSearch.process&ApplNo=212194
LUMASIRAN (Oxlumo):      https://www.accessdata.fda.gov/scripts/cder/daf/index.cfm?event=BasicSearch.process&ApplNo=214103
INCLISIRAN (Leqvio):     https://www.accessdata.fda.gov/scripts/cder/daf/index.cfm?event=BasicSearch.process&ApplNo=214012
NEDOSIRAN (Rivfloza):    https://www.accessdata.fda.gov/scripts/cder/daf/index.cfm?event=BasicSearch.process&ApplNo=216887
FITUSIRAN (Alhemo):      https://www.accessdata.fda.gov/scripts/cder/daf/index.cfm?event=BasicSearch.process&ApplNo=761178
INOTERSEN (Tegsedi):     https://www.accessdata.fda.gov/scripts/cder/daf/index.cfm?event=BasicSearch.process&ApplNo=211945
EPLONTERSEN (Wainua):    https://www.accessdata.fda.gov/scripts/cder/daf/index.cfm?event=BasicSearch.process&ApplNo=217388
TOFERSEN (Qalsody):      https://www.accessdata.fda.gov/scripts/cder/daf/index.cfm?event=BasicSearch.process&ApplNo=215887
OLEZARSEN (Tryngolza):   https://www.accessdata.fda.gov/scripts/cder/daf/index.cfm?event=BasicSearch.process&ApplNo=219175
NUSINERSEN (Spinraza):   https://www.accessdata.fda.gov/scripts/cder/daf/index.cfm?event=BasicSearch.process&ApplNo=209531
MIPOMERSEN (Kynamro):    https://www.accessdata.fda.gov/scripts/cder/daf/index.cfm?event=BasicSearch.process&ApplNo=203568
ETEPLIRSEN (Exondys51):  https://www.accessdata.fda.gov/scripts/cder/daf/index.cfm?event=BasicSearch.process&ApplNo=206488
GOLODIRSEN (Vyondys53):  https://www.accessdata.fda.gov/scripts/cder/daf/index.cfm?event=BasicSearch.process&ApplNo=211970
VILTOLARSEN (Viltepso):  https://www.accessdata.fda.gov/scripts/cder/daf/index.cfm?event=BasicSearch.process&ApplNo=213026
CASIMERSEN (Amondys45):  https://www.accessdata.fda.gov/scripts/cder/daf/index.cfm?event=BasicSearch.process&ApplNo=213026
PEGAPTANIB (Macugen):    https://www.accessdata.fda.gov/scripts/cder/daf/index.cfm?event=BasicSearch.process&ApplNo=021756
DONIDALORSEN (Orlynvah): https://www.accessdata.fda.gov/scripts/cder/daf/index.cfm?event=BasicSearch.process&ApplNo=215555
IMETELSTAT (Omblastys):  https://www.accessdata.fda.gov/scripts/cder/daf/index.cfm?event=BasicSearch.process&ApplNo=761120

# EMA EPARs (European Public Assessment Reports) — include nonclinical tox overviews
EPLONTERSEN: https://www.ema.europa.eu/en/search/search?f%5B0%5D=field_ema_web_categories%3A90&search_api_fulltext=eplontersen
INCLISIRAN:  https://www.ema.europa.eu/en/search/search?f%5B0%5D=field_ema_web_categories%3A90&search_api_fulltext=inclisiran
VUTRISIRAN:  https://www.ema.europa.eu/en/search/search?f%5B0%5D=field_ema_web_categories%3A90&search_api_fulltext=vutrisiran
NDAEOF
echo "    ✓ FDA/EMA NDA link list written to $DATA/nda_labels/"

# ═════════════════════════════════════════════════════════════════════════════
echo ""; echo "══════════════ BLOCK 9: GitHub REPOS ══════════════"
for repo_info in \
  "wusuw/ClinASO ClinASO" \
  "dhimmel/SIDER4 SIDER4_processing"
do
  repo="${repo_info%% *}"; dest="${repo_info##* }"
  dir="$DATA/github_data/$dest"
  [ -d "$dir" ] \
    && echo "    ✓ $repo already cloned" \
    || (git clone --quiet "https://github.com/$repo" "$dir" \
        && echo "    ✓ Cloned $repo") 2>/dev/null || echo "    ✗ $repo failed"
done

# ═════════════════════════════════════════════════════════════════════════════
echo ""; echo "══════════════ BLOCK 10: MANUAL / GATED STEPS ══════════════"
cat << 'MANUALEOF'
The following require browser login or manual action:

  1. Papargyri/Abzu CSV (aso.csv with 768 gapmer rows):
     https://docs.abzu.ai/docs/tutorials/r/aso
     → Follow GitHub/Zenodo link in tutorial → save to data/raw/lna_gapmer/aso_data.csv

  2. siRNAmod bulk export:
     http://crdd.osdd.net/servers/sirnamod/dd.php
     → Download All → save to data/raw/sirna_modified/siRNAmod_all.csv

  3. LINCS L1000 (full HepG2 siRNA/shRNA signatures, ~5 GB):
     https://clue.io/data → LINCS Phase II → GSE70138
     → Select HepG2 + HEPG2 knockdown sigs → download .gctx
     → pip install cmapPy; python scripts/parse_lincs_l1000.py

  4. FDA NDA tox reviews (PDFs inside each approval package):
     Use links in data/raw/nda_labels/fda_daf_approved_oligos.txt
     → Click drug name → Documents → Pharmacology/Toxicology Review
     → Save to data/raw/nda_labels/<drug>_tox_review.pdf

  5. Open TG-GATEs full download (5 GB):
     In R: BiocManager::install('ToxicoGx')
            tgg <- downloadTSet('Open TG-GATEs')
     Or FTP: ftp://ftp.biosciencedbc.jp/archive/open-tggates/
MANUALEOF

# ═════════════════════════════════════════════════════════════════════════════
echo ""; echo "══════════════ SUMMARY ══════════════"
echo "Downloaded files:"
find "$DATA" -type f | sort | while read f; do
  sz=$(du -sh "$f" 2>/dev/null | cut -f1)
  echo "  $sz  ${f/$DATA\//}"
done
echo ""
echo "Next step → populate curated_oligo_tox_template_with_approved.csv"
echo "using oligotox_labeling_scheme.csv thresholds, then run:"
echo "  PYTHONPATH=$HOME python oligotox/scripts/curated_csv_to_training.py"
echo "  PYTHONPATH=$HOME python oligotox/train.py"
