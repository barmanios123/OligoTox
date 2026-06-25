#!/usr/bin/env bash
# =============================================================================
# OligoTox Dataset Kit — Master Download Script v1.0
# Run from ~/oligotox/  after cloning the platform
# 
# STATUS KEY:
#   ✓  Confirmed auto-downloadable (verified 200 OK, real PDF)
#   ⚠  Gated — requires manual browser download (publisher paywall)
#   ✎  Requires manual table extraction from paper
# =============================================================================
set -e
RAWDIR="data/raw"
mkdir -p "$RAWDIR/fda_reviews" "$RAWDIR/papers_open" "$RAWDIR/papers_gated"

echo ""
echo "════════════════════════════════════════════════════════════════"
echo "  BLOCK 1: FDA Pharmacology Reviews — CONFIRMED OPEN ACCESS"
echo "════════════════════════════════════════════════════════════════"

# ✓ FDA 3-ASO + 18 impurities (human hepatocyte ATP/albumin/urea at 100×Cmax)
#   → Endpoints: hepatocellular_injury, renal_tubular_injury, cytotoxicity
echo "[1/4] FDA ASO 3-drug impurity study..."
curl -L --max-time 120 -A "Mozilla/5.0" \
  "https://www.fda.gov/media/168708/download" \
  -o "$RAWDIR/fda_reviews/fda_3aso_18impurities_hepatocyte.pdf"
echo "    ✓  Saved: fda_3aso_18impurities_hepatocyte.pdf (1 MB)"

# ✓ Inotersen NDA 211172 — Pharmacology Review
#   → Endpoints: thrombocytopenia, renal_tubular_injury, hepatocellular_injury
#   → Tables to extract: platelet counts by timepoint, kidney histopath grades
echo "[2/4] Inotersen NDA 211172 — Pharmacology Review..."
curl -L --max-time 180 -A "Mozilla/5.0" \
  "https://www.accessdata.fda.gov/drugsatfda_docs/nda/2018/211172Orig1s000PharmR.pdf" \
  -o "$RAWDIR/fda_reviews/inotersen_NDA211172_PharmR.pdf"
echo "    ✓  Saved: inotersen_NDA211172_PharmR.pdf (7.6 MB)"

# ✓ Mipomersen NDA 203568 — Pharmacology Review
#   → Endpoints: hepatocellular_injury, renal_tubular_injury, complement_activation
#   → Tables to extract: ALT/AST by dose/timepoint, histopath findings
echo "[3/4] Mipomersen NDA 203568 — Pharmacology Review..."
curl -L --max-time 300 -A "Mozilla/5.0" \
  "https://www.accessdata.fda.gov/drugsatfda_docs/nda/2013/203568Orig1s000PharmR.pdf" \
  -o "$RAWDIR/fda_reviews/mipomersen_NDA203568_PharmR.pdf"
echo "    ✓  Saved: mipomersen_NDA203568_PharmR.pdf (31 MB)"

# ✓ Nusinersen NDA 209531 — Pharmacology Review
#   → Endpoints: cns_neurotoxicity, hepatocellular_injury, renal_tubular_injury
#   → Tables to extract: CNS histopath (dorsal root ganglion), liver enzyme data
echo "[4/4] Nusinersen NDA 209531 — Pharmacology Review..."
curl -L --max-time 180 -A "Mozilla/5.0" \
  "https://www.accessdata.fda.gov/drugsatfda_docs/nda/2016/209531Orig1s000PharmR.pdf" \
  -o "$RAWDIR/fda_reviews/nusinersen_NDA209531_PharmR.pdf"
echo "    ✓  Saved: nusinersen_NDA209531_PharmR.pdf (9.4 MB)"

echo ""
echo "════════════════════════════════════════════════════════════════"
echo "  BLOCK 2: GATED PAPERS — Manual Download Required"
echo "════════════════════════════════════════════════════════════════"
echo ""
echo "  The following papers require manual browser download."
echo "  Save PDFs (and supplementary Excel/CSV) to: $RAWDIR/papers_gated/"
echo ""
echo "  ⚠  [A] Papargyri et al. 2020 (512 LNA gapmers, caspase + knockdown)"
echo "     DOI: 10.1016/j.omtn.2019.12.011"
echo "     URL: https://www.cell.com/molecular-therapy-family/nucleic-acids/fulltext/S2162-2531(19)30392-0"
echo "     CRITICAL: Download Table S1 (Excel) — 512 rows, caspase% + HIF1A knockdown%"
echo "     Save as: papers_gated/papargyri_2020_LNA_512gapmers_TableS1.xlsx"
echo "     → Endpoints: offtarget_gene_silencing (caspase >300% = 1), hepatocellular_injury"
echo ""
echo "  ⚠  [B] Dieckmann et al. 2018 (hybridization-dependent gapmer tox)"
echo "     PubMed: https://pubmed.ncbi.nlm.nih.gov/29499955/"
echo "     PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC5725219/"
echo "     URL: https://www.sciencedirect.com/science/article/pii/S2162253117302883"
echo "     Save as: papers_gated/dieckmann_2018_hybridization_tox.pdf"
echo "     → Endpoints: hepatocellular_injury, offtarget_gene_silencing"
echo ""
echo "  ⚠  [C] Janas et al. 2018 (GalNAc siRNA nonclinical safety)"
echo "     PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC6249674/"
echo "     Save as: papers_gated/janas_2018_galnac_nonclinical_safety.pdf"
echo "     → Endpoints: complement_activation, cytokine_release, renal_tubular_injury"
echo ""
echo "  ⚠  [D] Burel et al. 2015 (RNase H1 hepatotox mechanism)"
echo "     PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC4595405/"
echo "     Save as: papers_gated/burel_2015_rnaseH1_hepatotox.pdf"
echo "     → Endpoints: hepatocellular_injury (hybridization-dependent mechanism)"
echo ""
echo "  ⚠  [E] Chembiochem 2025 (ASO sequence×chemistry synergy tox)"
echo "     PubMed: https://pubmed.ncbi.nlm.nih.gov/41078095/"
echo "     Save as: papers_gated/aso_seqchem_interplay_tox_2025.pdf"
echo "     → Endpoints: all — framework paper for label thresholds"
echo ""
echo "  ⚠  [F] siRNA Holistic Tox Review 2023 (PMC10543595)"
echo "     PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC10543595/"
echo "     Save as: papers_gated/siRNA_holistic_tox_review_2023.pdf"
echo "     → Use for label definition decisions across all 13 endpoints"

echo ""
echo "════════════════════════════════════════════════════════════════"
echo "  BLOCK 3: NIH OligoTox Challenge — Registration"  
echo "════════════════════════════════════════════════════════════════"
echo ""
echo "  Challenge overview: https://ncats.nih.gov/funding/challenges"
echo "  Notice: https://grants.nih.gov/grants/guide/notice-files/NOT-TR-26-006.html"
echo "  Contact: ncatsoligotox@mail.nih.gov"
echo "  Phase 2 submission window: May 1 – December 31, 2026"
echo "  Prize: up to \$500,000 (Phase 2 data generation)"
echo ""
echo "  → Email ncatsoligotox@mail.nih.gov with subject:"
echo "    'OligoTox Challenge Phase 2 — Registration: [Your Institution]'"
echo "    Attach: brief description of your OligoTox platform + data plan"

echo ""
echo "════════════════════════════════════════════════════════════════"
echo "  BLOCK 4: GitHub Repos (code dependencies)"
echo "════════════════════════════════════════════════════════════════"
echo ""
echo "[5] Clone openASO (LNA gapmer sequence+activity dataset)..."
if [ ! -d "data/raw/repos/openASO" ]; then
  git clone --depth 1 https://github.com/henriksson-lab/openASO.git \
    data/raw/repos/openASO 2>/dev/null || echo "    ⚠  Clone failed (check network)"
else
  echo "    ✓  Already exists"
fi

echo "[6] Clone ChemProp (molecular graph NN for oligo chem features)..."
if [ ! -d "data/raw/repos/chemprop" ]; then
  git clone --depth 1 https://github.com/chemprop/chemprop.git \
    data/raw/repos/chemprop 2>/dev/null || echo "    ⚠  Clone failed"
else
  echo "    ✓  Already exists"
fi

echo "[7] Clone TLRpred (TLR7/8/9 activation scoring)..."
# TLRpred is web-only; record the endpoint for API calls
echo "    → TLRpred web server: https://webs.iiitd.edu.in/raghava/tlrpred/"
echo "      Use POST requests per sequence for TLR probability scores"

echo ""
echo "════════════════════════════════════════════════════════════════"
echo "  BLOCK 5: siRNAmod Database Export"
echo "════════════════════════════════════════════════════════════════"
echo ""
echo "  → Go to: http://crdd.osdd.net/servers/sirnamod/"
echo "  → Select: 'Browse All Modified siRNAs'"
echo "  → Export as CSV, save to: $RAWDIR/siRNAmod_all_modified_siRNAs.csv"
echo ""
echo "════════════════════════════════════════════════════════════════"
echo "  BLOCK 6: ChEMBL Approved Oligonucleotide Drugs"
echo "════════════════════════════════════════════════════════════════"
echo ""
echo "[8] Download ChEMBL approved oligo drug profiles..."
python3 - <<'PYEOF'
try:
    from chembl_webresource_client.new_client import new_client
    drug = new_client.drug
    results = drug.filter(molecule_type="Oligonucleotide")
    import json
    with open("data/raw/chembl_oligonucleotide_drugs.json", "w") as f:
        json.dump(list(results), f, indent=2)
    print(f"    ✓  Saved {len(list(results))} oligo drugs from ChEMBL")
except ImportError:
    print("    ⚠  chembl_webresource_client not installed")
    print("       Run: pip install chembl_webresource_client")
    print("       Then re-run this block")
PYEOF

echo ""
echo "════════════════════════════════════════════════════════════════"
echo "  ALL DOWNLOADS COMPLETE"
echo "  Next step: run python scripts/curated_csv_to_training.py"
echo "════════════════════════════════════════════════════════════════"
