#!/usr/bin/env python3
"""
oligotox_sider_parse.py
=======================
Parses SIDER 4.1 TSV files to extract adverse event ground-truth labels
for the 13 approved oligonucleotide drugs in your benchmark module.

Input:  data/raw/tier1_labelable/sider4_meddra_all_se.tsv.gz
        data/raw/tier1_labelable/sider4_drug_names.tsv
Output: data/raw/tier1_labelable/sider_oligo_drugs_ae_labels.csv

Usage:
  cd ~/oligotox
  python oligotox/scripts/sider_parse.py
"""

import os
import gzip
import csv

# Approved oligonucleotide drugs to extract from SIDER
# Keys are SIDER stitch IDs (CID-based) or drug name substrings
OLIGO_DRUGS = {
    "mipomersen":    {"modality": "ASO_gapmer", "chemistry": "2-MOE/PS"},
    "inotersen":     {"modality": "ASO_gapmer", "chemistry": "2-MOE/PS"},
    "nusinersen":    {"modality": "SSO",         "chemistry": "2-MOE/PS"},
    "inclisiran":    {"modality": "GalNAc-siRNA","chemistry": "2-OMe/PS"},
    "givosiran":     {"modality": "GalNAc-siRNA","chemistry": "2-OMe/PS"},
    "lumasiran":     {"modality": "GalNAc-siRNA","chemistry": "2-OMe/PS"},
    "vutrisiran":    {"modality": "GalNAc-siRNA","chemistry": "2-OMe/PS"},
    "patisiran":     {"modality": "LNP-siRNA",   "chemistry": "2-OMe/PS"},
    "eplontersen":   {"modality": "ASO_gapmer",  "chemistry": "2-MOE/PS/GalNAc"},
    "nedosiran":     {"modality": "GalNAc-siRNA","chemistry": "2-OMe/PS"},
}

# MedDRA terms → OligoTox 13-endpoint mapping
ENDPOINT_MAP = {
    "hepatocellular_injury":     ["hepatic", "liver", "hepatotoxic", "hepatocellular",
                                   "alanine aminotransferase", "aspartate aminotransferase",
                                   "transaminase", "hepatitis", "jaundice"],
    "renal_tubular_injury":      ["renal", "kidney", "nephrotoxic", "tubular", "proteinuria",
                                   "creatinine", "glomerular", "nephritis"],
    "thrombocytopenia":          ["thrombocytopenia", "platelet", "platelet count decreased"],
    "complement_activation":     ["complement", "urticaria", "angioedema"],
    "innate_immune_activation":  ["injection site reaction", "flu-like", "influenza",
                                   "pyrexia", "chills", "cytokine", "interferon"],
    "cytokine_release":          ["cytokine release", "infusion related", "fever",
                                   "hypotension", "flushing"],
    "cns_delayed_neurotoxicity": ["neurotoxic", "peripheral neuropathy", "neuropathy",
                                   "tremor", "dizziness", "headache", "seizure"],
    "offtarget_gene_silencing":  [],   # not directly in SIDER — use transcriptomics
    "offtarget_mirna_loading":   [],   # not in SIDER
}

def parse_sider(drug_names_path, se_path, out_path):
    """Extract oligo drug adverse events from SIDER."""
    print("Parsing SIDER 4.1...")

    # Load drug name → stitch ID mapping
    drug_to_stitch = {}
    if not os.path.exists(drug_names_path):
        print(f"  NOT FOUND: {drug_names_path}")
        print("  Run oligotox_dataset_download.sh first to download SIDER files.")
        return

    with open(drug_names_path, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                stitch_id, name = parts[0], parts[1].lower()
                for drug in OLIGO_DRUGS:
                    if drug in name:
                        drug_to_stitch[drug] = stitch_id

    print(f"  Matched {len(drug_to_stitch)} oligo drugs in SIDER: {list(drug_to_stitch.keys())}")

    if not drug_to_stitch:
        print("  WARNING: No oligo drugs matched in SIDER drug_names.tsv.")
        print("  Note: Newer oligos (inclisiran, vutrisiran) may not be in SIDER 4.1 (2015).")
        print("  Use UniTox (unitox_2418_drugs.csv) for newer approvals.")
        return

    # Parse adverse events
    if not os.path.exists(se_path):
        print(f"  NOT FOUND: {se_path}")
        return

    results = []
    stitch_ids_of_interest = set(drug_to_stitch.values())

    open_fn = gzip.open if se_path.endswith(".gz") else open
    with open_fn(se_path, "rt", encoding="utf-8") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            stitch_flat, stitch_stereo, umls_cid, meddra_term = parts[:4]
            if stitch_flat in stitch_ids_of_interest or stitch_stereo in stitch_ids_of_interest:
                # Map to OligoTox endpoint
                term_lower = meddra_term.lower()
                matched_endpoints = []
                for ep, keywords in ENDPOINT_MAP.items():
                    if any(kw in term_lower for kw in keywords):
                        matched_endpoints.append(ep)
                results.append({
                    "stitch_flat": stitch_flat,
                    "meddra_term": meddra_term,
                    "endpoints": "|".join(matched_endpoints) if matched_endpoints else "other",
                })

    print(f"  Extracted {len(results)} adverse event records for oligo drugs")

    with open(out_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["stitch_flat", "meddra_term", "endpoints"])
        writer.writeheader()
        writer.writerows(results)
    print(f"  Output: {out_path}")

if __name__ == "__main__":
    parse_sider(
        drug_names_path="data/raw/tier1_labelable/sider4_drug_names.tsv",
        se_path="data/raw/tier1_labelable/sider4_meddra_all_se.tsv.gz",
        out_path="data/raw/tier1_labelable/sider_oligo_drugs_ae_labels.csv",
    )
