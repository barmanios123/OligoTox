#!/usr/bin/env python3
"""
oligotox_geo_fetch.py
=====================
Fetches GEO series metadata + supplemental files for OligoTox training.
Requires: pip install GEOparse requests

Usage:
  cd ~/oligotox
  python oligotox/scripts/geo_fetch.py

Writes to: data/raw/tier3_mechanistic_omics/geo/
"""

import os
import sys
import json
import urllib.request
import urllib.error

GEO_ACCESSIONS = [
    # (accession, description, primary_endpoint, priority)
    ("GSE47334", "siRNA_offtarget_seed_match_HeLa_transcriptomics", "offtarget_gene_silencing", 1),
    ("GSE68580", "siRNA_RISC_innate_immune_HEK293", "innate_immune_activation", 1),
    ("GSE27638", "ASO_EGFR_hepatic_fibrosis_liver_transcriptomics", "hepatocellular_injury", 1),
    ("GSE152641", "innate_immune_two_gene_marker_COVID", "innate_immune_activation", 2),
    ("GSE171110", "LNP_lipid_nanoparticle_liver_immune_mouse", "complement_activation", 2),
]

OUT_DIR = "data/raw/tier3_mechanistic_omics/geo"
os.makedirs(OUT_DIR, exist_ok=True)

def fetch_geo_metadata(acc):
    """Fetch GEO series metadata via NCBI e-utilities."""
    url = (
        f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        f"?db=gds&term={acc}[Accession]&retmax=5&retmode=json"
    )
    try:
        with urllib.request.urlopen(url, timeout=30) as resp:
            data = json.loads(resp.read().decode())
        ids = data.get("esearchresult", {}).get("idlist", [])
        print(f"  {acc}: found {len(ids)} GDS record(s)")
        return ids
    except Exception as e:
        print(f"  {acc}: metadata fetch failed — {e}")
        return []

def fetch_geo_soft(acc):
    """Download GEO SOFT format (metadata + sample info)."""
    soft_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{acc[:-3]}nnn/{acc}/soft/{acc}_family.soft.gz"
    out_path = os.path.join(OUT_DIR, f"{acc}_family.soft.gz")
    if os.path.exists(out_path):
        print(f"  {acc}: SOFT file already exists, skipping")
        return True
    try:
        urllib.request.urlretrieve(soft_url, out_path)
        size_kb = os.path.getsize(out_path) // 1024
        print(f"  {acc}: SOFT downloaded ({size_kb} KB)")
        return True
    except Exception as e:
        print(f"  {acc}: SOFT download failed — {e}")
        # Try mini-SOFT
        mini_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{acc[:-3]}nnn/{acc}/miniml/{acc}_family.xml.tgz"
        try:
            urllib.request.urlretrieve(mini_url, out_path.replace(".soft.gz", "_miniml.xml.tgz"))
            print(f"  {acc}: MINiML fallback downloaded")
            return True
        except Exception as e2:
            print(f"  {acc}: MINiML also failed — {e2}")
            return False

def main():
    print("OligoTox GEO Fetcher")
    print(f"Output directory: {OUT_DIR}")
    print()

    results = []
    for acc, desc, endpoint, priority in sorted(GEO_ACCESSIONS, key=lambda x: x[3]):
        print(f"[P{priority}] {acc} — {desc}")
        ids = fetch_geo_metadata(acc)
        ok = fetch_geo_soft(acc)
        results.append({
            "accession": acc,
            "description": desc,
            "primary_endpoint": endpoint,
            "priority": priority,
            "gds_ids": ids,
            "downloaded": ok,
        })

    # Write manifest
    manifest_path = os.path.join(OUT_DIR, "geo_manifest.json")
    with open(manifest_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nManifest written: {manifest_path}")

    downloaded = sum(1 for r in results if r["downloaded"])
    print(f"Downloaded: {downloaded}/{len(results)}")
    print("\nNext step: use GEOparse to extract expression matrices from SOFT files.")
    print("  pip install GEOparse")
    print("  import GEOparse; gse = GEOparse.get_GEO(filepath='data/raw/.../GSExxx_family.soft.gz')")

if __name__ == "__main__":
    main()
