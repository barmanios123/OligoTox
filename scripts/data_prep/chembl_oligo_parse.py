#!/usr/bin/env python3
"""
oligotox_chembl_oligo_parse.py
==============================
Fetches approved oligonucleotide drug profiles from ChEMBL API
and builds a structured CSV with safety-relevant fields.

Output: data/raw/regulatory/chembl_oligo_drug_profiles.csv

Usage:
  cd ~/oligotox
  pip install chembl_webresource_client pandas
  python oligotox/scripts/chembl_oligo_parse.py
"""

import json
import os
import urllib.request

OUT_DIR = "data/raw/regulatory"
os.makedirs(OUT_DIR, exist_ok=True)

# ChEMBL IDs for approved oligonucleotide drugs (Phase 4)
OLIGO_CHEMBL_IDS = {
    "CHEMBL1201428": "Mipomersen",
    "CHEMBL2105710": "Nusinersen",
    "CHEMBL3545100": "Inotersen",
    "CHEMBL4298808": "Inclisiran",
    "CHEMBL4523659": "Givosiran",
    "CHEMBL4526018": "Lumasiran",
    "CHEMBL4593498": "Vutrisiran",
    "CHEMBL1201902": "Patisiran",
    "CHEMBL4802765": "Eplontersen",
    "CHEMBL4802099": "Nedosiran",
}

BASE_URL = "https://www.ebi.ac.uk/chembl/api/data"

def fetch_chembl(endpoint_path, params=""):
    url = f"{BASE_URL}/{endpoint_path}.json?{params}"
    try:
        with urllib.request.urlopen(url, timeout=30) as resp:
            return json.loads(resp.read().decode())
    except Exception as e:
        print(f"  API error for {url}: {e}")
        return None

def get_molecule(chembl_id):
    return fetch_chembl(f"molecule/{chembl_id}")

def get_mechanism(chembl_id):
    return fetch_chembl("mechanism", f"molecule_chembl_id={chembl_id}")

def get_atc(chembl_id):
    return fetch_chembl("atc_class", f"molecule_chembl_id={chembl_id}")

def main():
    print("ChEMBL Oligonucleotide Drug Parser")
    print(f"Fetching {len(OLIGO_CHEMBL_IDS)} approved oligo drugs...\n")

    rows = []
    for cid, name in OLIGO_CHEMBL_IDS.items():
        print(f"  {name} ({cid})...")
        mol = get_molecule(cid)
        mech = get_mechanism(cid)

        row = {
            "chembl_id": cid,
            "drug_name": name,
            "max_phase": None,
            "molecule_type": None,
            "oral": None,
            "mw_freebase": None,
            "mechanism_of_action": None,
            "target_chembl_id": None,
            "helm_notation": None,
            "sequence": None,
        }

        if mol:
            props = mol.get("molecule_properties") or {}
            row["max_phase"] = mol.get("max_phase")
            row["molecule_type"] = mol.get("molecule_type")
            row["oral"] = mol.get("oral")
            row["mw_freebase"] = props.get("mw_freebase")
            # HELM notation for oligos
            structures = mol.get("molecule_structures") or {}
            row["helm_notation"] = structures.get("helm")
            row["sequence"] = structures.get("canonical_smiles", "")[:80]

        if mech:
            mechs = mech.get("mechanisms", [])
            if mechs:
                row["mechanism_of_action"] = mechs[0].get("mechanism_of_action")
                row["target_chembl_id"] = mechs[0].get("target_chembl_id")

        rows.append(row)
        print(f"    max_phase={row['max_phase']}, type={row['molecule_type']}, moa={row['mechanism_of_action']}")

    # Write CSV
    import csv
    out_path = os.path.join(OUT_DIR, "chembl_oligo_drug_profiles.csv")
    if rows:
        with open(out_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
        print(f"\nOutput: {out_path} ({len(rows)} drugs)")

    # Also save raw JSON
    json_path = os.path.join(OUT_DIR, "chembl_oligo_drug_profiles_raw.json")
    with open(json_path, "w") as f:
        json.dump(rows, f, indent=2)
    print(f"Raw JSON: {json_path}")

if __name__ == "__main__":
    main()
