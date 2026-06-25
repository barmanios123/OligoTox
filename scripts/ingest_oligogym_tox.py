#!/usr/bin/env python3
"""
OligoTox ingest_oligogym_tox.py

Ingests toxicity-relevant datasets from OligoGym (already in data/oligogym/)
into a standardised format for integration with the OligoTox pipeline.

OligoGym schema: x, y, y_raw, targets, smiles, fasta [, cell_line, concentration]
  - y      = binarized label (already thresholded by OligoGym)
  - y_raw  = raw continuous readout (assay-specific units)
  - fasta  = nucleotide sequence

Output: data/oligogym/oligogym_tox_combined.csv

Run from repo root:
  PYTHONPATH="$HOME:$PWD:$PWD/src" python3 scripts/ingest_oligogym_tox.py
"""
import pandas as pd
import numpy as np
from pathlib import Path

OLIGOGYM_DIR = Path("data/oligogym")
OUT_PATH = OLIGOGYM_DIR / "oligogym_tox_combined.csv"

# Dataset definitions: (filename, endpoint_col_name, label_desc, binarize_threshold, higher_is_toxic)
# For MOE_Neurotox: y is 0 (safe) > 7 (toxic) — OligoGym pre-binarised; y=1 means toxic
# For Hagedorn_2022 (Neurotox LNA): y=1 means neurotoxic (calcium oscillation)
# For Papargyri_2020 (Cytotox LNA): y=1 means cytotoxic (caspase assay)
# For Ichihara_2007: y=1 means inhibitory activity
# For McQuisten_2007: y=1 means active
DATASET_CONFIG = [
    {
        "file": "MOE_Neurotox_1_processed.csv.gz",
        "endpoint": "endpoint_cns_neurotoxicity_invivo",
        "source": "MOE_Neurotox_2024",
        "modality": "ASO",
        "assay_type": "in_vivo",
        "readout": "FOB_score_3h",
        "label_polarity": "higher_toxic",  # y=1 = toxic (FOB > threshold)
    },
    {
        "file": "Hagedorn_2022_1_processed.csv.gz",
        "endpoint": "endpoint_cns_neurotoxicity_invitro",
        "source": "Hagedorn_2022",
        "modality": "ASO",
        "assay_type": "in_vitro",
        "readout": "calcium_oscillation_score",
        "label_polarity": "higher_toxic",
    },
    {
        "file": "Papargyri_2020_1_processed.csv.gz",
        "endpoint": "endpoint_cytotoxicity",
        "source": "Papargyri_2020_OligoGym",
        "modality": "ASO",
        "assay_type": "in_vitro",
        "readout": "caspase_minmax_scaled",
        "label_polarity": "higher_toxic",
    },
    {
        "file": "Ichihara_2007_1_processed.csv.gz",
        "endpoint": "endpoint_cell_inhibition",
        "source": "Ichihara_2007_1",
        "modality": "siRNA",
        "assay_type": "in_vitro",
        "readout": "pct_inhibition",
        "label_polarity": "higher_active",  # not toxicity per se
    },
    {
        "file": "Ichihara_2007_2_processed.csv.gz",
        "endpoint": "endpoint_cell_inhibition_2",
        "source": "Ichihara_2007_2",
        "modality": "siRNA",
        "assay_type": "in_vitro",
        "readout": "pct_inhibition",
        "label_polarity": "higher_active",
    },
    {
        "file": "McQuisten_2007_1_processed.csv.gz",
        "endpoint": "endpoint_target_inhibition",
        "source": "McQuisten_2007",
        "modality": "ASO",
        "assay_type": "in_vitro",
        "readout": "activity_range_0to1",
        "label_polarity": "higher_active",
    },
]

def extract_sequence(row):
    """Extract nucleotide sequence from fasta field."""
    fasta = str(row.get("fasta",""))
    if fasta and fasta != "nan":
        # fasta format: >header\nSEQUENCE or just sequence
        lines = fasta.strip().split("\n")
        seq_lines = [l for l in lines if not l.startswith(">")]
        return "".join(seq_lines).upper().replace(" ","")
    return ""

def gc_content(seq):
    n = len(seq)
    if n == 0: return np.nan
    return (seq.count("G") + seq.count("C")) / n

def length_nt(seq):
    return len(seq) if seq else np.nan

def main():
    all_dfs = []
    
    for cfg in DATASET_CONFIG:
        fpath = OLIGOGYM_DIR / cfg["file"]
        if not fpath.exists():
            print(f"MISSING: {fpath} — skipping")
            continue
        
        df = pd.read_csv(fpath)
        print(f"\nLoaded: {cfg['file']}")
        print(f"  Shape: {df.shape}")
        print(f"  Columns: {df.columns.tolist()}")
        print(f"  y distribution: {df['y'].value_counts().to_dict()}")
        
        # Build standardised record
        records = []
        for _, row in df.iterrows():
            seq = extract_sequence(row)
            gc = gc_content(seq)
            length = length_nt(seq)
            
            rec = {
                "source": cfg["source"],
                "modality": cfg["modality"],
                "assay_type": cfg["assay_type"],
                "readout": cfg["readout"],
                "sequence": seq,
                "length_nt": length,
                "gc_content": gc,
                "y_raw": row.get("y_raw", np.nan),
                "y_binary": int(row.get("y", 0)),
                cfg["endpoint"]: int(row.get("y", 0)),
            }
            records.append(rec)
        
        ds_df = pd.DataFrame(records)
        ds_df["dataset"] = cfg["file"].replace("_processed.csv.gz","")
        all_dfs.append(ds_df)
        print(f"  Processed: {len(ds_df)} records, {int(ds_df['y_binary'].sum())} positives")
    
    if not all_dfs:
        print("No datasets found! Check data/oligogym/ directory.")
        return
    
    combined = pd.concat(all_dfs, axis=0, ignore_index=True, sort=False)
    combined.to_csv(OUT_PATH, index=False)
    print(f"\n=== COMBINED OUTPUT ===")
    print(f"Saved: {OUT_PATH}")
    print(f"Shape: {combined.shape}")
    print(f"Sources: {combined['source'].value_counts().to_dict()}")
    print(f"\nEndpoint label summary:")
    for ep_col in [c for c in combined.columns if c.startswith("endpoint_")]:
        n_pos = int(combined[ep_col].fillna(0).sum())
        n_tot = combined[ep_col].notna().sum()
        print(f"  {ep_col}: {n_pos}/{n_tot} positives")

if __name__ == "__main__":
    main()
