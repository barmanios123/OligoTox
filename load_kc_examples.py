#!/usr/bin/env python3
"""
load_kc_examples.py
───────────────────
Merges the new 33 Kupffer cell activation positive cases into the existing
training dataset (training_dataset_v2.csv or training_dataset.csv) and
writes training_dataset_v3.csv ready for retraining.

This is a DROP-IN complement to augment_hepatic_labels.py.
Run AFTER augment_hepatic_labels.py (which created training_dataset_v2).

Usage:
    python load_kc_examples.py \
        --kc_csv   data/processed/kupffer_activation_30plus_examples.csv \
        --training data/processed/training_dataset_v2.csv \
        --out      data/processed/training_dataset_v3.csv

Then retrain:
    python train_hepatic_endpoints.py \
        --dataset data/processed/training_dataset_v3.csv \
        --outdir  models/
"""

import argparse
import pandas as pd
from pathlib import Path


def merge_kc_examples(kc_path: Path, training_path: Path, out_path: Path):
    print(f"Loading KC examples from {kc_path}")
    kc = pd.read_csv(kc_path)
    kc["sequence"] = kc["sequence"].str.upper().str.replace(" ", "").str.strip()
    print(f"  {len(kc)} KC rows  |  kupffer_activation=1: {kc.label_kupffer_activation.sum()}")

    print(f"Loading training dataset from {training_path}")
    df = pd.read_csv(training_path)
    df["sequence"] = df["sequence"].str.upper().str.strip()
    for col in ["label_hepatic_vacuolation", "label_kupffer_activation"]:
        if col not in df.columns:
            df[col] = 0
    print(f"  {len(df)} training rows")

    kc_seqs = set(kc["sequence"].values)
    train_seqs = set(df["sequence"].values)

    # Overwrite labels for matching sequences
    overlap = kc_seqs & train_seqs
    print(f"  Sequence overlap (label update): {len(overlap)}")
    for _, row in kc.iterrows():
        mask = df["sequence"] == row["sequence"]
        if mask.any():
            df.loc[mask, "label_kupffer_activation"] = row["label_kupffer_activation"]
            df.loc[mask, "label_hepatic_vacuolation"] = row["label_hepatic_vacuolation"]

    # Append novel KC sequences
    novel = kc[~kc["sequence"].isin(train_seqs)].copy()
    print(f"  Novel sequences appended: {len(novel)}")
    merged = pd.concat([df, novel], ignore_index=True)

    vac_pos = int(merged["label_hepatic_vacuolation"].sum())
    kup_pos = int(merged["label_kupffer_activation"].sum())
    print(f"\nMerged dataset: {len(merged)} total")
    print(f"  hepatic_vacuolation=1 : {vac_pos}")
    print(f"  kupffer_activation=1  : {kup_pos}")

    merged.to_csv(out_path, index=False)
    print(f"\nSaved → {out_path}")
    print("\nNext step: python train_hepatic_endpoints.py --dataset", out_path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--kc_csv",   default="data/processed/kupffer_activation_30plus_examples.csv")
    parser.add_argument("--training", default="data/processed/training_dataset_v2.csv")
    parser.add_argument("--out",      default="data/processed/training_dataset_v3.csv")
    args = parser.parse_args()
    merge_kc_examples(Path(args.kc_csv), Path(args.training), Path(args.out))


if __name__ == "__main__":
    main()
