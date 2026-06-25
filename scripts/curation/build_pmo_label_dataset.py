#!/usr/bin/env python3
"""
scripts/curation/build_pmo_label_dataset.py

Reads:   data/processed/pmo_label_training_rows.csv   (curator-encoded truth)
Writes:  data/processed/pmo_training_ready.csv         (numeric scores added)
         data/processed/pmo_endpoint_summary.csv        (per-endpoint positive rate)

Run:
    cd ~/oligotox
    PYTHONPATH=$HOME python scripts/curation/build_pmo_label_dataset.py
"""

import csv
import os
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
IN_PATH  = ROOT / "data" / "processed" / "pmo_label_training_rows.csv"
OUT_PATH = ROOT / "data" / "processed" / "pmo_training_ready.csv"
SUM_PATH = ROOT / "data" / "processed" / "pmo_endpoint_summary.csv"

SEVERITY_MAP = {"none": 0, "low": 1, "medium": 2, "high": 3}

ENDPOINT_COLS = [
    "endpoint_hepatocellular_injury",
    "endpoint_renal_tubular_injury",
    "endpoint_renal_accumulation",
    "endpoint_thrombocytopenia",
    "endpoint_coagulation_effect",
    "endpoint_complement_activation",
    "endpoint_innate_immune_activation",
    "endpoint_cytokine_release",
    "endpoint_hypersensitivity",
    "endpoint_cns_neurotoxicity",
    "endpoint_offtarget_gene_silencing",
    "endpoint_reproductive_dev_tox",
    "endpoint_injection_site_reaction",
]

SEVERITY_COLS = [
    "severity_liver",
    "severity_renal",
    "severity_hematology",
    "severity_immune",
]


def load_rows(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def add_numeric_scores(rows):
    out = []
    for r in rows:
        row = dict(r)
        # Severity → numeric score
        for sc in SEVERITY_COLS:
            organ = sc.replace("severity_", "")
            row[f"{organ}_score"] = SEVERITY_MAP.get(row.get(sc, "none"), 0)
        # Composite tox burden (sum of positive endpoints)
        row["tox_endpoint_count"] = sum(
            int(row.get(ep, 0)) for ep in ENDPOINT_COLS
        )
        # Boolean: any serious immune event
        row["any_serious_immune"] = int(
            int(row.get("endpoint_complement_activation", 0)) or
            int(row.get("endpoint_innate_immune_activation", 0)) or
            int(row.get("endpoint_cytokine_release", 0)) or
            int(row.get("endpoint_hypersensitivity", 0))
        )
        # Boolean: renal monitoring required
        row["renal_monitoring_required"] = int(row.get("warning_renal_monitoring", 0))
        out.append(row)
    return out


def write_summary(rows, path):
    """Per-endpoint positive rate across the labelled set."""
    summary = []
    for ep in ENDPOINT_COLS:
        vals = [int(r.get(ep, 0)) for r in rows]
        summary.append({
            "endpoint": ep,
            "n_positive": sum(vals),
            "n_total": len(vals),
            "positive_rate": round(sum(vals) / len(vals), 3) if vals else 0,
        })
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["endpoint", "n_positive", "n_total", "positive_rate"])
        writer.writeheader()
        writer.writerows(summary)
    return summary


def main():
    if not IN_PATH.exists():
        raise FileNotFoundError(f"Input not found: {IN_PATH}")

    rows = load_rows(IN_PATH)
    print(f"Loaded {len(rows)} rows from {IN_PATH.name}")

    enriched = add_numeric_scores(rows)

    # Write training-ready CSV
    out_fields = list(enriched[0].keys())
    with open(OUT_PATH, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=out_fields)
        writer.writeheader()
        writer.writerows(enriched)
    print(f"Written training-ready CSV: {OUT_PATH}")

    # Write endpoint summary
    summary = write_summary(enriched, SUM_PATH)
    print(f"\nEndpoint positive rates across {len(enriched)} PMO drugs:")
    print(f"{'Endpoint':<42} {'Pos':>5} {'Rate':>6}")
    print("-" * 56)
    for s in summary:
        print(f"{s['endpoint']:<42} {s['n_positive']:>5} {s['positive_rate']:>6.1%}")
    print(f"\nSummary CSV: {SUM_PATH}")


if __name__ == "__main__":
    main()
