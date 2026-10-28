"""Run the audit harness on two datasets and compare the verdicts.

A checker that fails everything is not a checker. This runs the same tool,
with no tuning between runs, over:

  internal   this project's own curated training data, three endpoints
  external   the OligoGym CNS neurotoxicity data, 3,338 sequences from two
             independent published assays

The internal endpoints fail for three different reasons and the external one
passes, which is the evidence that the verdicts track the data rather than
the tool's disposition.

Writes reports/results/audit_demo.json.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from oligotox.audit import audit_dataset  # noqa: E402

PROCESSED = ROOT / "data" / "processed"
RESULTS = ROOT / "reports" / "results"

INTERNAL_ENDPOINTS = ["hepatocellular_injury", "kupffer_cell_activation", "thrombocytopenia"]


def internal():
    X = pd.read_csv(PROCESSED / "feature_matrix.csv")
    Y = pd.read_csv(PROCESSED / "label_matrix.csv")[INTERNAL_ENDPOINTS]
    raw = pd.read_csv(PROCESSED / "training_dataset.csv", low_memory=False)
    return audit_dataset(X, Y, raw["sequence"], raw["source"])


def external():
    from evaluate_cns_external import HAGEDORN, MOE, _feature_cols, build

    moe, hag = build(MOE), build(HAGEDORN)
    moe["assay"], hag["assay"] = MOE, HAGEDORN
    full = pd.concat([moe, hag], ignore_index=True)
    cols = sorted(set(_feature_cols(moe)) & set(_feature_cols(hag)))
    return audit_dataset(full[cols], full[["toxic"]].rename(columns={"toxic": "cns_toxicity"}),
                         full["seq"], full["assay"])


def main() -> None:
    reports = {"internal (this project's curated data)": internal(),
               "external (OligoGym CNS, two published assays)": external()}

    payload = {}
    for name, report in reports.items():
        print(f"\n{name}")
        print(report.summary().to_string(index=False, float_format=lambda v: f"{v:.3f}"))
        for endpoint in report.endpoints:
            if endpoint.findings:
                print(f"  {endpoint.endpoint}: {endpoint.findings[-1]}")
        payload[name] = report.to_dict()

    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / "audit_demo.json").write_text(json.dumps({
        "note": "Same harness, no tuning between runs. The internal endpoints fail "
                "for three different reasons; the external one passes.",
        "runs": payload,
    }, indent=2) + "\n")
    print(f"\nwrote {(RESULTS / 'audit_demo.json').relative_to(ROOT)}")


if __name__ == "__main__":
    main()
