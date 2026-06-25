from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
EVAL_DIR = ROOT / "artifacts" / "eval"
OUT_PATH = ROOT / "artifacts" / "config" / "endpoint_reliability.json"

REPORT = EVAL_DIR / "train_eval_report.csv"
COMPARISON = EVAL_DIR / "train_eval_comparison.csv"


def classify_reliability(row_random, row_grouped):
  """
  Heuristic reliability tier:
    - HIGH: enough positives, enough positive groups, and grouped PR/ROC not terrible
    - MEDIUM: decent counts but grouped metrics degraded / undefined
    - LOW: very low support or grouped metrics unusable
  """
  total_pos = row_random["total_positives"]
  total_neg = row_random["total_negatives"]
  total_groups = row_random["total_groups"]
  groups_with_pos = row_random["groups_with_positives"]

  pr_random = row_random["pr_auc"]
  roc_random = row_random["roc_auc"]
  brier_random = row_random["brier"]

  pr_group = row_grouped["pr_auc"]
  roc_group = row_grouped["roc_auc"]
  brier_group = row_grouped["brier"]

  # Basic support rules
  if total_pos < 5 or total_neg < 5:
    return "LOW", "Very few positives/negatives; treat as exploratory only."

  if groups_with_pos <= 2:
    return "LOW", "Positives confined to very few groups; generalization is highly uncertain."

  # If grouped metrics are missing entirely → fragile
  if np.isnan(pr_group) or np.isnan(roc_group) or np.isnan(brier_group):
    return "MEDIUM", "Reasonable label support but grouped metrics undefined; use with caution."

  # Degradation thresholds
  delta_pr = pr_group - pr_random
  delta_roc = roc_group - roc_random
  delta_brier = brier_group - brier_random

  # Thresholds tuned for your scale
  if pr_group >= 0.6 and roc_group >= 0.7 and brier_group <= 0.30 and delta_pr > -0.2:
    return "HIGH", "Grouped performance acceptable relative to random; model generalizes across groups."
  else:
    return "MEDIUM", "Model has signal but grouped performance is lower; treat as moderate reliability."


def main():
  report = pd.read_csv(REPORT)
  comp = pd.read_csv(COMPARISON)

  rel = {}
  for endpoint in sorted(report["endpoint"].unique()):
    row_r = report[(report["endpoint"] == endpoint) & (report["split_type"] == "random")].iloc[0]
    row_g = report[(report["endpoint"] == endpoint) & (report["split_type"] == "grouped")].iloc[0]

    level, note = classify_reliability(row_r, row_g)

    rel[endpoint] = {
      "reliability_level": level,
      "note": note,
      "support": {
        "total_positives": int(row_r["total_positives"]),
        "total_negatives": int(row_r["total_negatives"]),
        "total_groups": int(row_r["total_groups"]),
        "groups_with_positives": int(row_r["groups_with_positives"]),
      },
      "metrics_random": {
        "pr_auc": None if np.isnan(row_r["pr_auc"]) else float(row_r["pr_auc"]),
        "roc_auc": None if np.isnan(row_r["roc_auc"]) else float(row_r["roc_auc"]),
        "brier": None if np.isnan(row_r["brier"]) else float(row_r["brier"]),
      },
      "metrics_grouped": {
        "pr_auc": None if np.isnan(row_g["pr_auc"]) else float(row_g["pr_auc"]),
        "roc_auc": None if np.isnan(row_g["roc_auc"]) else float(row_g["roc_auc"]),
        "brier": None if np.isnan(row_g["brier"]) else float(row_g["brier"]),
      },
    }

  OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
  OUT_PATH.write_text(json.dumps(rel, indent=2))
  print("Saved endpoint reliability to:", OUT_PATH)


if __name__ == "__main__":
  main()
