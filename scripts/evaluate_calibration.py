"""Is a predicted 0.7 actually 70%?

Every toxicity model is a CalibratedClassifierCV, but calibration was never
reported. For a safety tool that matters as much as ranking: an analyst who
reads 0.30 as "roughly a one-in-three chance" is relying on calibration, not
discrimination, and the two can come apart badly at low base rates.

Reports, per endpoint, on the same sequence-grouped out-of-fold predictions
the published PR-AUCs come from:

  brier            mean squared error of the probabilities
  brier_baseline   the same for always predicting the base rate
  brier_skill      1 - brier/brier_baseline; <= 0 means no better than
                   quoting the base rate
  calibration_mae  mean |predicted - observed| over probability bins
  mean_predicted   average predicted probability
  observed_rate    the actual positive rate

A well calibrated model has mean_predicted close to observed_rate and a
small calibration_mae. CNS toxicity is included via its own external model
and feature list.

Writes reports/results/calibration.csv and calibration.png.
"""
from __future__ import annotations

import json
from pathlib import Path

import joblib
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from sklearn.metrics import brier_score_loss  # noqa: E402
from sklearn.model_selection import StratifiedGroupKFold, StratifiedKFold, cross_val_predict  # noqa: E402

from oligotox.training.train import _make_model  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
PROCESSED = ROOT / "data" / "processed"
RESULTS = ROOT / "reports" / "results"
CONFIG = ROOT / "artifacts" / "config"
N_BINS = 5
SEED = 42


def _calibration_mae(y, p, n_bins=N_BINS):
    """Mean |predicted - observed| over equal-count bins, weighted by size."""
    order = np.argsort(p)
    bins = np.array_split(order, min(n_bins, max(1, len(p) // 20)))
    total, weight = 0.0, 0
    for idx in bins:
        if len(idx) == 0:
            continue
        total += abs(p[idx].mean() - y[idx].mean()) * len(idx)
        weight += len(idx)
    return total / weight if weight else float("nan")


def _row(endpoint, y, prob):
    base = float(y.mean())
    brier = float(brier_score_loss(y, prob))
    brier_base = float(np.mean((y - base) ** 2))
    return {
        "endpoint": endpoint,
        "n": int(len(y)),
        "positives": int(y.sum()),
        "observed_rate": base,
        "mean_predicted": float(prob.mean()),
        "brier": brier,
        "brier_baseline": brier_base,
        "brier_skill": float(1 - brier / brier_base) if brier_base > 0 else float("nan"),
        "calibration_mae": float(_calibration_mae(y, prob)),
    }


def internal_rows() -> tuple[list[dict], dict]:
    X = pd.read_csv(PROCESSED / "feature_matrix.csv")
    Y = pd.read_csv(PROCESSED / "label_matrix.csv")
    raw = pd.read_csv(PROCESSED / "training_dataset.csv", low_memory=False)
    summary = json.loads((RESULTS / "training_summary_after.json").read_text())
    groups = raw["sequence"].astype(str).str.upper().str.replace("U", "T", regex=False).values

    rows, curves = [], {}
    for endpoint, stats in summary.items():
        if stats["status"] != "ok":
            continue
        # Only the rows measured for this endpoint, matching the trainer.
        keep = pd.to_numeric(Y[endpoint], errors="coerce").notna().values
        y = pd.to_numeric(Y[endpoint][keep], errors="coerce").astype(int).values
        cv = StratifiedGroupKFold(n_splits=stats["n_splits"])
        prob = cross_val_predict(_make_model(), X[keep], y, cv=cv,
                                 groups=groups[keep],
                                 method="predict_proba")[:, 1]
        rows.append(_row(endpoint, y, prob))
        curves[endpoint] = (y, prob)
    return rows, curves


def cns_row() -> tuple[dict | None, tuple | None]:
    """CNS has its own external training set and feature list."""
    try:
        from evaluate_cns_external import HAGEDORN, MOE, _feature_cols, _model, build
    except ImportError:
        return None, None
    moe, hag = build(MOE), build(HAGEDORN)
    shared = sorted(set(_feature_cols(moe)) & set(_feature_cols(hag)))
    full = pd.concat([moe, hag], ignore_index=True)
    y = full["toxic"].values
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
    prob = cross_val_predict(_model(), full[shared].astype(float), y, cv=cv,
                             method="predict_proba")[:, 1]
    return _row("cns_toxicity", y, prob), (y, prob)


def plot(curves: dict, evidence: dict, out: Path) -> None:
    validated = [e for e in curves if evidence.get(e, {}).get("tier") == "validated"]
    if not validated:
        validated = list(curves)[:3]
    fig, ax = plt.subplots(figsize=(6.2, 5.2), dpi=160)
    ax.plot([0, 1], [0, 1], color="#999", lw=1, ls="--", label="perfect calibration")
    colors = ["#2a78d6", "#eb6834", "#2e9e5b", "#8457c4"]
    for (endpoint, color) in zip(validated, colors):
        y, p = curves[endpoint]
        order = np.argsort(p)
        bins = np.array_split(order, min(N_BINS, max(1, len(p) // 20)))
        xs = [p[i].mean() for i in bins if len(i)]
        ys = [y[i].mean() for i in bins if len(i)]
        ax.plot(xs, ys, marker="o", color=color, lw=1.6,
                label=endpoint.replace("_", " "))
    ax.set_xlabel("mean predicted probability")
    ax.set_ylabel("observed positive rate")
    ax.set_title("Calibration of the validated endpoints", loc="left")
    ax.legend(frameon=False, fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out)


def main() -> None:
    rows, curves = internal_rows()
    cns, cns_curve = cns_row()
    if cns:
        rows.append(cns)
        curves["cns_toxicity"] = cns_curve

    evidence = json.loads((CONFIG / "endpoint_evidence.json").read_text())["endpoints"]
    for row in rows:
        row["tier"] = evidence.get(row["endpoint"], {}).get("tier", "unknown")

    df = pd.DataFrame(rows).sort_values("brier_skill", ascending=False)
    RESULTS.mkdir(parents=True, exist_ok=True)
    df.to_csv(RESULTS / "calibration.csv", index=False)
    plot(curves, evidence, RESULTS / "calibration.png")

    cols = ["endpoint", "tier", "observed_rate", "mean_predicted", "brier",
            "brier_skill", "calibration_mae"]
    print(df[cols].to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    print(f"\nwrote {(RESULTS / 'calibration.csv').relative_to(ROOT)} and calibration.png")


if __name__ == "__main__":
    main()
