"""Plot per-endpoint PR-AUC before and after the leakage fix.

Reads reports/results/training_summary_{before,after}.json and writes
reports/results/pr_auc_by_endpoint.png.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "reports" / "results"

SURFACE = "#fcfcfb"
TEXT = "#0b0b0b"
TEXT_2 = "#52514e"
GRID = "#e4e3df"
SERIES = [
    ("Leak-free features, sequence-grouped CV", "#2a78d6"),
    ("Leak-free features, row-level CV", "#eb6834"),
    ("Before fix (leaky features, row-level CV)", "#1baf7a"),
]


def main() -> None:
    before = json.loads((RESULTS / "training_summary_before.json").read_text())
    after = json.loads((RESULTS / "training_summary_after.json").read_text())

    endpoints = [e for e, v in after.items() if v.get("pr_auc") is not None]
    endpoints.sort(key=lambda e: after[e]["pr_auc"])

    values = [
        [after[e]["pr_auc"] for e in endpoints],
        [after[e].get("random_split_pr_auc") or 0.0 for e in endpoints],
        [before[e].get("pr_auc") or 0.0 for e in endpoints],
    ]
    base = [after[e]["positive_rate"] for e in endpoints]
    labels = [
        f"{e.replace('_', ' ')}  (n+ seq = {after[e]['positive_sequences']})"
        for e in endpoints
    ]

    fig, ax = plt.subplots(figsize=(9, 0.62 * len(endpoints) + 1.6), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    y = np.arange(len(endpoints))
    h = 0.26
    for i, ((name, color), vals) in enumerate(zip(SERIES, values)):
        ax.barh(y + (1 - i) * h, vals, height=h - 0.03, color=color, label=name, zorder=2)
    ax.scatter(base, y, marker="|", s=260, color=TEXT, linewidths=2, zorder=3,
               label="Base rate (random guess)")

    ax.set_yticks(y, labels, color=TEXT, fontsize=9)
    ax.set_xlim(0, 1)
    ax.set_xlabel("PR-AUC (out-of-fold)", color=TEXT_2)
    ax.tick_params(axis="x", colors=TEXT_2)
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="x", color=GRID, zorder=0)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    fig.legend(loc="upper left", bbox_to_anchor=(0.01, 0.955), ncol=2, frameon=False,
               fontsize=8, labelcolor=TEXT_2)
    fig.suptitle(
        "Per-endpoint PR-AUC: blue is the honest estimate for an unseen oligo sequence",
        x=0.01, ha="left", color=TEXT, fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    out = RESULTS / "pr_auc_by_endpoint.png"
    fig.savefig(out, dpi=160, facecolor=SURFACE)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
