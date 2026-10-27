"""Plot model vs baseline PR-AUC with 95% bootstrap CIs.

Reads reports/results/baselines.csv and writes reports/results/baselines.png.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "reports" / "results"

SURFACE = "#fcfcfb"
TEXT = "#0b0b0b"
TEXT_2 = "#52514e"
GRID = "#e4e3df"
MODELS = [
    ("Logistic (", "#2a78d6", "o"),
    ("Gradient boosting (", "#eb6834", "D"),
    ("Data source only", "#1baf7a", "s"),
    ("GC content only", "#eda100", "^"),
    ("Length only", "#e87ba4", "v"),
]


def _resolve(df, prefix):
    """Model labels carry the feature count, so match on the prefix."""
    hits = [m for m in df.model.unique() if m.startswith(prefix)]
    return hits[0] if hits else prefix

def main() -> None:
    df = pd.read_csv(RESULTS / "baselines.csv")
    logistic = df[df.model == MODELS[0][0]].set_index("endpoint")
    endpoints = logistic.sort_values("positive_sequences").index.tolist()

    fig, ax = plt.subplots(figsize=(9, 0.7 * len(endpoints) + 1.8), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    y = np.arange(len(endpoints))
    step = 0.15
    for i, (name, color, marker) in enumerate(MODELS):
        sub = df[df.model == name].set_index("endpoint").loc[endpoints]
        yy = y + (2 - i) * step
        ax.hlines(yy, sub.ci_low, sub.ci_high, color=color, linewidth=2, zorder=2)
        ax.scatter(sub.pr_auc, yy, color=color, marker=marker, s=46, zorder=3,
                   edgecolors=SURFACE, linewidths=1.5, label=name)
    ax.scatter(logistic.loc[endpoints, "base_rate"], y, marker="|", s=420, color=TEXT,
               linewidths=2, zorder=4, label="Base rate (random guess)")

    labels = [
        f"{e.replace('_', ' ')}  (n+ seq = {int(logistic.loc[e, 'positive_sequences'])})"
        for e in endpoints
    ]
    ax.set_yticks(y, labels, color=TEXT, fontsize=9)
    ax.set_xlim(0, 1)
    ax.set_xlabel("PR-AUC, sequence-grouped CV (bars: 95% bootstrap CI over sequences)", color=TEXT_2)
    ax.tick_params(axis="x", colors=TEXT_2)
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="x", color=GRID, zorder=0)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    fig.legend(loc="upper left", bbox_to_anchor=(0.01, 0.955), ncol=3, frameon=False,
               fontsize=8, labelcolor=TEXT_2)
    fig.suptitle(
        "Model vs simple baselines: most intervals overlap, and data source alone wins on Kupffer",
        x=0.01, ha="left", color=TEXT, fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    out = RESULTS / "baselines.png"
    fig.savefig(out, dpi=160, facecolor=SURFACE)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
