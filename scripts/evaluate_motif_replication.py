"""Does the published gapmer hepatotoxicity motif replicate?

Kasuya et al. (NAR 2022;50(13):7224) report that TGC and TCC counts track
hepatotoxicity in LNA gapmers, and the Papargyri rows in this corpus agree
strongly. Until the patent ingestion there was no second in-vivo source to
check it against; there are now 256 sequences from 30 Ionis patents.

This is worth running rather than assuming, because a motif that holds in
one source and not another is a source effect wearing a mechanism's clothes,
and adding it as a feature would teach the model to recognise the source.

    python scripts/evaluate_motif_replication.py

Writes reports/results/motif_replication.json.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
PROCESSED = ROOT / "data" / "processed"
OUT = ROOT / "reports" / "results" / "motif_replication.json"
ENDPOINT = "endpoint_hepatocellular_injury"
MOTIFS = ("TGC", "TCC")
N_BOOT = 2000


def _counts(seq: pd.Series, motif: str) -> pd.Series:
    return seq.str.count(f"(?={motif})")


def _diff_ci(values: np.ndarray, labels: np.ndarray, groups: np.ndarray,
             rng: np.random.Generator) -> tuple[float, float, float]:
    """Positive-minus-negative mean, resampling whole sequences."""
    point = values[labels == 1].mean() - values[labels == 0].mean()
    uniq = np.unique(groups)
    idx_by = {g: np.flatnonzero(groups == g) for g in uniq}
    draws = []
    for _ in range(N_BOOT):
        drawn = rng.choice(uniq, size=len(uniq), replace=True)
        idx = np.concatenate([idx_by[g] for g in drawn])
        lab = labels[idx]
        if lab.min() == lab.max():
            continue
        draws.append(values[idx][lab == 1].mean() - values[idx][lab == 0].mean())
    lo, hi = np.percentile(draws, [2.5, 97.5])
    return float(point), float(lo), float(hi)


def main() -> int:
    raw = pd.read_csv(PROCESSED / "training_dataset.csv", low_memory=False)
    df = raw[raw[ENDPOINT].notna()].copy()
    df["seq"] = (df["sequence"].astype(str).str.upper()
                 .str.replace("U", "T", regex=False))
    labels = df[ENDPOINT].astype(int).values
    rng = np.random.default_rng(42)

    report: dict[str, dict] = {}
    for motif in MOTIFS:
        values = _counts(df["seq"], motif).values.astype(float)
        per_source = {}
        for source, sub in df.groupby("source"):
            if sub[ENDPOINT].nunique() < 2 or len(sub) < 20:
                continue
            v = _counts(sub["seq"], motif).values.astype(float)
            lab = sub[ENDPOINT].astype(int).values
            point, lo, hi = _diff_ci(v, lab, sub["seq"].values, rng)
            per_source[source] = {
                "sequences": int(sub["seq"].nunique()),
                "rows": int(len(sub)),
                "mean_in_positives": float(v[lab == 1].mean()),
                "mean_in_negatives": float(v[lab == 0].mean()),
                "difference": point,
                "ci": [lo, hi],
                "direction": ("higher in toxic" if point > 0
                              else "higher in non-toxic"),
                "clears_zero": bool(lo > 0 or hi < 0),
            }
        point, lo, hi = _diff_ci(values, labels, df["seq"].values, rng)
        directions = {s["direction"] for s in per_source.values()}
        report[motif] = {
            "pooled": {"difference": point, "ci": [lo, hi]},
            "per_source": per_source,
            "replicates": bool(len(directions) == 1
                               and all(s["clears_zero"]
                                       for s in per_source.values())),
        }

    summary = {
        "claim": "Kasuya et al., NAR 2022;50(13):7224: TGC and TCC counts "
                 "track hepatotoxicity in LNA gapmers.",
        "question": "Does it hold in a second in-vivo source?",
        "rule": "It replicates only if the direction agrees across sources "
                "and each source's own interval clears zero. A motif that "
                "holds in one source is a source effect, and adding it as a "
                "feature would teach the model to recognise that source.",
        "motifs": report,
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(summary, indent=2) + "\n")

    for motif, spec in report.items():
        print(f"\n{motif}: pooled {spec['pooled']['difference']:+.3f} "
              f"[{spec['pooled']['ci'][0]:+.3f}, {spec['pooled']['ci'][1]:+.3f}]"
              f"   replicates: {spec['replicates']}")
        for source, s in spec["per_source"].items():
            print(f"   {source:18s} {s['sequences']:4d} seq  "
                  f"{s['difference']:+.3f} [{s['ci'][0]:+.3f}, {s['ci'][1]:+.3f}]"
                  f"  {s['direction']}")
    print(f"\nwrote {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
