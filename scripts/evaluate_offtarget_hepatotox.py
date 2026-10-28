"""Does off-target complementarity burden predict hepatotoxicity?

This is the one feature with a published mechanism behind it. Burel et al.
and Kasuya et al. showed gapmer hepatotoxicity is RNase H1 cleaving
unintended transcripts the oligo is complementary to, particularly long
pre-mRNAs. If that mechanism is right and the burden is computable, it
should predict.

Hepatocellular injury is the endpoint that most needs it: it ranks within
its own source, does not calibrate, and scores 0.40 on an independent
screen. A mechanistic feature is a different kind of help from more rows.

Evaluated the same way as everything else: sequence-grouped folds, against
the same features without the off-target block, so the question is whether
it adds anything rather than whether it correlates.

**The background is whatever transcripts are cached**, which is a few
hundred genes rather than the transcriptome. Counts are therefore relative
between oligos, not absolute off-target risk, and a thin background
understates burden for every oligo roughly equally. If the feature helps on
a partial background it is worth computing properly; if it does nothing, a
larger background is unlikely to rescue it.

    python scripts/evaluate_offtarget_hepatotox.py

Writes reports/results/offtarget_hepatotox.csv and .json.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score
from sklearn.model_selection import StratifiedGroupKFold, cross_val_predict

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from oligotox.features.offtarget import burden, load_background  # noqa: E402
from oligotox.training.train import _make_model  # noqa: E402

PROCESSED = ROOT / "data" / "processed"
TRANSCRIPTS = ROOT / "data" / "external" / "transcripts"
RESULTS = ROOT / "reports" / "results"
ENDPOINTS = ["hepatocellular_injury", "offtarget_gene_silencing", "kupffer_cell_activation"]
N_BOOT = 1000
SEED = 42


def _oof(X: pd.DataFrame, y: np.ndarray, groups: np.ndarray, n_splits: int) -> np.ndarray:
    cv = StratifiedGroupKFold(n_splits=n_splits)
    return cross_val_predict(_make_model(), X, y, cv=cv, groups=groups,
                             method="predict_proba")[:, 1]


def _boot(y, p, groups, rng, n=N_BOOT):
    unique = np.unique(groups)
    scores = []
    for _ in range(n):
        picked = rng.choice(unique, size=len(unique), replace=True)
        idx = np.concatenate([np.flatnonzero(groups == g) for g in picked])
        if len(np.unique(y[idx])) < 2:
            continue
        scores.append(average_precision_score(y[idx], p[idx]))
    if not scores:
        return (float("nan"), float("nan"))
    return tuple(float(v) for v in np.percentile(scores, [2.5, 97.5]))


def main() -> None:
    # Prefer whole-transcriptome counts when they exist. A background of a few
    # cached target genes is 0.1 Mb against a transcriptome of 1,280 Mb, and a
    # null measured there says nothing about the mechanism.
    precomputed = PROCESSED / "offtarget_transcriptome.csv"
    use_precomputed = precomputed.exists()
    if not use_precomputed:
        try:
            background, n_transcripts, megabases = load_background(TRANSCRIPTS)
        except (FileNotFoundError, ValueError) as exc:
            raise SystemExit(
                f"{exc}\nRun scripts/offtarget_transcriptome.py for a real background, "
                "or scripts/fetch_transcripts.py for a partial one."
            )

    X = pd.read_csv(PROCESSED / "feature_matrix.csv")
    Y = pd.read_csv(PROCESSED / "label_matrix.csv")
    raw = pd.read_csv(PROCESSED / "training_dataset.csv", low_memory=False)
    summary = json.loads((RESULTS / "training_summary_after.json").read_text())
    groups = raw["sequence"].astype(str).str.upper().str.replace("U", "T", regex=False).values


    unique = pd.unique(groups)
    if use_precomputed:
        table = pd.read_csv(precomputed).set_index("sequence")
        print(f"background: whole transcriptome, "
              f"{json.loads((RESULTS / 'offtarget_transcriptome.json').read_text())['megabases_scanned']:,} Mb")
        feature_cols = [c for c in table.columns if c.startswith("offtarget_")]
        blank = {c: 0.0 for c in feature_cols}
        cache = {s: table.loc[s, feature_cols].to_dict() if s in table.index else blank
                 for s in unique}
        missing = sum(1 for s in unique if s not in table.index)
        if missing:
            print(f"  {missing} of {len(unique)} sequences had no transcriptome count")
    else:
        print(f"background: {n_transcripts} cached transcripts, {megabases:.2f} Mb "
              "(partial; counts are relative only)")
        cache = {s: burden(s, background, n_transcripts, megabases).as_features()
                 for s in unique}
    off = pd.DataFrame([cache[s] for s in groups], index=X.index).astype(float)
    off_cols = list(off.columns)
    print(f"off-target features: {off_cols}")
    spread = off[[c for c in off.columns if "seed" in c][0]]
    print(f"seed hits per Mb: min {spread.min():.1f} median {spread.median():.1f} "
          f"max {spread.max():.1f}\n")

    rng = np.random.default_rng(SEED)
    rows = []
    for endpoint in ENDPOINTS:
        stats = summary.get(endpoint, {})
        if stats.get("status") != "ok":
            continue
        y = pd.to_numeric(Y[endpoint], errors="coerce").fillna(0).astype(int).values
        n_splits = stats["n_splits"]

        base = _oof(X, y, groups, n_splits)
        with_off = _oof(pd.concat([X, off], axis=1), y, groups, n_splits)
        only_off = _oof(off, y, groups, n_splits)

        base_pr = float(average_precision_score(y, base))
        with_pr = float(average_precision_score(y, with_off))
        only_pr = float(average_precision_score(y, only_off))
        lo, hi = _boot(y, with_off, groups, rng)

        rows.append({
            "endpoint": endpoint,
            "base_rate": round(float(y.mean()), 4),
            "design_only": round(base_pr, 4),
            "design_plus_offtarget": round(with_pr, 4),
            "offtarget_only": round(only_pr, 4),
            "gain": round(with_pr - base_pr, 4),
            "ci_low": round(lo, 4), "ci_high": round(hi, 4),
        })
        print(f"  {endpoint:26s} design={base_pr:.3f}  +offtarget={with_pr:.3f}  "
              f"offtarget_only={only_pr:.3f}  gain={with_pr - base_pr:+.3f}")

    best = max(rows, key=lambda r: r["gain"]) if rows else None
    verdict = ("off-target burden adds signal" if best and best["gain"] > 0.01
               else "off-target burden adds nothing on this background")
    print(f"\n  {verdict}")

    RESULTS.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(RESULTS / "offtarget_hepatotox.csv", index=False)
    (RESULTS / "offtarget_hepatotox.json").write_text(json.dumps({
        "background": ("whole human transcriptome (Ensembl cDNA)" if use_precomputed
                       else "partial: cached target transcripts only"),
        "caveat": ("Counts are absolute occurrences across Ensembl cDNA."
                   if use_precomputed else
                   "Partial background; counts are relative between oligos only."),
        "mechanism": "Burel et al. and Kasuya et al.: gapmer hepatotoxicity is RNase H1 "
                     "cleavage of complementary unintended transcripts.",
        "verdict": verdict,
        "results": rows,
    }, indent=2) + "\n")
    print(f"wrote {(RESULTS / 'offtarget_hepatotox.csv').relative_to(ROOT)}")


if __name__ == "__main__":
    main()
