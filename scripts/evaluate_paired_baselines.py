"""Is the model better than the best simple baseline, tested as a difference?

This is the gate's simple-baseline check. `derive_evidence_tiers.py` reads
what it writes, and an endpoint missing from the output fails that check, so
this runs inside `run_pipeline.py` before the tiers are derived rather than
by hand.

What it replaced: a comparison of the model's marginal 95% interval lower
bound to the baseline's point estimate. The two models are scored on the
same rows with the same folds, so their errors are correlated, and that is
not a test of the difference between them. It is far more conservative than
one, and conservative in a direction nobody chose.

This resamples whole sequences, scores both models on each resample, and
takes the distribution of the difference. The check passes when that
distribution excludes zero.

The change was not a loosening. On the run that introduced it:

    hepatocellular_injury     +0.119 [+0.036, +0.194]  P 0.996   promoted
    offtarget_gene_silencing  +0.355 [-0.075, +0.688]  P 0.917   demoted

The second had passed the marginal check and was being presented, on 12
positive sequences. A test that moves endpoints in both directions is a
different gate, not a weaker one.

    python scripts/evaluate_paired_baselines.py

Writes reports/results/paired_baselines.json.
"""
from __future__ import annotations

import json
import sys
import zlib
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score
from sklearn.model_selection import StratifiedGroupKFold, cross_val_predict

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from oligotox.training.train import _make_model, _prepare_training_matrix  # noqa: E402

PROCESSED = ROOT / "data" / "processed"
RESULTS = ROOT / "reports" / "results"
OUT = RESULTS / "paired_baselines.json"
N_BOOT = 1500
SIMPLE = {"GC content only": ["gc_content"], "length only": ["seq_length"]}


def main() -> int:
    X_all = pd.read_csv(PROCESSED / "feature_matrix.csv")
    Y = pd.read_csv(PROCESSED / "label_matrix.csv")
    raw = pd.read_csv(PROCESSED / "training_dataset.csv", low_memory=False)
    summary = json.loads((RESULTS / "training_summary_after.json").read_text())

    X_all = X_all.drop(
        columns=[c for c in X_all.columns if c.startswith(("label_", "endpoint_"))],
        errors="ignore")
    all_groups = (raw["sequence"].astype(str).str.upper()
                  .str.replace("U", "T", regex=False).values)

    report = {}
    print(f"{'endpoint':26s} {'model':>7s} {'simple':>7s} "
          f"{'paired difference':>24s} {'P':>6s}")
    for endpoint, spec in summary.items():
        if spec.get("status") != "ok" or endpoint not in Y.columns:
            continue
        keep = pd.to_numeric(Y[endpoint], errors="coerce").notna().values
        y = pd.to_numeric(Y[endpoint][keep], errors="coerce").astype(int).values
        groups = all_groups[keep]
        X = _prepare_training_matrix(X_all[keep].reset_index(drop=True))
        n_splits = spec["n_splits"]

        def oof(matrix):
            cv = StratifiedGroupKFold(n_splits=n_splits, shuffle=True,
                                      random_state=42)
            return cross_val_predict(_make_model(), matrix, y,
                                     cv=cv.split(matrix, y, groups),
                                     method="predict_proba")[:, 1]

        model = oof(X)
        scored = {name: oof(X[cols]) for name, cols in SIMPLE.items()}
        best = max(scored, key=lambda k: average_precision_score(y, scored[k]))
        baseline = scored[best]

        uniq = np.unique(groups)
        idx_by = {g: np.flatnonzero(groups == g) for g in uniq}
        # Seeded per endpoint, so this interval does not depend on what else
        # was measured first. See NEGATIVE_RESULTS entry 21.
        rng = np.random.default_rng([42, zlib.crc32(endpoint.encode())])
        diffs = []
        for _ in range(N_BOOT):
            drawn = rng.choice(uniq, size=len(uniq), replace=True)
            idx = np.concatenate([idx_by[g] for g in drawn])
            if len(np.unique(y[idx])) < 2:
                continue
            diffs.append(average_precision_score(y[idx], model[idx])
                         - average_precision_score(y[idx], baseline[idx]))
        diffs = np.array(diffs)
        lo, hi = (float(v) for v in np.percentile(diffs, [2.5, 97.5]))
        report[endpoint] = {
            "model_pr_auc": float(average_precision_score(y, model)),
            "best_simple": best,
            "simple_pr_auc": float(average_precision_score(y, baseline)),
            "difference": float(diffs.mean()),
            "ci": [lo, hi],
            "p_better": float((diffs > 0).mean()),
            "clears_paired": bool(lo > 0),
        }
        r = report[endpoint]
        print(f"{endpoint:26s} {r['model_pr_auc']:7.3f} {r['simple_pr_auc']:7.3f}"
              f"  {r['difference']:+.3f} [{lo:+.3f},{hi:+.3f}] {r['p_better']:6.3f}"
              + ("  PASS" if r["clears_paired"] else ""))

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps({
        "question": "Does the model beat the best simple baseline, tested as "
                    "a paired difference on the same resampled sequences?",
        "why": "The gate compares the model's marginal CI lower bound to the "
               "baseline's point estimate, which is not a test of the "
               "difference and is far more conservative than one.",
        "status": "this is the gate's simple-baseline check; "
                  "derive_evidence_tiers.py reads it, and an endpoint "
                  "missing from here fails that check.",
        "endpoints": report,
    }, indent=2) + "\n")
    print(f"\nwrote {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
