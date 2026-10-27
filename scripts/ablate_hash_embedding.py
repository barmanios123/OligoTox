"""Does the hash 'sequence embedding' contribute anything?

16 of the 55 features were `embed_stub_hash_rna_00..15`, a sha256 of the
sequence. The roadmap described the platform as having RNA foundation model
embeddings integrated; what reached the models was this hash.

A hash carries no structure a model can generalise from, but with 74 unique
sequences it can still act as an identity key, so it is worth measuring
rather than arguing about. This scores every evaluated endpoint with and
without the block, on the same sequence-grouped folds.

Result when first run: removing the block improved 9 of 12 endpoints, mean
+0.009 PR-AUC, off-target gene silencing by +0.073. The block is now dropped
in training (HASH_EMBEDDING_PREFIX in oligotox.training.train).

Needs a feature matrix that still contains the block, so run it against a
checkout from before that change.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score
from sklearn.model_selection import StratifiedGroupKFold, cross_val_predict

from oligotox.training.train import _make_model

ROOT = Path(__file__).resolve().parent.parent
PROCESSED = ROOT / "data" / "processed"
RESULTS = ROOT / "reports" / "results"
PREFIX = "embed_stub_hash"


def main() -> None:
    X = pd.read_csv(PROCESSED / "feature_matrix.csv")
    Y = pd.read_csv(PROCESSED / "label_matrix.csv")
    raw = pd.read_csv(PROCESSED / "training_dataset.csv", low_memory=False)
    summary = json.loads((RESULTS / "training_summary_after.json").read_text())

    hashed = [c for c in X.columns if c.startswith(PREFIX)]
    if not hashed:
        raise SystemExit(
            "No hash features in the matrix. Training already drops them, so run "
            "this against a checkout from before that change."
        )

    groups = raw["sequence"].astype(str).str.upper().str.replace("U", "T", regex=False).values
    kept = [c for c in X.columns if c not in hashed]

    rows = []
    for endpoint, stats in summary.items():
        if stats["status"] != "ok":
            continue
        y = pd.to_numeric(Y[endpoint], errors="coerce").fillna(0).astype(int).values
        cv = StratifiedGroupKFold(n_splits=stats["n_splits"])
        scored = {}
        for label, cols in (("with_hash", list(X.columns)), ("no_hash", kept)):
            prob = cross_val_predict(_make_model(), X[cols], y, cv=cv, groups=groups,
                                     method="predict_proba")[:, 1]
            scored[label] = float(average_precision_score(y, prob))
        rows.append({"endpoint": endpoint, **scored,
                     "delta_without_hash": scored["no_hash"] - scored["with_hash"]})

    df = pd.DataFrame(rows).sort_values("with_hash", ascending=False)
    df.to_csv(RESULTS / "hash_embedding_ablation.csv", index=False)
    delta = df["delta_without_hash"].values
    print(df.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    print(f"\nmean delta {delta.mean():+.4f}; better without the hash on "
          f"{int((delta > 0).sum())}/{len(delta)} endpoints")


if __name__ == "__main__":
    main()
