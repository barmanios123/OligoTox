#!/usr/bin/env python3
"""train_v4_fast.py

Canonical v4 training script. Trains lightweight RF models with joblib serialization.
All models saved to artifacts/models_v4/ as .joblib files.
This script replaces all legacy pkl artifacts.

Usage:
    cd ~/oligotox
    PYTHONPATH="$HOME:$PWD:$PWD/src" python3 scripts/train_v4_fast.py
"""
from __future__ import annotations
import glob
import os
import sys
from itertools import product
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

REPO = Path(__file__).resolve().parents[1]
MODEL_DIR = REPO / "artifacts" / "models_v4"
REPORT_DIR = REPO / "reports"
MODEL_DIR.mkdir(parents=True, exist_ok=True)
REPORT_DIR.mkdir(parents=True, exist_ok=True)

# ── Hyperparameters ──────────────────────────────────────────────────────────
RF_PARAMS = dict(
    n_estimators=100,
    max_depth=15,
    max_features="sqrt",
    min_samples_leaf=2,
    n_jobs=1,           # single-threaded avoids macOS fork issues
    random_state=42,
    class_weight="balanced",
)
MIN_POS = 5            # minimum positives to attempt training
CV_FOLDS = 5


# ── Feature engineering ──────────────────────────────────────────────────────
def _kmer(seq: str, k: int) -> dict:
    nts = "ACGTU"
    kmers = ["".join(p) for p in product(nts, repeat=k)]
    n = len(seq)
    return {f"k{k}_{km}": seq.count(km) / max(n - k + 1, 1) for km in kmers}


def featurise(seq: str) -> dict:
    seq = seq.upper().replace(" ", "").replace("-", "")
    n = len(seq)
    gc = (seq.count("G") + seq.count("C")) / max(n, 1)
    t = max(n // 3, 1)
    f = {
        "length_nt": float(n),
        "gc_content": gc,
        "gc_5p": (seq[:t].count("G") + seq[:t].count("C")) / max(t, 1),
        "gc_mid": (seq[t:2*t].count("G") + seq[t:2*t].count("C")) / max(t, 1),
        "gc_3p": (seq[2*t:].count("G") + seq[2*t:].count("C")) / max(len(seq[2*t:]), 1),
    }
    f.update(_kmer(seq, 2))
    f.update(_kmer(seq, 3))
    return f


def featurise_df(seqs: pd.Series) -> pd.DataFrame:
    return pd.DataFrame([featurise(s) for s in seqs.fillna("")]).fillna(0.0)


# ── Load training data ────────────────────────────────────────────────────────
def load_data() -> pd.DataFrame:
    candidates = [
        REPO / "data" / "curated" / "curated_v4.csv",
        REPO / "data" / "curated" / "curated_v3.csv",
        REPO / "data" / "processed" / "training_data_with_features.csv",
    ]
    for p in candidates:
        if p.exists():
            print(f"Loading {p}")
            return pd.read_csv(p)
    raise FileNotFoundError(f"No training data found. Checked: {candidates}")


# ── Train one endpoint ─────────────────────────────────────────────────────────
def train_endpoint(X: pd.DataFrame, y: pd.Series, ep: str) -> dict:
    n_pos = int(y.sum())
    n_neg = int((y == 0).sum())
    if n_pos < MIN_POS:
        return {"endpoint": ep, "n": len(y), "n_pos": n_pos,
                "n_neg": n_neg, "auroc": None, "status": "skipped"}

    model = RandomForestClassifier(**RF_PARAMS)

    # CV evaluation
    skf = StratifiedKFold(n_splits=min(CV_FOLDS, n_pos), shuffle=True, random_state=42)
    oof_probs, oof_labels = [], []
    for tr, te in skf.split(X, y):
        model.fit(X.iloc[tr], y.iloc[tr])
        oof_probs.extend(model.predict_proba(X.iloc[te])[:, 1])
        oof_labels.extend(y.iloc[te])

    try:
        auroc = round(float(roc_auc_score(oof_labels, oof_probs)), 4)
    except Exception:
        auroc = None

    # Fit full model
    model.fit(X, y)

    # Save
    ep_short = ep.replace("endpoint_", "")
    out_path = MODEL_DIR / f"{ep_short}.joblib"
    joblib.dump(model, out_path, compress=3)

    return {"endpoint": ep, "n": len(y), "n_pos": n_pos, "n_neg": n_neg,
            "auroc": auroc, "status": "trained", "model_path": str(out_path)}


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    # Remove stale models
    for ext in ("*.pkl", "*.joblib"):
        for f in glob.glob(str(MODEL_DIR / ext)):
            os.remove(f)
            print(f"  Removed {f}")

    df = load_data()
    print(f"Loaded {len(df)} rows, columns: {list(df.columns[:8])}...")

    # Detect sequence column
    seq_col = next((c for c in df.columns if "sequence" in c.lower()), None)
    if seq_col is None:
        raise ValueError("No sequence column found")
    print(f"Sequence column: {seq_col}")

    # Detect endpoint columns
    ep_cols = [c for c in df.columns if c.startswith("endpoint_")
               and not c.endswith("_pmo") and not c.endswith("_score")]
    print(f"Endpoints ({len(ep_cols)}): {ep_cols}")

    # Featurise
    print("Featurising sequences...")
    X = featurise_df(df[seq_col])
    print(f"Feature matrix: {X.shape}")

    results = []
    for ep in ep_cols:
        y = pd.to_numeric(df[ep], errors="coerce").dropna()
        y = y[y.isin([0, 1])]
        X_ep = X.loc[y.index]
        res = train_endpoint(X_ep, y, ep)
        status = res.get("status", "?")
        auroc_str = f"AUROC={res['auroc']:.3f}" if res["auroc"] else "(skipped)"
        print(f"  {ep:<45s}  n_pos={res['n_pos']:4d}  {auroc_str}  [{status}]")
        results.append(res)

    res_df = pd.DataFrame(results)
    report_path = REPORT_DIR / "eval_v4_fast.csv"
    res_df.to_csv(report_path, index=False)
    print(f"\nSaved evaluation report: {report_path}")

    valid = res_df.dropna(subset=["auroc"])
    if len(valid):
        print(f"Mean AUROC: {valid['auroc'].mean():.3f}")
        print(f"Median AUROC: {valid['auroc'].median():.3f}")
        print(f"Best:  {valid.iloc[0]['endpoint']} = {valid['auroc'].max():.3f}")
        print(f"Worst: {valid.iloc[-1]['endpoint']} = {valid['auroc'].min():.3f}")

    print(f"\nModels written to: {MODEL_DIR}")
    for f in sorted(MODEL_DIR.glob("*.joblib")):
        print(f"  {f.name}  ({os.path.getsize(f)//1024}KB)")


if __name__ == "__main__":
    main()
