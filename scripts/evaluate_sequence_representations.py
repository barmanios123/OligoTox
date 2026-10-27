"""Do richer sequence representations beat the hand-built features?

The strategy roadmap lists three RNA foundation models (Nucleotide
Transformer, RiNALMo, BiRNA-BERT) under "what we have". Adapter files exist
under src/oligotox/embeddings, but nothing calls them, and what actually
reached the models was a sha256 hash that an ablation showed was worse than
nothing (scripts/ablate_hash_embedding.py).

Until recently there was no point testing a learned representation on 74
sequences. The CNS data changed that: 3,338 sequences across two independent
assays, with a cross-assay split that cannot be gamed.

This compares, on that split:

  Hand features      the 17 design features the shipped CNS model uses
  k-mer counts       normalised 3-mer and 4-mer frequencies, a real
                     representation of sequence content with no hand design
  k-mer + hand       both together
  Hash stub          the sha256 block, as a control that should do nothing

Pretrained RNA language model embeddings are not included: no weights are
available in this environment, and downloading them is a prerequisite rather
than a result. The k-mer arm answers the question those embeddings are meant
to answer, which is whether learned sequence content beats hand design here.

Writes reports/results/sequence_representations.csv.
"""
from __future__ import annotations

import hashlib
import json
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict

from evaluate_cns_external import HAGEDORN, MOE, SEED, _boot_auc, _feature_cols, _model, build

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "reports" / "results"


def kmer_features(seq: str, k: int) -> dict:
    seq = seq.strip().upper().replace("U", "T")
    counts = {"".join(p): 0 for p in product("ACGT", repeat=k)}
    total = max(len(seq) - k + 1, 1)
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        if kmer in counts:
            counts[kmer] += 1
    return {f"k{k}_{m}": c / total for m, c in counts.items()}


def hash_features(seq: str, dim: int = 16) -> dict:
    """The removed stub, kept here only as a control."""
    out = {}
    for i in range(dim):
        digest = hashlib.sha256(f"{seq}|{i}".encode()).hexdigest()
        out[f"hash_{i:02d}"] = int(digest[:8], 16) / 0xFFFFFFFF
    return out


def with_representations(df: pd.DataFrame) -> pd.DataFrame:
    seqs = df["seq"]
    k3 = pd.DataFrame([kmer_features(s, 3) for s in seqs], index=df.index)
    k4 = pd.DataFrame([kmer_features(s, 4) for s in seqs], index=df.index)
    hsh = pd.DataFrame([hash_features(s) for s in seqs], index=df.index)
    return pd.concat([df, k3, k4, hsh], axis=1)


def main() -> None:
    rng = np.random.default_rng(SEED)
    moe, hag = (with_representations(build(s)) for s in (MOE, HAGEDORN))

    hand = [c for c in _feature_cols(build(MOE))]
    k3 = [c for c in moe.columns if c.startswith("k3_")]
    k4 = [c for c in moe.columns if c.startswith("k4_")]
    hsh = [c for c in moe.columns if c.startswith("hash_")]

    sets = {
        "Hand features": hand,
        "3-mer counts": k3,
        "3-mer + 4-mer": k3 + k4,
        "k-mer + hand": k3 + k4 + hand,
        "Hash stub (control)": hsh,
    }

    rows = []
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
    for label, cols in sets.items():
        # Within-assay, on the larger set.
        y = moe["toxic"].values
        prob = cross_val_predict(_model(), moe[cols].astype(float), y, cv=cv,
                                 method="predict_proba")[:, 1]
        lo, hi = _boot_auc(y, prob, rng)
        rows.append({"evaluation": f"within {MOE} 5-fold", "representation": label,
                     "n_features": len(cols), "roc_auc": float(roc_auc_score(y, prob)),
                     "ci_low": lo, "ci_high": hi})

        # Cross-assay, both directions. This is the one that counts.
        for train, test, name in ((moe, hag, f"{MOE} -> {HAGEDORN}"),
                                  (hag, moe, f"{HAGEDORN} -> {MOE}")):
            model = _model().fit(train[cols].astype(float), train["toxic"].values)
            p = model.predict_proba(test[cols].astype(float))[:, 1]
            yt = test["toxic"].values
            lo, hi = _boot_auc(yt, p, rng)
            rows.append({"evaluation": f"cross-assay {name}", "representation": label,
                         "n_features": len(cols), "roc_auc": float(roc_auc_score(yt, p)),
                         "ci_low": lo, "ci_high": hi})

    df = pd.DataFrame(rows)
    RESULTS.mkdir(parents=True, exist_ok=True)
    df.to_csv(RESULTS / "sequence_representations.csv", index=False)
    (RESULTS / "sequence_representations.json").write_text(json.dumps({
        "note": "Pretrained RNA language model embeddings were not evaluated: no "
                "weights are available in this environment. The k-mer arms test "
                "whether learned sequence content beats hand-designed features.",
        "results": rows,
    }, indent=2) + "\n")

    for evaluation, sub in df.groupby("evaluation"):
        print(f"\n{evaluation}")
        print(sub[["representation", "n_features", "roc_auc", "ci_low", "ci_high"]]
              .to_string(index=False, float_format=lambda v: f"{v:.3f}"))


if __name__ == "__main__":
    main()
