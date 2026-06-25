"""
Integration 5: DataSAIL Leakage-Free Grouped Splits
====================================================
Three strategies in order of rigor:
  1. scaffold (seed-based grouping, sklearn)  — always available
  2. group (target gene / drug series, sklearn) — always available
  3. datasail (cluster-based, pip install datasail) — gold standard

Ref: Gorgulla et al., DataSAIL, Bioinformatics 2024; DOI:10.1093/bioinformatics/btae447
     PMC 2025: https://pmc.ncbi.nlm.nih.gov/articles/PMC11217521/
"""
from __future__ import annotations
import numpy as np
import logging
from typing import List, Optional, Tuple
from difflib import SequenceMatcher

logger = logging.getLogger(__name__)


def scaffold_sequence_split(
    sequences: List[str], y: np.ndarray,
    test_size: float = 0.2, val_size: float = 0.1,
    seed_len: int = 7, random_state: int = 42,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Group by seed region (positions 2-8) — prevents within-seed-family leakage."""
    np.random.seed(random_state)

    def get_seed(seq):
        seq = seq.upper().replace("U", "T")
        return seq[1:1 + seed_len] if len(seq) >= 1 + seed_len else seq

    seed_groups = [get_seed(s) for s in sequences]
    unique_seeds = list(set(seed_groups))
    np.random.shuffle(unique_seeds)

    n = len(unique_seeds)
    n_test = max(1, int(n * test_size))
    n_val  = max(1, int(n * val_size))

    test_seeds = set(unique_seeds[:n_test])
    val_seeds  = set(unique_seeds[n_test:n_test + n_val])

    train_idx = np.array([i for i, s in enumerate(seed_groups) if s not in test_seeds and s not in val_seeds])
    val_idx   = np.array([i for i, s in enumerate(seed_groups) if s in val_seeds])
    test_idx  = np.array([i for i, s in enumerate(seed_groups) if s in test_seeds])

    logger.info(f"ScaffoldSplit ({len(unique_seeds)} seed groups): train={len(train_idx)}, val={len(val_idx)}, test={len(test_idx)}")
    return train_idx, val_idx, test_idx


def group_shuffle_split(
    X: np.ndarray, y: np.ndarray, groups: List[str],
    test_size: float = 0.2, val_size: float = 0.1,
    random_state: int = 42,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Group by target_gene or drug_name — sklearn GroupShuffleSplit."""
    from sklearn.model_selection import GroupShuffleSplit
    gss_test = GroupShuffleSplit(1, test_size=test_size, random_state=random_state)
    gss_val  = GroupShuffleSplit(1, test_size=val_size / (1 - test_size), random_state=random_state + 1)

    tv_idx, test_idx = next(gss_test.split(X, y, groups=groups))
    X_tv, y_tv = X[tv_idx], y[tv_idx]
    groups_tv = [groups[i] for i in tv_idx]

    tr_rel, va_rel = next(gss_val.split(X_tv, y_tv, groups=groups_tv))
    logger.info(f"GroupShuffleSplit: train={len(tr_rel)}, val={len(va_rel)}, test={len(test_idx)}")
    return tv_idx[tr_rel], tv_idx[va_rel], test_idx


def datasail_split(
    sequences: List[str], y: np.ndarray,
    test_size: float = 0.2, val_size: float = 0.1,
    similarity_threshold: float = 0.6, technique: str = "I1",
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    DataSAIL full leakage-reducing split. Falls back to scaffold split.
    Install: pip install datasail
    """
    try:
        from datasail.sail import datasail
        assignments = datasail(
            techniques=[technique],
            splits=[1 - test_size - val_size, val_size, test_size],
            names=["train", "val", "test"],
            e_type="S",
            e_data={str(i): s for i, s in enumerate(sequences)},
            similarity=similarity_threshold,
        )
        result = assignments[technique]
        train_idx = np.array([i for i, (_, split) in enumerate(result.items()) if split == "train"])
        val_idx   = np.array([i for i, (_, split) in enumerate(result.items()) if split == "val"])
        test_idx  = np.array([i for i, (_, split) in enumerate(result.items()) if split == "test"])
        logger.info(f"DataSAIL ({technique}): train={len(train_idx)}, val={len(val_idx)}, test={len(test_idx)}")
        return train_idx, val_idx, test_idx
    except ImportError:
        logger.warning("DataSAIL not installed (pip install datasail). Using scaffold split.")
        return scaffold_sequence_split(sequences, y, test_size=test_size, val_size=val_size)
    except Exception as e:
        logger.warning(f"DataSAIL failed ({e}). Using scaffold split.")
        return scaffold_sequence_split(sequences, y, test_size=test_size, val_size=val_size)


def check_leakage(
    train_seqs: List[str], test_seqs: List[str],
    seed_len: int = 7, sim_threshold: float = 0.8,
) -> dict:
    """Detect seed-level and sequence-level leakage between splits."""
    def seed(seq):
        seq = seq.upper().replace("U", "T")
        return seq[1:1 + seed_len] if len(seq) >= 1 + seed_len else seq

    train_seeds = set(seed(s) for s in train_seqs)
    test_seeds  = set(seed(s) for s in test_seqs)
    overlap = train_seeds & test_seeds

    leakage_pairs = []
    for ts in test_seqs[:50]:
        for tr in train_seqs[:50]:
            sim = SequenceMatcher(None, ts, tr).ratio()
            if sim >= sim_threshold:
                leakage_pairs.append((ts[:10], tr[:10], sim))

    report = {
        "train_size": len(train_seqs),
        "test_size":  len(test_seqs),
        "shared_seeds": len(overlap),
        "shared_seed_fraction": round(len(overlap) / max(len(test_seeds), 1), 3),
        "high_sim_pairs": len(leakage_pairs),
        "leakage_risk": "HIGH" if (len(overlap) > 0 or len(leakage_pairs) > 0) else "LOW",
    }
    if report["leakage_risk"] == "HIGH":
        logger.warning(
            f"LEAKAGE: {len(overlap)} shared seeds, {len(leakage_pairs)} high-sim pairs. "
            f"Use datasail_split() for publication-grade evaluation."
        )
    else:
        logger.info("Leakage check passed: LOW risk.")
    return report


def evaluate_endpoint(y_true, y_prob, endpoint_name="endpoint") -> dict:
    from sklearn.metrics import roc_auc_score, average_precision_score, brier_score_loss
    if len(np.unique(y_true)) < 2:
        return {"endpoint": endpoint_name, "pr_auc": float("nan"),
                "roc_auc": float("nan"), "brier": float("nan")}
    return {
        "endpoint": endpoint_name,
        "pr_auc":  float(average_precision_score(y_true, y_prob)),
        "roc_auc": float(roc_auc_score(y_true, y_prob)),
        "brier":   float(brier_score_loss(y_true, y_prob)),
        "n_pos":   int(y_true.sum()),
        "n_total": int(len(y_true)),
    }


if __name__ == "__main__":
    seqs = ["ACGUACGUACGUACGUACGU","GCGCGCGCGCGCGCGCGCGC",
            "UUUUACGUACGUACGUACGU","TTTTCGTTCGTTCGTTCGTT",
            "AAAAACGUACGUACGUACGU","CCCCGCGCGCGCGCGCGCGC"]
    y = np.array([1, 0, 1, 0, 1, 0])
    tr, va, te = scaffold_sequence_split(seqs, y)
    print(f"train={tr}, val={va}, test={te}")
    report = check_leakage([seqs[i] for i in tr], [seqs[i] for i in te])
    print(f"Leakage: {report}")
