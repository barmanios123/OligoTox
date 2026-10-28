"""The checks must work on data that has nothing to do with any one project.

These import `evidence_audit` directly, never through OligoTox, so that a
dependency creeping back in fails here rather than at someone else's install.
"""
import numpy as np
import pandas as pd
import pytest

from evidence_audit import audit_dataset, audit_endpoint


@pytest.fixture
def synthetic():
    rng = np.random.default_rng(0)
    n_seq, reps = 60, 4
    seqs = [f"SEQ{i:03d}" for i in range(n_seq)]
    signal = rng.normal(size=n_seq)
    y_seq = (signal > 0.6).astype(int)
    rows = []
    for i, s in enumerate(seqs):
        for _ in range(reps):
            rows.append({
                "sequence": s,
                "gc_content": 0.4 + 0.02 * rng.normal(),
                "seq_length": 20.0,
                "real_feature": signal[i] + 0.3 * rng.normal(),
                "source": "study_a" if i < n_seq // 2 else "study_b",
                "y": y_seq[i],
            })
    return pd.DataFrame(rows)


def _features(df):
    return df.drop(columns=["sequence", "source", "y"])


def test_counts_sequences_not_rows(synthetic):
    a = audit_endpoint(_features(synthetic), synthetic["y"], synthetic["sequence"])
    assert a.rows == 240 and a.unique_groups == 60
    assert any("unique groups" in f for f in a.findings)


def test_grouping_changes_the_answer():
    """A per-sequence key is memorisable across replicates but generalises to nothing.

    Ungrouped folds let the model look it up; grouped folds do not. This is
    the hash-embedding failure in miniature.
    """
    rng = np.random.default_rng(7)
    n_seq, reps = 60, 4
    key = rng.normal(size=n_seq)          # unique per sequence, unrelated to y
    y_seq = rng.integers(0, 2, size=n_seq)
    rows = []
    for i in range(n_seq):
        for _ in range(reps):
            rows.append({
                "sequence": f"S{i:03d}",
                "identity_key": key[i],    # constant within a sequence
                "gc_content": 0.5 + 0.01 * rng.normal(),
                "seq_length": 20.0,
                "y": int(y_seq[i]),
            })
    df = pd.DataFrame(rows)
    a = audit_endpoint(df.drop(columns=["sequence", "y"]), df["y"], df["sequence"])
    # Memorisation shows up only when replicates span the split.
    assert a.ungrouped_pr_auc > a.grouped_pr_auc + 0.05


def test_a_leaked_label_is_caught(synthetic):
    df = synthetic.copy()
    df["curated_annotation"] = df["y"] * 1.0  # a relabelled copy of the target
    a = audit_endpoint(_features(df), df["y"], df["sequence"])
    assert "curated_annotation" in a.leaky_features
    assert any("annotations" in f for f in a.findings)


def test_source_confounding_is_caught():
    """When the source IS the label, the audit must say so."""
    rng = np.random.default_rng(1)
    rows = []
    for i in range(80):
        toxic = int(i >= 40)
        rows.append({
            "sequence": f"S{i:03d}",
            "gc_content": 0.5 + 0.01 * rng.normal(),
            "seq_length": 20.0,
            "noise": rng.normal(),
            # Positives live in exactly one source.
            "source": "curated" if toxic else "systematic",
            "y": toxic,
        })
    df = pd.DataFrame(rows)
    a = audit_endpoint(df.drop(columns=["sequence", "source", "y"]), df["y"],
                       df["sequence"], df["source"])
    assert a.source_only_pr_auc is not None
    assert a.verdict == "source confounded"


def test_too_few_positives_is_refused():
    rng = np.random.default_rng(2)
    df = pd.DataFrame({
        "sequence": [f"S{i}" for i in range(60)],
        "gc_content": rng.normal(0.5, 0.01, 60),
        "seq_length": 20.0,
        "y": [1] * 4 + [0] * 56,
    })
    a = audit_endpoint(df.drop(columns=["sequence", "y"]), df["y"], df["sequence"])
    assert a.verdict == "insufficient data"
    assert a.positive_groups == 4


def test_rna_and_dna_spellings_group_together():
    df = pd.DataFrame({
        "sequence": ["AUGC", "ATGC", "GGGG", "GGGG"],
        "gc_content": [0.5, 0.5, 1.0, 1.0],
        "y": [1, 1, 0, 0],
    })
    a = audit_endpoint(df.drop(columns=["sequence", "y"]), df["y"], df["sequence"],
                       nucleotide_groups=True)
    assert a.unique_groups == 2, "AUGC and ATGC are the same sequence"


def test_audit_dataset_covers_every_label(synthetic):
    df = synthetic.copy()
    df["y2"] = 1 - df["y"]
    report = audit_dataset(_features(df.drop(columns=["y2"])), df[["y", "y2"]],
                           df["sequence"], df["source"])
    assert {e.endpoint for e in report.endpoints} == {"y", "y2"}
    assert len(report.summary()) == 2


def test_group_keys_are_left_alone_unless_asked():
    """U is not T outside a nucleotide sequence.

    The folding that makes an RNA and a DNA spelling one group would turn the
    subject id "USER-1" into "TSER-1". It has to be opt-in, or the generic
    case is silently wrong.
    """
    rng = np.random.default_rng(3)
    rows = [{"group": g, "gc_content": rng.normal(), "seq_length": 20.0, "y": i % 2}
            for i, g in enumerate(["USER-1", "TSER-1"] * 20)]
    frame = pd.DataFrame(rows)
    features = frame[["gc_content", "seq_length"]]

    plain = audit_endpoint(features, frame["y"], frame["group"])
    assert plain.unique_groups == 2, "distinct ids must stay distinct"

    folded = audit_endpoint(features, frame["y"], frame["group"],
                            nucleotide_groups=True)
    assert folded.unique_groups == 1, "with folding on, U and T collapse"



def test_a_calibration_verdict_needs_a_margin_not_a_sign():
    """A threshold at exactly zero makes the verdict a coin flip on noise.

    The same dataset scored Brier skill +0.002 on one machine and -0.009 on
    another, which flipped it between "passes" and "ranks but uncalibrated"
    while nothing about the model had changed. The threshold has to sit far
    enough from zero that run-to-run variation cannot cross it.
    """
    from evidence_audit.checks import MIN_BRIER_SKILL, EndpointAudit, _verdict

    assert MIN_BRIER_SKILL > 0, "a sign test is not a calibration test"

    def verdict_at(skill: float) -> str:
        audit = EndpointAudit(
            endpoint="e", rows=400, unique_groups=400, positives=200,
            positive_groups=200, base_rate=0.5,
        )
        audit.grouped_pr_auc = 0.80
        audit.ungrouped_pr_auc = 0.80
        audit.brier_skill = skill
        return _verdict(audit)

    # Barely-positive skill is not calibration.
    assert verdict_at(0.002) == "ranks but uncalibrated"
    assert verdict_at(-0.009) == "ranks but uncalibrated"
    # Comfortably above the margin it passes.
    assert verdict_at(MIN_BRIER_SKILL + 0.1) == "passes"
