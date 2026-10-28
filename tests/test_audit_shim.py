"""OligoTox's audit import path must keep nucleotide grouping on.

The checks moved to the standalone evidence-audit package, where folding U
onto T is opt-in because it corrupts any non-sequence key. This project's
group keys are always oligo sequences, so the default has to be flipped back
here. If it is not, an RNA and a DNA spelling of one oligo count as two
independent sequences and every grouped fold in the repository silently
leaks a replicate across the split.
"""
import pandas as pd

from oligotox.audit import audit_dataset, audit_endpoint


def _frame():
    return pd.DataFrame({
        "sequence": ["AUGC", "ATGC", "GGGG", "GGGG"],
        "gc_content": [0.5, 0.5, 1.0, 1.0],
        "y": [1, 1, 0, 0],
    })


def test_the_shim_folds_rna_onto_dna_without_being_asked():
    df = _frame()
    audit = audit_endpoint(df[["gc_content"]], df["y"], df["sequence"])
    assert audit.unique_groups == 2, "AUGC and ATGC are one oligo"


def test_the_caller_can_still_turn_it_off():
    df = _frame()
    audit = audit_endpoint(df[["gc_content"]], df["y"], df["sequence"],
                           nucleotide_groups=False)
    assert audit.unique_groups == 3


def test_audit_dataset_folds_too():
    df = _frame()
    report = audit_dataset(df[["gc_content"]], df[["y"]], df["sequence"])
    assert report.endpoints[0].unique_groups == 2
