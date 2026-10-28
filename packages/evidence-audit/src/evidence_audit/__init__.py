"""Check whether a model's score on a dataset means what it appears to mean.

A PR-AUC is a number you can always compute and rarely interpret. These
checks ask the questions it hides, and each one exists because it caught a
real false result (recorded in OligoTox's docs/NEGATIVE_RESULTS.md):

    leakage          is a feature a relabelled copy of the target?
    grouping         how many independent subjects are behind the row count?
    source_baseline  does a model seeing only which dataset a row came from
                     do just as well?
    within_source    does anything survive inside a single source?
    simple_baselines does one crude number match the whole feature matrix?
    calibration      are the probabilities probabilities, or only ranks?

Nothing here is specific to oligonucleotides, or to biology. The inputs are a
feature frame, a binary label, a grouping key per row, and optionally a
source per row. Anywhere rows repeat a subject and arrive from several
studies, the same six questions apply.

    from evidence_audit import audit_dataset
    report = audit_dataset(X, y, groups, sources)

or:

    evidence-audit data.csv --label outcome --group subject_id --source study
"""
from .checks import AuditReport, EndpointAudit, audit_dataset, audit_endpoint

__all__ = ["audit_dataset", "audit_endpoint", "AuditReport", "EndpointAudit"]
__version__ = "0.1.0"
