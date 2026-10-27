"""Audit any oligo toxicity dataset and model, not just this project's.

The checks in this package are the ones that caught every failure recorded in
docs/NEGATIVE_RESULTS.md. They need nothing from OligoTox: a feature matrix,
a label column, a sequence column to group on, and optionally a source
column. Point them at someone else's data and they answer the same questions.

    from oligotox.audit import audit_dataset
    report = audit_dataset(X, y, sequences, sources)

or from the command line:

    python -m oligotox.audit data.csv --label hepatotox --sequence seq --source study
"""
from .checks import AuditReport, EndpointAudit, audit_dataset, audit_endpoint

__all__ = ["audit_dataset", "audit_endpoint", "AuditReport", "EndpointAudit"]
