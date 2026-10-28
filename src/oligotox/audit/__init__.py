"""OligoTox's view of the standalone `evidence-audit` package.

The checks used to live here. They turned out to need nothing from OligoTox
and nothing from biology, so they were moved to packages/evidence-audit to be
installable and usable on their own. This module stays as the import path the
repository's scripts and docs already use, and it sets the one default that
is specific to this project: group keys are nucleotide sequences, so U and T
spellings of one oligo count as one group.

    from oligotox.audit import audit_dataset
    report = audit_dataset(X, y, sequences, sources)

The generic tool is `evidence-audit` on the command line, or
`pip install packages/evidence-audit` to use it elsewhere.
"""
from __future__ import annotations

from typing import Any, Sequence

try:
    from evidence_audit import AuditReport, EndpointAudit
    from evidence_audit import audit_dataset as _audit_dataset
    from evidence_audit import audit_endpoint as _audit_endpoint
except ImportError as exc:  # pragma: no cover - the repo installs it
    raise ImportError(
        "oligotox.audit needs the evidence-audit package, which lives in this "
        "repository at packages/evidence-audit. Install it with "
        "`pip install -e packages/evidence-audit`."
    ) from exc


def audit_dataset(*args: Any, **kwargs: Any) -> AuditReport:
    """`evidence_audit.audit_dataset`, with nucleotide grouping on."""
    kwargs.setdefault("nucleotide_groups", True)
    return _audit_dataset(*args, **kwargs)


def audit_endpoint(*args: Any, **kwargs: Any) -> EndpointAudit:
    """`evidence_audit.audit_endpoint`, with nucleotide grouping on."""
    kwargs.setdefault("nucleotide_groups", True)
    return _audit_endpoint(*args, **kwargs)


__all__ = ["audit_dataset", "audit_endpoint", "AuditReport", "EndpointAudit"]
