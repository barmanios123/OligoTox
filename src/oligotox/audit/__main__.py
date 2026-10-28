"""`python -m oligotox.audit` runs the standalone auditor on sequence data.

Kept so the command the README documents keeps working. It forwards to
`evidence-audit` with --nucleotide set, since this project's group keys are
oligo sequences.
"""
from __future__ import annotations

import sys

from evidence_audit.__main__ import main as _main


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if "--nucleotide" not in argv:
        argv.append("--nucleotide")
    return _main(argv)


if __name__ == "__main__":
    sys.exit(main())
