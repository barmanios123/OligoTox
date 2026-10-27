"""Audit someone else's oligo toxicity dataset from the command line.

    python -m oligotox.audit data.csv --label hepatotox --sequence seq

    python -m oligotox.audit data.csv \\
        --label tox_a --label tox_b \\
        --sequence sequence --source study --json out.json

Everything not named by --label, --sequence, --source or --exclude is treated
as a feature. The exit status is 1 if any endpoint fails a check, so this can
gate a pipeline.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

from .checks import audit_dataset

PASSING = {"passes"}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m oligotox.audit",
        description="Check an oligo toxicity dataset for leakage, replicate splits, "
                    "source confounding, trivial baselines and miscalibration.",
    )
    parser.add_argument("data", type=Path, help="CSV with features, labels and a sequence column")
    parser.add_argument("--label", action="append", required=True,
                        help="label column; repeat for several endpoints")
    parser.add_argument("--sequence", required=True,
                        help="column to group folds on, so replicates never span a split")
    parser.add_argument("--source", default=None,
                        help="column naming the originating dataset, for the source baseline")
    parser.add_argument("--exclude", action="append", default=[],
                        help="column to ignore; repeat as needed")
    parser.add_argument("--simple", action="append", default=None,
                        help="feature to use as a one-column baseline "
                             "(default: gc_content, seq_length)")
    parser.add_argument("--json", type=Path, default=None, help="write the full report here")
    args = parser.parse_args(argv)

    df = pd.read_csv(args.data, low_memory=False)
    missing = [c for c in [*args.label, args.sequence, args.source]
               if c and c not in df.columns]
    if missing:
        parser.error(f"column(s) not in {args.data}: {missing}")

    reserved = {*args.label, args.sequence, args.source, *args.exclude} - {None}
    features = df.drop(columns=[c for c in reserved if c in df.columns])
    features = features.select_dtypes(include="number")
    if features.empty:
        parser.error("no numeric feature columns left after exclusions")

    report = audit_dataset(
        X=features,
        Y=df[args.label],
        groups=df[args.sequence],
        sources=df[args.source] if args.source else None,
        simple_features=tuple(args.simple) if args.simple else ("gc_content", "seq_length"),
    )

    summary = report.summary()
    print(summary.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    for endpoint in report.endpoints:
        if endpoint.findings:
            print(f"\n{endpoint.endpoint}  [{endpoint.verdict}]")
            for finding in endpoint.findings:
                print(f"  - {finding}")

    if args.json:
        args.json.write_text(json.dumps(report.to_dict(), indent=2) + "\n")
        print(f"\nwrote {args.json}")

    failed = [e.endpoint for e in report.endpoints if e.verdict not in PASSING]
    if failed:
        print(f"\n{len(failed)} of {len(report.endpoints)} endpoint(s) did not pass: "
              f"{', '.join(failed)}")
        return 1
    print(f"\nall {len(report.endpoints)} endpoint(s) passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
