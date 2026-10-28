"""Convert the ASO Atlas pickle to parquet without executing its code.

ASO Atlas is 188,521 RNase H gapmer ASOs with qRT-PCR knockdown, extracted
from 417 USPTO patents, Apache-2.0 licensed:

    https://github.com/barneyhill/aso_atlas
    Accurately modeling RNase H-mediated antisense oligonucleotide efficacy.
    Mol Ther Nucleic Acids 2026, doi:10.1016/j.omtn.2026.103004

It ships as `data/aso_atlas.pkl`, and that is the problem this script exists
for. **Unpickling runs code.** An opcode scan of the file shows it
reconstructs two classes from the authors' own package:

    src.formats.Chemistry
    modifications.Modification

Nothing in the scan looks hostile, and no `os`, `subprocess` or `eval`
appears. But loading it normally means importing their package, and `src` is
generic enough to collide with this repository's own `src/` directory. A
dataset should not require running anyone's code to read.

So the loader below refuses to import anything outside numpy, pandas and the
standard library, and replaces every foreign class with an inert stub that
records its state and does nothing else. The DataFrame comes through; the
behaviour does not.

    # once, with the file downloaded to its own directory
    python scripts/ingest_aso_atlas.py path/to/aso_atlas.pkl

Writes data/external/aso_atlas.parquet and a summary. After that nothing in
the pipeline touches the pickle again.

Treat the downloaded file as untrusted: put it in its own empty directory and
pass the path as an argument, which is why this script takes one rather than
reaching for a fixed location.
"""
from __future__ import annotations

import argparse
import json
import pickle
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "data" / "external"
RESULTS = ROOT / "reports" / "results"

# Everything else becomes a stub.
ALLOWED_ROOTS = {"numpy", "pandas", "builtins", "collections", "_codecs", "copyreg"}


class _InertStub:
    """Stands in for a foreign class. Records state, executes nothing."""

    __slots__ = ("_args", "_state")

    def __init__(self, *args, **_kwargs):
        self._args = args
        self._state = None

    def __setstate__(self, state):
        self._state = state

    def describe(self):
        state = self._state if self._state is not None else self._args
        if isinstance(state, dict):
            return ";".join(f"{k}={v}" for k, v in sorted(state.items()))
        return str(state)

    def __repr__(self):
        return f"<stub {self.describe()[:60]}>"


class RestrictedUnpickler(pickle.Unpickler):
    """Allows data types through; refuses to import foreign code."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.stubbed: set[str] = set()

    def find_class(self, module: str, name: str):
        if module.split(".")[0] in ALLOWED_ROOTS:
            return super().find_class(module, name)
        self.stubbed.add(f"{module}.{name}")
        return type(f"Stub_{name}", (_InertStub,), {})


def load(path: Path):
    with path.open("rb") as handle:
        unpickler = RestrictedUnpickler(handle)
        obj = unpickler.load()
    return obj, sorted(unpickler.stubbed)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("pickle_path", type=Path,
                        help="downloaded aso_atlas.pkl, in its own directory")
    parser.add_argument("--out", type=Path, default=OUT / "aso_atlas.parquet")
    args = parser.parse_args(argv)

    if not args.pickle_path.exists():
        raise SystemExit(f"not found: {args.pickle_path}")

    import pandas as pd

    frame, stubbed = load(args.pickle_path)
    if not isinstance(frame, pd.DataFrame):
        raise SystemExit(f"expected a DataFrame, got {type(frame).__name__}")

    print(f"loaded {len(frame):,} rows x {len(frame.columns)} columns")
    print(f"foreign classes stubbed rather than imported: {stubbed or 'none'}")

    # Stubs cannot be written to parquet, so flatten them to their state.
    flattened = frame.copy()
    for column in flattened.columns:
        sample = flattened[column].dropna().head(1)
        if len(sample) and isinstance(sample.iloc[0], _InertStub):
            flattened[column] = flattened[column].map(
                lambda v: v.describe() if isinstance(v, _InertStub) else v
            )
            print(f"  flattened stub column: {column}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    flattened.to_parquet(args.out, index=False)

    summary = {
        "source": "https://github.com/barneyhill/aso_atlas",
        "paper_doi": "10.1016/j.omtn.2026.103004",
        "licence": "Apache-2.0",
        "rows": int(len(flattened)),
        "columns": [str(c) for c in flattened.columns],
        "foreign_classes_stubbed": stubbed,
        "note": "Converted with a restricted unpickler; no code from the archive "
                "was executed. The parquet is what the pipeline reads.",
    }
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / "aso_atlas_ingest.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"\nwrote {args.out.relative_to(ROOT)} and the ingest summary")
    print("Declare the readout direction in data/registry/sources.json before training on it.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
