"""Count off-target complementarity against the whole human transcriptome.

The earlier test used the cached target transcripts as a background: 26
transcripts, 0.12 Mb against a transcriptome of roughly 1,280 Mb. At that
size the counts are noise, so the null it produced said nothing about the
mechanism.

This streams the Ensembl cDNA FASTA once, matching every oligo's probes in a
single pass with an Aho-Corasick automaton. Searching 200-odd patterns
separately over 1.28 Gb would take hours; one automaton pass takes minutes.

Three probes per oligo, each a reverse complement of a different part:

    full      the whole oligo, a perfect off-target site
    gap       the central DNA gap, which is what RNase H1 needs to cleave
    seed      positions 2-8, the miRNA-like seed

Patterns can straddle the boundary between streamed chunks, so a tail of
(longest pattern - 1) characters is carried across. Without it, matches at
chunk edges are silently lost.

    python scripts/offtarget_transcriptome.py path/to/cdna.fa.gz

Writes data/processed/offtarget_transcriptome.csv, one row per unique
training sequence.
"""
from __future__ import annotations

import argparse
import gzip
import json
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from oligotox.features.offtarget import DEFAULT_GAP, SEED, reverse_complement  # noqa: E402

PROCESSED = ROOT / "data" / "processed"
RESULTS = ROOT / "reports" / "results"
CHUNK = 4_000_000  # characters of sequence per pass


def probes_for(sequence: str) -> dict[str, str]:
    """The three reverse-complement probes, skipping any too short to mean much."""
    oligo = sequence.strip().upper().replace("U", "T")
    out = {"full": reverse_complement(oligo)}
    gap_start, gap_end = DEFAULT_GAP
    if len(oligo) > gap_end:
        out["gap"] = reverse_complement(oligo[gap_start:gap_end])
    seed = oligo[SEED[0]:SEED[1]]
    if len(seed) >= 6:
        out["seed"] = reverse_complement(seed)
    return out


def build_automaton(sequences: list[str]):
    import ahocorasick

    automaton = ahocorasick.Automaton()
    pattern_owners: dict[str, list[tuple[str, str]]] = defaultdict(list)
    for sequence in sequences:
        for kind, probe in probes_for(sequence).items():
            pattern_owners[probe].append((sequence, kind))

    for probe, owners in pattern_owners.items():
        automaton.add_word(probe, (probe, owners))
    automaton.make_automaton()
    longest = max(len(p) for p in pattern_owners)
    return automaton, longest, len(pattern_owners)


def stream_counts(fasta: Path, automaton, longest: int):
    """One pass over the FASTA, counting every pattern occurrence."""
    counts: dict[str, int] = defaultdict(int)
    total_bases = 0
    buffer = ""
    opener = gzip.open if fasta.suffix == ".gz" else open

    def scan(text: str):
        for _, (probe, _owners) in automaton.iter(text):
            counts[probe] += 1

    with opener(fasta, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                # A pattern must not match across two transcripts.
                scan(buffer)
                total_bases += len(buffer)
                buffer = ""
                continue
            buffer += line.strip().upper()
            if len(buffer) >= CHUNK:
                scan(buffer)
                total_bases += len(buffer) - (longest - 1)
                buffer = buffer[-(longest - 1):]
    scan(buffer)
    total_bases += len(buffer)
    return counts, total_bases


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("fasta", type=Path, help="Ensembl cDNA FASTA, .fa or .fa.gz")
    args = parser.parse_args(argv)
    if not args.fasta.exists():
        raise SystemExit(f"not found: {args.fasta}")

    import pandas as pd

    raw = pd.read_csv(PROCESSED / "training_dataset.csv", low_memory=False)
    sequences = sorted({
        s.upper().replace("U", "T")
        for s in raw["sequence"].dropna().astype(str)
        if s and set(s.upper().replace("U", "T")) <= set("ACGT")
    })
    print(f"{len(sequences)} unique training sequences")

    automaton, longest, n_patterns = build_automaton(sequences)
    print(f"{n_patterns} distinct probes, longest {longest} nt")
    print("streaming the transcriptome (one pass)...")

    counts, total_bases = stream_counts(args.fasta, automaton, longest)
    megabases = total_bases / 1_000_000
    print(f"scanned {megabases:,.1f} Mb\n")

    rows = []
    for sequence in sequences:
        probes = probes_for(sequence)
        row = {"sequence": sequence}
        for kind in ("full", "gap", "seed"):
            probe = probes.get(kind)
            hits = counts.get(probe, 0) if probe else 0
            row[f"offtarget_{kind}_hits"] = hits
            row[f"offtarget_{kind}_per_mb"] = round(hits / megabases, 6)
        rows.append(row)

    frame = pd.DataFrame(rows)
    PROCESSED.mkdir(parents=True, exist_ok=True)
    frame.to_csv(PROCESSED / "offtarget_transcriptome.csv", index=False)

    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / "offtarget_transcriptome.json").write_text(json.dumps({
        "source": str(args.fasta.name),
        "megabases_scanned": round(megabases, 1),
        "sequences": len(sequences),
        "probes": n_patterns,
        "note": "Whole-transcriptome background. Counts are absolute occurrences of "
                "each reverse-complement probe across Ensembl cDNA.",
        "summary": {
            kind: {
                "min": int(frame[f"offtarget_{kind}_hits"].min()),
                "median": float(frame[f"offtarget_{kind}_hits"].median()),
                "max": int(frame[f"offtarget_{kind}_hits"].max()),
            }
            for kind in ("full", "gap", "seed")
        },
    }, indent=2) + "\n")

    for kind in ("full", "gap", "seed"):
        column = frame[f"offtarget_{kind}_hits"]
        print(f"  {kind:5s} hits: min {column.min():6d}  median {column.median():9.0f}  "
              f"max {column.max():8d}")
    print(f"\nwrote {(PROCESSED / 'offtarget_transcriptome.csv').relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
