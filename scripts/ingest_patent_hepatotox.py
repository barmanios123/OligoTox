"""Add the patent liver readouts to the training corpus, with the label rule
written down.

`hepatocellular_injury` is trained on 40 unique sequences, 19 of them
positive, and 579 of its 607 rows come from one source. That is the source
confounding its own evidence entry describes, and no amount of modelling
fixes it. This adds 444 study-by-sequence observations over 256 sequences
from 30 US patents, extracted by scripts/extract_patent_hepatotox.py.

Turning ALT into a toxic/not label is a modelling decision, so it is here
rather than implied:

    label  = 1 when the compound's ALT reaches ALT_FOLD x its own study's
             control, taken at the worst dose that study tested
    group  = one row per (publication, table, sequence), so two studies of
             one compound stay two observations and the sequence-grouped CV
             keeps them on the same side of a fold

Fold change against the study's own control, not an absolute ALT, because
control ALT across these patents runs from 9 to 57 U/l and an absolute cut
would label whole studies by their baseline. 3x is the conventional cut for
a hepatocellular signal.

**What is assumed and what is not.** ps_fraction is set to 1.0: these are
Ionis ASO gapmers and a phosphorothioate backbone is uniform across the
class. lna_fraction is left empty, because the patent gives a gapmer motif
such as 5-10-5 but not always what the wings are made of, and the trainer
imputes an empty feature to 0.0. That means every patent row reads as
"no LNA", which is a source indicator rather than a measurement, and it is
the first thing to suspect if this source appears to help. The run prints
the source-only baseline for exactly that reason.

    python scripts/ingest_patent_hepatotox.py --dry-run
    python scripts/ingest_patent_hepatotox.py

Writes data/processed/training_dataset.csv and a manifest beside it.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
PATENTS = ROOT / "data" / "external" / "patent_hepatotox.csv"
TRAINING = ROOT / "data" / "processed" / "training_dataset.csv"
MANIFEST = ROOT / "reports" / "results" / "patent_ingest.json"

SOURCE = "ionis_patents"
LABEL = "endpoint_hepatocellular_injury"
ALT_FOLD = 3.0

# A cEt or LNA wing is a different molecule from an MOE wing, and the
# patents that state it say so in their own words.
LNA_LIKE = ("cet", "lna", "bna")


def _modality(chemistry: str | float, motif: str | float) -> str:
    """Everything here is a gapmer; the wings are what differ."""
    text = str(chemistry or "").lower()
    if any(token in text for token in LNA_LIKE):
        return "LNA_gapmer"
    return "ASO_gapmer"


def build(patents: pd.DataFrame) -> pd.DataFrame:
    """One labelled row per study per sequence."""
    keys = ["publication", "table", "sequence"]
    worst = patents.loc[patents.groupby(keys)["alt_fold_vs_control"].idxmax()]

    out = pd.DataFrame({
        "candidate_id": (worst["publication"] + "_t" + worst["table"].astype(str)
                         + "_" + worst["isis_number"].astype(str)),
        "sequence": worst["sequence"],
        "modality": [_modality(c, m) for c, m
                     in zip(worst["chemistry"], worst["motif"])],
        # Declared, not measured: see the module docstring.
        "ps_fraction": 1.0,
        "lna_fraction": pd.NA,
        "source": SOURCE,
        "seq_length": worst["sequence"].str.len().astype(float),
        "gc_content": worst["sequence"].apply(
            lambda s: (s.count("G") + s.count("C")) / len(s)),
        LABEL: (worst["alt_fold_vs_control"] >= ALT_FOLD).astype(float),
        "notes": (worst["publication"] + " table " + worst["table"].astype(str)
                  + ", ALT " + worst["alt"].astype(str) + " vs control "
                  + worst["control_alt"].astype(str)
                  + " (" + worst["alt_fold_vs_control"].round(2).astype(str) + "x)"),
    })
    return out.reset_index(drop=True)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--dry-run", action="store_true",
                        help="report what would change and write nothing")
    parser.add_argument("--remove", action="store_true",
                        help="take the patent rows back out again")
    args = parser.parse_args(argv)

    if not PATENTS.exists():
        print(f"no extraction at {PATENTS.relative_to(ROOT)}; run "
              "scripts/extract_patent_hepatotox.py first")
        return 1

    training = pd.read_csv(TRAINING, low_memory=False)
    existing = training[training["source"] != SOURCE]
    print(f"corpus without patents: {len(existing)} row(s)")

    if args.remove:
        existing.to_csv(TRAINING, index=False)
        print(f"removed {len(training) - len(existing)} patent row(s)")
        return 0

    rows = build(pd.read_csv(PATENTS))
    labelled = rows[LABEL].notna()
    before = existing[existing[LABEL].notna()]
    merged = pd.concat([existing, rows], ignore_index=True)
    after = merged[merged[LABEL].notna()]

    summary = {
        "source": SOURCE,
        "label_rule": f"ALT >= {ALT_FOLD}x the study's own control, at the "
                      "worst dose that study tested",
        "grouping": "one row per (publication, table, sequence)",
        "rows_added": int(len(rows)),
        "positive_rows_added": int(rows[LABEL].sum()),
        "sequences": {
            "before": int(before["sequence"].nunique()),
            "after": int(after["sequence"].nunique()),
        },
        "positive_sequences": {
            "before": int(before[before[LABEL] == 1]["sequence"].nunique()),
            "after": int(after[after[LABEL] == 1]["sequence"].nunique()),
        },
        "base_rate": {
            "before": round(float(before[LABEL].mean()), 4),
            "after": round(float(after[LABEL].mean()), 4),
        },
        "modality": rows["modality"].value_counts().to_dict(),
        "assumed": {
            "ps_fraction": "set to 1.0 for every patent row",
            "lna_fraction": "left empty; the trainer imputes 0.0, so it reads "
                            "as a source indicator and the source-only "
                            "baseline is the check on it",
        },
        "warning": "hepatocellular_injury already fails to transfer to the "
                   "NAR 2022 screen. Whether this source helps is settled by "
                   "that check and by the source-only baseline, not by the "
                   "number of rows added.",
    }
    print(json.dumps(summary, indent=2))

    if args.dry_run:
        print("\ndry run, nothing written")
        return 0

    assert labelled.all(), "every patent row must carry a label"
    merged.to_csv(TRAINING, index=False)
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    MANIFEST.write_text(json.dumps(summary, indent=2) + "\n")
    print(f"\nwrote {TRAINING.relative_to(ROOT)} and "
          f"{MANIFEST.relative_to(ROOT)}")
    print("now re-run: python scripts/run_pipeline.py --trigger patent-ingest")
    return 0


if __name__ == "__main__":
    sys.exit(main())
