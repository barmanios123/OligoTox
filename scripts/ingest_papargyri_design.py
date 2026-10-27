"""Recover the gapmer design annotation the pipeline never ingested.

data/external/papargyri_2020_LNA_512gapmers_TableS1.xlsx is a legacy .xls
despite the extension, which is why nothing had opened it. It holds 772 rows
and 165 columns for the same gapmers that supply 772 of the 1,614 training
rows, including design information the training pipeline never used:

  lnadna        per-position chemistry, L for LNA and D for DNA, so the real
                gapmer architecture rather than a guessed lna_fraction
  sequence      case-encoded, uppercase LNA and lowercase DNA
  tm_*          measured melting temperatures
  gapmer.length, number.lna, number.dna, gap, length.gap, flank5, flank3

Everywhere else in this repo, chemistry is either absent (the external CNS
and potency sets) or was inconsistently filled upstream and recomputed from
sequence. This is the one place with the real design, and it is design-time
information, available before an oligo is made.

Deliberately NOT ingested: cas_avg and kd_avg, the caspase and knockdown
readouts. They are measured after dosing, they were removed from training as
leakage, and they must stay out.

Writes data/processed/papargyri_design_features.csv and a coverage report.
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
WORKBOOK = ROOT / "data" / "external" / "papargyri_2020_LNA_512gapmers_TableS1.xlsx"
PROCESSED = ROOT / "data" / "processed"
RESULTS = ROOT / "reports" / "results"

# Measured after dosing; these are the leakage columns and never leave here.
FORBIDDEN = {"cas_1", "cas_2", "cas_3", "cas_avg", "kd_avg"}

DESIGN_COLUMNS = [
    "gapmer", "sequence", "lnadna", "gapmer.length", "number.lna", "number.dna",
    "gap", "length.gap", "flank5", "flank3", "tm_up", "tm_down", "tm_dt", "tm_avg",
]


def design_from_pattern(pattern: str) -> dict:
    """Gapmer architecture from the per-position L/D string."""
    p = str(pattern).strip().upper()
    if not p or set(p) - {"L", "D"}:
        return {}
    n = len(p)
    first_d = p.find("D")
    last_d = p.rfind("D")
    return {
        "design_length": float(n),
        "design_lna_count": float(p.count("L")),
        "design_lna_fraction": p.count("L") / n,
        "design_flank5_len": float(first_d) if first_d >= 0 else float(n),
        "design_flank3_len": float(n - last_d - 1) if last_d >= 0 else 0.0,
        "design_gap_len": float(last_d - first_d + 1) if first_d >= 0 else 0.0,
        # A classic gapmer is LNA flanks around a contiguous DNA gap.
        "design_is_gapmer": float(first_d > 0 and last_d < n - 1 and "L" not in p[first_d:last_d + 1]),
    }


def main() -> None:
    if not WORKBOOK.exists():
        raise SystemExit(f"missing {WORKBOOK}")
    raw = pd.read_excel(WORKBOOK, sheet_name="data", engine="xlrd")

    leaked = FORBIDDEN & set(raw.columns)
    keep = [c for c in DESIGN_COLUMNS if c in raw.columns]
    df = raw[keep].copy()

    df["sequence_norm"] = (
        df["sequence"].astype(str).str.upper().str.replace("U", "T", regex=False).str.strip()
    )
    design = pd.DataFrame([design_from_pattern(p) for p in df.get("lnadna", [])], index=df.index)
    out = pd.concat([df, design], axis=1)

    assert not (FORBIDDEN & set(out.columns)), "a post-dose readout reached the output"

    PROCESSED.mkdir(parents=True, exist_ok=True)
    out.to_csv(PROCESSED / "papargyri_design_features.csv", index=False)

    # How much of the training set does this cover?
    training = pd.read_csv(PROCESSED / "training_dataset.csv", low_memory=False)
    train_seqs = set(
        training["sequence"].astype(str).str.upper().str.replace("U", "T", regex=False)
    )
    matched = sorted(train_seqs & set(out["sequence_norm"]))

    report = {
        "workbook": str(WORKBOOK.relative_to(ROOT)),
        "format_note": "legacy .xls despite the .xlsx extension, which is why it was never opened",
        "rows": int(len(out)),
        "design_columns_recovered": keep + list(design.columns),
        "post_dose_columns_excluded": sorted(leaked),
        "training_sequences_total": len(train_seqs),
        "training_sequences_matched": len(matched),
        "coverage": len(matched) / len(train_seqs) if train_seqs else 0.0,
        "usable_as_features": (
            "design_* columns are design-time and legitimate. They give the real gapmer "
            "architecture for the matched sequences, where the pipeline currently infers "
            "lna_fraction from the sequence alone."
        ),
    }
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / "papargyri_design_coverage.json").write_text(json.dumps(report, indent=2) + "\n")

    print(f"rows {len(out)}, design columns {len(keep) + len(design.columns)}")
    print(f"excluded post-dose readouts: {sorted(leaked)}")
    print(f"covers {len(matched)} of {len(train_seqs)} training sequences "
          f"({report['coverage']:.1%})")
    if "design_is_gapmer" in out:
        print(f"classic gapmer architecture: {int(out['design_is_gapmer'].sum())} of {len(out)}")
        print(f"median LNA fraction {out['design_lna_fraction'].median():.2f}, "
              f"median gap {out['design_gap_len'].median():.0f} nt")


if __name__ == "__main__":
    main()
