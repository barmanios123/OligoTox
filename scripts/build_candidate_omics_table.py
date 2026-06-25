from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from oligotox.omics.training_join import build_candidate_level_omics_table


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default="data/processed/omics/omics_training_table.csv")
    ap.add_argument("--output", default="data/processed/omics/candidate_level_omics_features.csv")
    args = ap.parse_args()

    df = build_candidate_level_omics_table(
        omics_training_table_path=args.input,
        output_path=args.output,
    )

    print("=" * 80)
    print("CANDIDATE-LEVEL OMICS TABLE BUILT")
    print("=" * 80)
    print(f"Rows: {len(df)}")
    print(f"Cols: {len(df.columns)}")
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()
