from pathlib import Path
import argparse
import pandas as pd

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--group-col", default="sequence_family")
    args = ap.parse_args()

    df = pd.read_csv(args.input)
    bad = (
        df.groupby(args.group_col)["split"]
        .nunique()
        .reset_index(name="n_splits")
        .query("n_splits > 1")
    )

    print(f"rows={len(df)} groups={df[args.group_col].nunique()}")
    print(df["split"].value_counts(dropna=False).to_string())
    print(f"leaky_groups={len(bad)}")
    if len(bad):
        print(bad.head(20).to_string(index=False))
        raise SystemExit(1)

if __name__ == "__main__":
    main()
