from pathlib import Path
import argparse
import pandas as pd
import numpy as np

def make_group_splits(df, group_col, seed=42, train_frac=0.7, val_frac=0.15):
    rng = np.random.default_rng(seed)
    groups = df[group_col].fillna("UNGROUPED").astype(str).unique().tolist()
    rng.shuffle(groups)

    n = len(groups)
    n_train = max(1, int(round(n * train_frac)))
    n_val = max(1, int(round(n * val_frac))) if n >= 3 else 0

    train_groups = set(groups[:n_train])
    val_groups = set(groups[n_train:n_train + n_val])
    test_groups = set(groups[n_train + n_val:])

    if not test_groups and val_groups:
        moved = next(iter(val_groups))
        val_groups.remove(moved)
        test_groups.add(moved)
    if not test_groups and train_groups:
        moved = next(iter(train_groups))
        train_groups.remove(moved)
        test_groups.add(moved)

    split = []
    for g in df[group_col].fillna("UNGROUPED").astype(str):
        if g in train_groups:
            split.append("train")
        elif g in val_groups:
            split.append("val")
        else:
            split.append("test")
    return split

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--group-col", default="sequence_family")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    df = pd.read_csv(args.input)
    if args.group_col not in df.columns:
        raise ValueError(f"group column '{args.group_col}' not found in {args.input}")

    df["split"] = make_group_splits(df, args.group_col, seed=args.seed)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)

    print(f"wrote {out}")
    print(df["split"].value_counts(dropna=False).to_string())
    overlap = (
        df.groupby(args.group_col)["split"]
        .nunique()
        .reset_index(name="n_splits")
        .query("n_splits > 1")
    )
    print(f"group leakage count: {len(overlap)}")

if __name__ == "__main__":
    main()
