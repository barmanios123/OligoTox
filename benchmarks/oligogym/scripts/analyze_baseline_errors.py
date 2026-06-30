from pathlib import Path
import argparse
import json
import pandas as pd

def add_basic_features(df):
    seq = df["sequence"].astype(str).str.upper()
    df = df.copy()
    df["seq_len"] = seq.str.len()
    df["gc_frac"] = seq.str.count("G").add(seq.str.count("C")) / df["seq_len"].clip(lower=1)
    df["u_count"] = seq.str.count("U")
    df["t_count"] = seq.str.count("T")
    df["cpg_count"] = seq.str.count("CG")
    return df

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--predictions", required=True)
    ap.add_argument("--output-prefix", required=True)
    args = ap.parse_args()

    df = pd.read_csv(args.predictions)
    required = {"sequence", "label", "pred_label", "pred_prob"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    df = add_basic_features(df)
    df["error_type"] = "correct"
    df.loc[(df["label"] == 1) & (df["pred_label"] == 0), "error_type"] = "false_negative"
    df.loc[(df["label"] == 0) & (df["pred_label"] == 1), "error_type"] = "false_positive"
    df["abs_error"] = (df["label"] - df["pred_prob"]).abs()

    summary = (
        df.groupby("error_type")
        .agg(
            n=("sequence", "size"),
            mean_prob=("pred_prob", "mean"),
            mean_len=("seq_len", "mean"),
            mean_gc=("gc_frac", "mean"),
            mean_cpg=("cpg_count", "mean"),
        )
        .reset_index()
    )

    hardest = df.sort_values("abs_error", ascending=False).head(50)
    fps = df[df["error_type"] == "false_positive"].sort_values("pred_prob", ascending=False).head(50)
    fns = df[df["error_type"] == "false_negative"].sort_values("pred_prob", ascending=True).head(50)

    out_prefix = Path(args.output_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    summary_csv = out_prefix.with_suffix(".error_summary.csv")
    hardest_csv = out_prefix.with_suffix(".hardest_cases.csv")
    fp_csv = out_prefix.with_suffix(".false_positives.csv")
    fn_csv = out_prefix.with_suffix(".false_negatives.csv")
    json_path = out_prefix.with_suffix(".error_summary.json")

    summary.to_csv(summary_csv, index=False)
    hardest.to_csv(hardest_csv, index=False)
    fps.to_csv(fp_csv, index=False)
    fns.to_csv(fn_csv, index=False)

    with open(json_path, "w") as f:
        json.dump(summary.to_dict(orient="records"), f, indent=2)

    print("Wrote:")
    print(" -", summary_csv)
    print(" -", hardest_csv)
    print(" -", fp_csv)
    print(" -", fn_csv)
    print(" -", json_path)
    print(summary.to_string(index=False))

if __name__ == "__main__":
    main()
