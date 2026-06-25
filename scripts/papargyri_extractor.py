from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd
import numpy as np


def pick_column(df, candidates):
    lower_map = {c.lower().strip(): c for c in df.columns}
    for cand in candidates:
        key = cand.lower().strip()
        if key in lower_map:
            return lower_map[key]
    for c in df.columns:
        cl = c.lower().strip()
        for cand in candidates:
            if cand.lower().strip() in cl:
                return c
    return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    path = Path(args.input)
    print(f"Loading: {path}")

    try:
        if path.suffix.lower() == ".xls":
            df = pd.read_excel(path, engine="xlrd")
        else:
            try:
                df = pd.read_excel(path, engine="openpyxl")
            except Exception:
                df = pd.read_excel(path)
    except Exception as e:
        raise SystemExit(f"ERROR loading Excel: {e}")

    print(f"Raw table shape: {df.shape}")
    print(f"Columns: {list(df.columns[:10])} ...")

    gapmer_col = pick_column(df, ["gapmer", "sequence", "oligo", "name"])
    cas_cols = [c for c in df.columns if str(c).lower().startswith("cas_")]

    if not cas_cols:
        raise SystemExit("ERROR: No caspase columns found.")

    work = df.copy()

    for c in cas_cols:
        work[c] = pd.to_numeric(work[c], errors="coerce")

    if "cas_avg" in work.columns:
        work["cas_avg_numeric"] = pd.to_numeric(work["cas_avg"], errors="coerce")
    else:
        work["cas_avg_numeric"] = work[cas_cols].mean(axis=1, skipna=True)

    if gapmer_col is None:
        work["gapmer_id"] = [f"papargyri_gapmer_{i+1}" for i in range(len(work))]
    else:
        work["gapmer_id"] = work[gapmer_col].astype(str).fillna("").replace("nan", "")

    work["candidate_id"] = ["papargyri_" + str(i + 1) for i in range(len(work))]
    work["compound"] = work["gapmer_id"].replace("", np.nan).fillna(work["candidate_id"])
    work["source_type"] = "papargyri_table_s1"
    work["study_id"] = "Papargyri_2020_Table_S1"
    work["species"] = "human_in_vitro"
    work["route"] = "in_vitro"
    work["modality"] = "ASO"
    work["delivery_vehicle"] = "none"
    work["repeat_dose"] = 0

    work["caspase_toxic_flag"] = np.where(work["cas_avg_numeric"] > 300, 1, 0)
    work["caspase_safe_flag"] = np.where(work["cas_avg_numeric"].notna() & (work["cas_avg_numeric"] <= 300), 1, 0)

    work["offtarget_gene_silencing"] = work["caspase_toxic_flag"]
    work["hepatocellular_injury"] = np.where(work["cas_avg_numeric"] > 400, 1, 0)

    endpoint_cols = [
        "hepatocellular_injury",
        "renal_tubular_injury",
        "renal_accumulation",
        "thrombocytopenia",
        "coagulation_effect",
        "complement_activation",
        "innate_immune_activation",
        "cytokine_release",
        "hypersensitivity",
        "cns_neurotoxicity",
        "offtarget_gene_silencing",
        "reproductive_dev_tox",
        "injection_site_reaction",
    ]

    for col in endpoint_cols:
        if col not in work.columns:
            work[col] = 0

    out_cols = [
        "candidate_id",
        "compound",
        "source_type",
        "study_id",
        "species",
        "route",
        "modality",
        "delivery_vehicle",
        "repeat_dose",
        "gapmer_id",
        "cas_avg_numeric",
        "caspase_toxic_flag",
        "caspase_safe_flag",
    ] + endpoint_cols

    out = work[out_cols].copy()
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, index=False)

    toxic_n = int(out["caspase_toxic_flag"].sum())
    safe_n = int(out["caspase_safe_flag"].sum())

    print(f"\n✓  Processed {len(out)} gapmers")
    print(f"   Caspase TOXIC (>300%): {toxic_n} ({toxic_n/len(out)*100:.1f}%)")
    print(f"   Caspase SAFE (<=300%): {safe_n} ({safe_n/len(out)*100:.1f}%)")
    print(f"   Saved to: {out_path}")


if __name__ == "__main__":
    main()
