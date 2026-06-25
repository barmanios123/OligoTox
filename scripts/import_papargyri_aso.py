import pandas as pd
from pathlib import Path

CURATED_PATH = Path("data/raw/curated_oligo_tox_template.csv")
ASO_PATH = Path("data/raw/lna_gapmer/aso_data.csv")  # adjust filename if needed

# Threshold from your labeling scheme: caspase >= 300% = toxic
CASPASE_TOX_THRESHOLD = 300.0

def main():
    if not ASO_PATH.exists():
        raise SystemExit(f"Missing {ASO_PATH} — download Papargyri aso.csv there first.")

    curated = pd.read_csv(CURATED_PATH)
    aso = pd.read_csv(ASO_PATH)

    # TODO: adjust column names to match aso_data.csv
    # Example guesses:
    #  - 'Sequence' for the ASO sequence
    #  - 'CaspasePercent' for caspase activity (% of control)
    seq_col = "Sequence"
    casp_col = "CaspasePercent"

    missing_cols = [c for c in (seq_col, casp_col) if c not in aso.columns]
    if missing_cols:
        raise SystemExit(f"aso_data.csv is missing expected columns: {missing_cols}")

    new_rows = []
    for idx, row in aso.iterrows():
        seq = str(row[seq_col]).strip().upper()
        if not seq:
            continue

        casp = float(row[casp_col])

        # Map caspase to endpoints using your scheme:
        # caspase >= 300% → offtarget_gene_silencing = 1, hepatocellular_injury = 1
        offtarget = int(casp >= CASPASE_TOX_THRESHOLD)
        hepatotox = int(casp >= CASPASE_TOX_THRESHOLD)

        oligo_id = f"OTX-PAP-{idx+1:04d}"

        new_row = {c: "" for c in curated.columns}

        new_row["oligo_id"] = oligo_id
        new_row["sequence"] = seq
        new_row["modality"] = "ASO_gapmer"
        new_row["vehicle"] = "naked"
        new_row["route"] = "subcutaneous"
        new_row["dose_mg_kg"] = ""
        new_row["repeat_dose"] = ""
        new_row["source_paper"] = "papargyri_2019"
        new_row["series_id"] = "papargyri_2019"
        new_row["curation_notes"] = f"Papargyri 2019 gapmer; caspase={casp:.1f}% of control"

        # endpoints
        new_row["offtarget_gene_silencing"] = offtarget
        new_row["hepatocellular_injury"] = hepatotox

        new_rows.append(new_row)

    print(f"Prepared {len(new_rows)} Papargyri rows")

    if not new_rows:
        print("No rows prepared; nothing to do.")
        return

    out_df = pd.concat([curated, pd.DataFrame(new_rows)], ignore_index=True)
    out_df.to_csv(CURATED_PATH, index=False)
    print(f"Wrote {len(out_df)} rows back to {CURATED_PATH}")

if __name__ == "__main__":
    main()
