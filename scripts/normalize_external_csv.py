import argparse
from pathlib import Path
import pandas as pd

CANONICAL_PATH = Path("data/raw/aso_tox_ingest_template.csv")

SYNONYMS = {
    "oligo_id": ["oligo_id", "id", "candidate_id", "compound_id"],
    "sequence": ["sequence", "seq", "oligo_sequence", "Sequence"],
    "target_gene": ["target_gene", "gene", "target", "TargetGene"],
    "source_paper": ["source_paper", "source", "paper", "study"],
    "series_id": ["series_id", "series", "cohort", "dataset"],
    "curation_notes": ["curation_notes", "notes", "comment", "annotation"],
    "modality": ["modality", "type", "oligo_type"],
    "vehicle": ["vehicle", "delivery_vehicle", "formulation"],
    "route": ["route", "administration_route"],
    "dose_mg_kg": ["dose_mg_kg", "dose", "dose_mgkg"],
    "repeat_dose": ["repeat_dose", "repeat", "multi_dose"],
    "platelet_count_baseline": ["platelet_count_baseline", "platelet_baseline"],
    "platelet_count_treated": ["platelet_count_treated", "platelet_treated"],
    "pf4_ps_elisa_od": ["pf4_ps_elisa_od", "pf4_od"],
    "c3a_pg_ml": ["c3a_pg_ml", "c3a"],
    "c5a_pg_ml": ["c5a_pg_ml", "c5a"],
    "sc5b9_ng_ml": ["sc5b9_ng_ml", "sc5b9"],
    "alt_iu_l": ["alt_iu_l", "alt"],
    "ast_iu_l": ["ast_iu_l", "ast"],
    "bun_mg_dl": ["bun_mg_dl", "bun"],
    "creatinine_mg_dl": ["creatinine_mg_dl", "creatinine", "cr"],
    "cytokine_il6": ["cytokine_il6", "il6"],
    "cytokine_tnfa": ["cytokine_tnfa", "tnfa", "tnf_alpha"],
    "cytokine_ifna": ["cytokine_ifna", "ifna", "ifn_alpha"],
    "offtarget_gene_silencing": ["offtarget_gene_silencing", "too_toxic", "offtarget", "caspase_toxic"],
    "hepatocellular_injury": ["hepatocellular_injury", "hepatotox", "liver_toxicity"],
}

def find_source_col(df, canonical):
    for cand in SYNONYMS.get(canonical, []):
        if cand in df.columns:
            return cand
    return None

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Path to external CSV")
    parser.add_argument("--output", default=str(CANONICAL_PATH), help="Path to canonical CSV")
    parser.add_argument("--source-name", default="", help="Default source_paper if missing")
    parser.add_argument("--series-name", default="", help="Default series_id if missing")
    args = parser.parse_args()

    ext = pd.read_csv(args.input)
    canonical = pd.read_csv(args.output)

    out = pd.DataFrame(columns=canonical.columns)

    for col in out.columns:
        src_col = find_source_col(ext, col)
        if src_col:
            out[col] = ext[src_col]
        else:
            out[col] = ""

    if args.source_name:
        out["source_paper"] = out["source_paper"].replace("", args.source_name)
        if "source" in out.columns:
            out["source"] = out["source"].replace("", args.source_name)

    if args.series_name:
        out["series_id"] = out["series_id"].replace("", args.series_name)

    combined = pd.concat([canonical, out], ignore_index=True)
    combined.to_csv(args.output, index=False)

    print(f"Imported {len(out)} rows from {args.input}")
    print(f"Wrote combined canonical CSV to {args.output}")

if __name__ == "__main__":
    main()
