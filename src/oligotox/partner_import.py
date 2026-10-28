
import pandas as pd
import numpy as np
import os, json, hashlib, datetime, warnings
warnings.filterwarnings("ignore")

REQUIRED_COLS  = ["candidate_id", "sequence", "modality"]
VALID_MODALITIES = {"ASO_PS","LNA_gapmer","siRNA","GalNAc_siRNA","PMO",
                    "splice_switching_ASO","mRNA_LNP","saRNA","other"}
ENDPOINT_COLS  = [
    "endpoint_hepatocellular_injury","endpoint_renal_tubular_injury",
    "endpoint_thrombocytopenia","endpoint_complement_activation",
    "endpoint_innate_immune_activation","endpoint_cytokine_release",
    "endpoint_injection_site_reaction","endpoint_offtarget_gene_silencing",
    "endpoint_cns_delayed_toxicity","endpoint_lnp_immune_effects",
    "endpoint_reproductive_developmental_risk","endpoint_genotoxicity",
    "endpoint_ontarget_excessive_knockdown","endpoint_transcriptome_perturbation",
    "endpoint_lncrna_network_disruption"
]


def validate_partner_dataset(df: pd.DataFrame, partner_id: str) -> dict:
    """Validate a partner-uploaded DataFrame against the OligoTox schema.
    Returns a validation report dict with status, errors, warnings, and stats."""
    errors, warnings_list, stats = [], [], {}

    # 1. Required columns
    for col in REQUIRED_COLS:
        if col not in df.columns:
            errors.append(f"MISSING_REQUIRED_COLUMN: {col}")

    if errors:
        return {"status": "FAILED", "partner_id": partner_id,
                "errors": errors, "warnings": warnings_list, "stats": stats}

    # 2. Row count
    stats["n_rows"] = len(df)
    if len(df) < 10:
        warnings_list.append(f"SMALL_DATASET: only {len(df)} rows. Model training requires >= 20 rows with any endpoint labeled.")

    # 3. Sequence validation
    valid_chars = set("ATGCUatgcu")
    bad_seqs = df["sequence"].apply(lambda s: not set(str(s).upper().replace(" ","")) <= valid_chars | {"-","N"} if pd.notna(s) else False)
    stats["sequences_null"]    = int(df["sequence"].isna().sum())
    stats["sequences_invalid"] = int(bad_seqs.sum())
    if stats["sequences_null"] > 0:
        warnings_list.append(f"NULL_SEQUENCES: {stats['sequences_null']} rows missing sequence")
    if stats["sequences_invalid"] > 0:
        warnings_list.append(f"INVALID_SEQUENCES: {stats['sequences_invalid']} rows with non-ATGCU characters")

    # 4. Modality validation
    bad_mod = df[~df["modality"].isin(VALID_MODALITIES)]
    stats["modality_distribution"] = df["modality"].value_counts().to_dict()
    if len(bad_mod) > 0:
        warnings_list.append(f"UNKNOWN_MODALITIES: {set(bad_mod['modality'].unique())} — will be mapped to 'other'")
        df.loc[~df["modality"].isin(VALID_MODALITIES), "modality"] = "other"

    # 5. Endpoint coverage
    ep_found = [ep for ep in ENDPOINT_COLS if ep in df.columns]
    stats["endpoints_present"]   = ep_found
    stats["endpoints_missing"]   = [ep for ep in ENDPOINT_COLS if ep not in df.columns]
    ep_label_counts = {}
    for ep in ep_found:
        n_labeled = int(df[ep].notna().sum())
        n_pos     = int((df[ep] == 1).sum()) if n_labeled > 0 else 0
        ep_label_counts[ep] = {"n_labeled": n_labeled, "n_positive": n_pos,
                               "trainable": n_pos >= 10 and (n_labeled - n_pos) >= 10}
    stats["endpoint_label_counts"] = ep_label_counts
    trainable_eps = [ep for ep, v in ep_label_counts.items() if v["trainable"]]
    stats["trainable_endpoints"] = trainable_eps
    if not trainable_eps:
        warnings_list.append("NO_TRAINABLE_ENDPOINTS: no endpoint has >= 10 positives AND >= 10 negatives")

    # 6. Duplicate check
    n_dup = int(df.duplicated(subset=["candidate_id"]).sum())
    if n_dup > 0:
        warnings_list.append(f"DUPLICATE_CANDIDATE_IDS: {n_dup} duplicate candidate_id values found")

    # 7. Dataset fingerprint (for audit trail, does NOT expose data)
    fingerprint = hashlib.sha256(
        "||".join(sorted(df["sequence"].dropna().astype(str))).encode()
    ).hexdigest()[:16]
    stats["dataset_fingerprint"] = fingerprint
    stats["timestamp"] = datetime.datetime.now(datetime.UTC).isoformat().replace("+00:00", "Z")

    status = "FAILED" if errors else ("WARNING" if warnings_list else "PASSED")
    return {
        "status": status, "partner_id": partner_id,
        "errors": errors, "warnings": warnings_list, "stats": stats
    }


def schema_map(df: pd.DataFrame, column_mapping: dict = None) -> pd.DataFrame:
    """Remap partner column names to OligoTox canonical schema.
    column_mapping example: {"compound_id": "candidate_id", "seq": "sequence"}"""
    if column_mapping:
        df = df.rename(columns=column_mapping)
    # Add missing endpoint cols as NaN
    for ep in ENDPOINT_COLS:
        if ep not in df.columns:
            df[ep] = np.nan
    return df


def register_partner_import(validation_report: dict, run_registry_path: str,
                             partner_id: str, file_path: str = None) -> str:
    """Log a validated partner import to the run registry. Returns the run_id."""
    import uuid
    run_id = f"partner_import_{partner_id}_{datetime.datetime.now(datetime.UTC).strftime('%Y%m%d_%H%M%S')}_{uuid.uuid4().hex[:8]}"
    record = {
        "run_id":       run_id,
        "job_name":     "partner_dataset_import",
        "partner_id":   partner_id,
        "status":       "succeeded" if validation_report["status"] != "FAILED" else "failed",
        "trigger":      "partner_upload",
        "started_at":   validation_report["stats"].get("timestamp"),
        "finished_at":  validation_report["stats"].get("timestamp"),
        "file_path":    file_path,
        "parameters":   {"n_rows": validation_report["stats"].get("n_rows"),
                         "trainable_endpoints": validation_report["stats"].get("trainable_endpoints"),
                         "fingerprint": validation_report["stats"].get("dataset_fingerprint")},
        "artifacts":    [],
        "error":        "; ".join(validation_report["errors"]) if validation_report["errors"] else None
    }
    os.makedirs(run_registry_path, exist_ok=True)
    record_path = os.path.join(run_registry_path, f"{run_id}.json")
    with open(record_path, "w") as f: json.dump(record, f, indent=2)
    return run_id
