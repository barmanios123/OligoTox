#!/usr/bin/env python3
"""
OligoTox Operator Console & Sequence Analysis Platform
Production-ready Streamlit app wired to the run_oligotox inference backend.

This UI distinguishes:
1. Trained backend outputs from run_oligotox
2. Optional PK/PD bridge outputs
3. Mechanistic proxy overlays used for interpretation/triage
"""

import json
import math
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Tuple

import pandas as pd
import streamlit as st

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / "src"))

from oligotox.config.paths import (  # noqa: E402
    PROCESSED_DATA_DIR,
    MODEL_ARTIFACTS_DIR,
    METRICS_ARTIFACTS_DIR,
)
from oligotox.inference.api import run_oligotox  # noqa: E402
from oligotox.ui_state import log_run, get_runs, get_last_run  # noqa: E402

try:
    from oligotox.pkpd.bridge import compute_pkpd_enriched_features  # noqa: E402

    PKPD_AVAILABLE = True
except Exception:
    PKPD_AVAILABLE = False

    def compute_pkpd_enriched_features(*args, **kwargs):
        return {}


APP_TITLE = "🧬 OligoTox: Oligonucleotide Toxicity Prediction Platform"
APP_CAPTION = (
    "Mechanism-aware toxicity prediction and safety decision support for ASOs, "
    "siRNAs, and RNA-targeting therapeutics"
)

TAB_NAMES = [
    "Sequence Analysis",
    "Overview",
    "Curation",
    "Runs",
    "Models",
    "Predictions",
    "Evaluation",
    "Comparator",
    "Assays",
    "Regulatory Report",
]


def normalize_sequence(seq: str) -> str:
    return seq.upper().replace(" ", "").replace("\n", "").replace("\t", "").replace("T", "U")


def validate_sequence(seq: str) -> Tuple[bool, str]:
    if not seq:
        return False, "Sequence is empty."
    allowed = set("ACGUN")
    bad = sorted(set(ch for ch in seq if ch not in allowed))
    if bad:
        return False, f"Invalid characters found: {', '.join(bad)}. Allowed: A, C, G, U, N, and DNA T (auto-normalized to U)."
    if len(seq) < 8:
        return False, "Sequence is too short for meaningful interpretation."
    if len(seq) > 80:
        return False, "Sequence is too long for this interface. Please use <= 80 nt."
    return True, ""


def confidence_label(prob: float) -> str:
    dist = abs(float(prob) - 0.5)
    if dist >= 0.40:
        return "High"
    if dist >= 0.25:
        return "Moderate"
    return "Low"


def uncertainty_band(prob: float) -> float:
    return round(1.0 - abs(float(prob) - 0.5) * 2.0, 4)


def safe_float(x: Any, default: float = 0.0) -> float:
    try:
        if x is None or (isinstance(x, float) and math.isnan(x)):
            return default
        return float(x)
    except Exception:
        return default


def build_pkpd_metrics(
    sequence: str,
    modality: str,
    chemistry: str,
    dose_mg_kg: float | None,
    schedule: str | None,
) -> Dict[str, float]:
    if not PKPD_AVAILABLE:
        return {}
    if dose_mg_kg is None or dose_mg_kg <= 0:
        return {}
    if not schedule:
        return {}
    try:
        return compute_pkpd_enriched_features(
            sequence=sequence,
            modality=modality,
            chemistry=chemistry,
            dose_mg_kg=float(dose_mg_kg),
            schedule=schedule,
        )
    except Exception:
        return {}


def detect_modality_class(modality: str, chemistry: str, conjugate: str, sequence: str) -> str:
    m = (modality or "").lower()
    c = (chemistry or "").lower()
    g = (conjugate or "").lower()
    n = len(sequence or "")

    if "sirna" in m and "galnac" in g:
        return "GalNAc-siRNA"
    if "sirna" in m:
        return "siRNA"
    if "aso" in m and ("lna" in c or "moe" in c):
        return "Gapmer/ASO"
    if "aso" in m:
        return "ASO"
    if "lncrna" in m:
        return "lncRNA-targeting oligo"
    if n in (19, 20, 21, 22) and "galnac" in g:
        return "GalNAc-siRNA"
    if n in (19, 20, 21, 22):
        return "siRNA-like"
    return "Other ONT"


def engine_dataset_registry(modality_class: str) -> Dict[str, Any]:
    base = {
        "tox_engine": "OligoTox endpoint classifiers (backend via run_oligotox)",
        "mechanistic_layer": "D1/D2/D3 attribution + interpretable rule/proxy layer",
        "benchmark_validation": ["OligoGym", "openASO", "PFRED", "scASO"],
        "regulatory_basis": ["OSWG off-target workflow", "FDA ONT nonclinical safety framing"],
        "backend_status": "Trained OligoTox toxicity backend wired; external efficacy engine family names may reflect architecture/roadmap unless separately integrated.",
    }

    if modality_class == "GalNAc-siRNA":
        base.update(
            {
                "efficacy_engine": "GalNAc-siRNA scoring context",
                "pkpd_bridge": "GalNAc PK/PD translational bridge",
                "primary_benchmarks": ["OligoGym", "Huesken siRNA", "PBMC immunostimulation"],
            }
        )
    elif modality_class in ("siRNA", "siRNA-like"):
        base.update(
            {
                "efficacy_engine": "siRNA-oriented ONT scoring context",
                "pkpd_bridge": "Optional delivery-aware bridge",
                "primary_benchmarks": ["OligoGym", "Huesken siRNA", "PBMC immunostimulation"],
            }
        )
    elif modality_class == "Gapmer/ASO":
        base.update(
            {
                "efficacy_engine": "ASO/gapmer-oriented ONT scoring context",
                "pkpd_bridge": "ASO liver/kidney exposure bridge",
                "primary_benchmarks": ["openASO", "PFRED", "LNA gapmer cytotox", "scASO"],
            }
        )
    else:
        base.update(
            {
                "efficacy_engine": "General ONT scoring stack",
                "pkpd_bridge": "Optional exposure bridge",
                "primary_benchmarks": ["OligoGym", "openASO", "PFRED"],
            }
        )
    return base


def compute_intermediate_safety_metrics(
    sequence: str,
    chemistry: str,
    conjugate: str,
    modality: str,
) -> Dict[str, float]:
    seq = (sequence or "").upper().replace("T", "U")
    n = len(seq) if seq else 0
    gc = ((seq.count("G") + seq.count("C")) / n) if n else 0.0
    u_frac = (seq.count("U") / n) if n else 0.0
    cpg_count = seq.count("CG")
    gu_count = seq.count("GU")
    au_count = seq.count("AU")

    chemistry_upper = (chemistry or "").upper()
    conjugate_upper = (conjugate or "").upper()

    ps_fraction = 1.0 if chemistry_upper == "PS" else 0.5 if "MIX" in chemistry_upper else 0.0
    lna_flag = 1.0 if "LNA" in chemistry_upper else 0.0
    moe_flag = 1.0 if "MOE" in chemistry_upper else 0.0
    ome_flag = 1.0 if "OMe".upper() in chemistry_upper else 0.0
    galnac_flag = 1.0 if "GALNAC" in conjugate_upper else 0.0

    hybridization_risk = min(1.0, round(0.35 * gc + 0.12 * cpg_count + 0.08 * gu_count + 0.05 * au_count, 4))
    tlr78_risk = min(1.0, round(0.55 * u_frac + 0.12 * (n >= 19) + 0.10 * cpg_count + 0.10 * ps_fraction - 0.12 * ome_flag, 4))
    platelet_ps_risk = min(1.0, round(0.55 * ps_fraction + 0.10 * (n >= 16) + 0.10 * tlr78_risk, 4))
    hepatic_gapmer_risk = min(1.0, round(0.30 * lna_flag + 0.18 * gc + 0.18 * hybridization_risk + 0.08 * ps_fraction, 4))
    renal_accumulation_proxy = min(1.0, round(0.42 * ps_fraction + 0.18 * (n >= 18) + 0.10 * moe_flag + 0.05 * galnac_flag, 4))

    return {
        "length_nt": n,
        "gc_fraction": round(gc, 4),
        "u_fraction": round(u_frac, 4),
        "cpg_count": cpg_count,
        "gu_count": gu_count,
        "au_count": au_count,
        "ps_fraction_proxy": round(ps_fraction, 4),
        "lna_flag": lna_flag,
        "moe_flag": moe_flag,
        "ome_flag": ome_flag,
        "galnac_flag": galnac_flag,
        "hybridization_risk_proxy": hybridization_risk,
        "tlr78_risk_proxy": tlr78_risk,
        "platelet_ps_risk_proxy": platelet_ps_risk,
        "hepatic_gapmer_risk_proxy": hepatic_gapmer_risk,
        "renal_accumulation_proxy": renal_accumulation_proxy,
    }


def build_assay_recommendations(result: Dict[str, Any], mech: Dict[str, float], modality_class: str) -> List[str]:
    recs: List[str] = []

    if mech.get("tlr78_risk_proxy", 0) >= 0.45:
        recs.append("PBMC cytokine release / innate immune activation panel")
    if mech.get("platelet_ps_risk_proxy", 0) >= 0.45:
        recs.append("Platelet activation and PF4-associated thrombocytopenia follow-up")
    if mech.get("hepatic_gapmer_risk_proxy", 0) >= 0.45:
        recs.append("Primary hepatocyte / hepatocellular injury follow-up")
    if mech.get("renal_accumulation_proxy", 0) >= 0.45:
        recs.append("Renal proximal tubule accumulation / injury follow-up")
    if modality_class in ("siRNA", "siRNA-like", "GalNAc-siRNA"):
        recs.append("Seed-mediated off-target transcriptomics review")
    if modality_class in ("ASO", "Gapmer/ASO"):
        recs.append("RNase H1 / hybridization off-target assessment")

    dominant = str(result.get("dominant_endpoint", "")).lower()
    if "complement" in dominant:
        recs.append("Complement activation biomarker follow-up")
    if "coag" in dominant or "platelet" in dominant or "thrombo" in dominant:
        recs.append("Coagulation panel and platelet safety review")
    if "renal" in dominant or "kidney" in dominant:
        recs.append("Kidney injury biomarkers and renal histopath follow-up")
    if "hep" in dominant or "liver" in dominant:
        recs.append("ALT / AST and hepatocyte stress biomarker follow-up")

    omics = result.get("omics", {})
    if isinstance(omics, dict):
        for rec in omics.get("recommended_assays", []) or []:
            if rec:
                recs.append(str(rec))

    deduped = []
    for r in recs:
        if r not in deduped:
            deduped.append(r)
    return deduped[:10]


def generate_regulatory_report(
    seq_clean: str,
    modality_class: str,
    result: Dict[str, Any],
    mech: Dict[str, float],
    registry: Dict[str, Any],
    recs: List[str],
) -> str:
    dominant = result.get("dominant_endpoint", "Not reported")
    return f"""
# OligoTox Regulatory-Style Safety Summary

## 1. Candidate Context
- Sequence: {seq_clean}
- Sequence length: {mech.get('length_nt', 'NA')} nt
- Modality class: {modality_class}
- Dominant endpoint: {dominant}

## 2. Mechanistic Decomposition
- D1 (on-target pharmacology): review target/context-specific liability
- D2 (hybridization off-target proxy): {mech.get('hybridization_risk_proxy', 'NA')}
- D3 (chemistry/delivery proxy): PS={mech.get('ps_fraction_proxy', 'NA')}, GalNAc={mech.get('galnac_flag', 'NA')}

## 3. Intermediate Risk Signals
- TLR7/8 innate immune proxy: {mech.get('tlr78_risk_proxy', 'NA')}
- Platelet / thrombocytopenia proxy: {mech.get('platelet_ps_risk_proxy', 'NA')}
- Hepatic gapmer proxy: {mech.get('hepatic_gapmer_risk_proxy', 'NA')}
- Renal accumulation proxy: {mech.get('renal_accumulation_proxy', 'NA')}

## 4. Engine and Benchmark Provenance
- Toxicity engine: {registry.get('tox_engine', 'NA')}
- Efficacy engine: {registry.get('efficacy_engine', 'NA')}
- PK/PD bridge: {registry.get('pkpd_bridge', 'NA')}
- Benchmarks: {", ".join(registry.get('primary_benchmarks', []))}
- Backend status note: {registry.get('backend_status', 'NA')}

## 5. Recommended Follow-Up Assays
{chr(10).join(f"- {x}" for x in recs) if recs else "- No additional assays suggested"}

## 6. Regulatory Framing
- OSWG-style interpretation: integrate in silico complementarity, cell context, in vitro verification, biological role assessment, and unavoidable off-target management
- FDA-style interpretation: align candidate risk with nonclinical safety package planning and modality-specific follow-up

## 7. Important Interpretation Note
- Trained model outputs come from the OligoTox backend.
- Intermediate mechanistic signals shown in the UI are interpretation/triage overlays and are not themselves validated clinical endpoints.
""".strip()


def flatten_endpoint_scores(result: Dict[str, Any]) -> pd.DataFrame:
    rows = []

    if "per_endpoint" in result and isinstance(result["per_endpoint"], dict):
        for endpoint, info in result["per_endpoint"].items():
            if isinstance(info, dict):
                prob = info.get("prob", info.get("risk_score", info.get("score", 0.5)))
                label = info.get("label", "")
                row = {"Endpoint": endpoint, "Risk Score": safe_float(prob, 0.5), "Label": label}
                for k, v in info.items():
                    if k not in row:
                        row[k] = v
                rows.append(row)
            else:
                rows.append({"Endpoint": endpoint, "Risk Score": safe_float(info, 0.5)})
    elif "risk_scores" in result and isinstance(result["risk_scores"], dict):
        for endpoint, prob in result["risk_scores"].items():
            rows.append({"Endpoint": endpoint, "Risk Score": safe_float(prob, 0.5)})

    if not rows:
        return pd.DataFrame()

    df = pd.DataFrame(rows)
    if "Risk Score" in df.columns:
        df = df.sort_values("Risk Score", ascending=False).reset_index(drop=True)
    return df


def compare_candidates_df(candidates: List[Dict[str, str]]) -> pd.DataFrame:
    rows = []
    for c in candidates:
        mech = compute_intermediate_safety_metrics(
            c.get("sequence", ""),
            c.get("chemistry", ""),
            c.get("conjugate", ""),
            c.get("modality", ""),
        )
        rows.append(
            {
                "Name": c.get("name", "candidate"),
                "Modality": c.get("modality", ""),
                "Chemistry": c.get("chemistry", ""),
                "Conjugate": c.get("conjugate", ""),
                "Length (nt)": mech.get("length_nt", 0),
                "GC Fraction": mech.get("gc_fraction", 0),
                "Hybridization Risk": mech.get("hybridization_risk_proxy", 0),
                "TLR7/8 Risk": mech.get("tlr78_risk_proxy", 0),
                "Platelet/PS Risk": mech.get("platelet_ps_risk_proxy", 0),
                "Hepatic Gapmer Risk": mech.get("hepatic_gapmer_risk_proxy", 0),
                "Renal Accumulation": mech.get("renal_accumulation_proxy", 0),
            }
        )
    return pd.DataFrame(rows)


def summarize_runs_df(runs: List[Dict[str, Any]]) -> pd.DataFrame:
    rows = []
    for r in runs:
        result = r.get("result", {}) if isinstance(r.get("result"), dict) else {}
        mech = compute_intermediate_safety_metrics(
            r.get("sequence", ""),
            r.get("chemistry", ""),
            r.get("conjugate", ""),
            r.get("modality", ""),
        )
        rows.append(
            {
                "Run ID": r.get("run_id", ""),
                "Timestamp": r.get("timestamp", ""),
                "Sequence": r.get("sequence", ""),
                "Length": len(r.get("sequence", "") or ""),
                "Modality": r.get("modality", ""),
                "Chemistry": r.get("chemistry", ""),
                "Conjugate": r.get("conjugate", ""),
                "Dominant Endpoint": result.get("dominant_endpoint", ""),
                "Hybridization Proxy": mech.get("hybridization_risk_proxy", 0.0),
                "TLR7/8 Proxy": mech.get("tlr78_risk_proxy", 0.0),
                "Platelet/PS Proxy": mech.get("platelet_ps_risk_proxy", 0.0),
            }
        )
    return pd.DataFrame(rows)


def render_pkpd_panel(result: Dict[str, Any], fallback_pkpd: Dict[str, float] | None = None) -> None:
    pkpd_metrics = result.get("pkpd_metrics", {}) if isinstance(result, dict) else {}
    if not pkpd_metrics and fallback_pkpd:
        pkpd_metrics = fallback_pkpd
    if not pkpd_metrics:
        return

    st.markdown("### Exposure Metrics (if dose/schedule provided)")
    c1, c2, c3 = st.columns(3)
    c1.metric("Liver AUC (pmol·g⁻¹·h)", f"{safe_float(pkpd_metrics.get('liver_auc_pmol_g_h')):,.1f}")
    c2.metric("Liver Cmax (pmol/g)", f"{safe_float(pkpd_metrics.get('liver_cmax_pmol_g')):,.1f}")
    c3.metric("Liver Half-life (h)", f"{safe_float(pkpd_metrics.get('liver_half_life_h')):,.1f}")

    extra_cols = st.columns(2)
    if "peak_knockdown_pct" in pkpd_metrics:
        extra_cols[0].metric("Peak Knockdown (%)", f"{safe_float(pkpd_metrics.get('peak_knockdown_pct')):.1f}")
    if "exposure_score" in pkpd_metrics:
        extra_cols[1].metric("Exposure Score", f"{safe_float(pkpd_metrics.get('exposure_score')):.2f}")

    with st.expander("Raw PK/PD payload"):
        st.json(pkpd_metrics)

    st.caption("PK/PD metrics are shown when dose and schedule are provided and PK/PD bridging is available.")


def render_confidence_panel(result: Dict[str, Any]) -> None:
    st.markdown("### Prediction Confidence")
    conf_rows = []
    endpoint_df = flatten_endpoint_scores(result)

    if not endpoint_df.empty:
        for _, row in endpoint_df.iterrows():
            prob = safe_float(row.get("Risk Score", 0.5), 0.5)
            conf_rows.append(
                {
                    "Endpoint": row.get("Endpoint", ""),
                    "Risk Score": round(prob, 4),
                    "Confidence": confidence_label(prob),
                    "Uncertainty Band": uncertainty_band(prob),
                }
            )

    if conf_rows:
        st.dataframe(pd.DataFrame(conf_rows), use_container_width=True)
        st.caption(
            "Current confidence is a probability-distance proxy. Replace with calibrated uncertainty "
            "intervals when model-level variance or ensemble dispersion is available."
        )
    else:
        st.info("No endpoint-level probability data available for confidence display.")


def render_modality_result_card(modality_class: str, registry: Dict[str, Any]) -> None:
    st.markdown("### Modality & Engine Context")
    c1, c2 = st.columns(2)
    with c1:
        st.info(
            f"**Modality class:** {modality_class}\n\n"
            f"**Toxicity engine:** {registry.get('tox_engine', 'N/A')}\n\n"
            f"**Efficacy / sequence engine:** {registry.get('efficacy_engine', 'N/A')}"
        )
    with c2:
        st.info(
            f"**PK/PD bridge:** {registry.get('pkpd_bridge', 'N/A')}\n\n"
            f"**Mechanistic layer:** {registry.get('mechanistic_layer', 'N/A')}\n\n"
            f"**Primary benchmarks:** {', '.join(registry.get('primary_benchmarks', []))}"
        )
    st.caption(registry.get("backend_status", ""))
    with st.expander("Benchmark & provenance support"):
        st.json(registry)


def render_intermediate_outputs(mech: Dict[str, float]) -> None:
    st.markdown("### Intermediate Safety Calculations")
    top = st.columns(5)
    top[0].metric("GC Fraction", mech.get("gc_fraction", 0))
    top[1].metric("U Fraction", mech.get("u_fraction", 0))
    top[2].metric("CpG Count", mech.get("cpg_count", 0))
    top[3].metric("Hybridization Risk", mech.get("hybridization_risk_proxy", 0))
    top[4].metric("TLR7/8 Risk", mech.get("tlr78_risk_proxy", 0))

    bottom = st.columns(4)
    bottom[0].metric("Platelet/PS Risk", mech.get("platelet_ps_risk_proxy", 0))
    bottom[1].metric("Hepatic Gapmer Risk", mech.get("hepatic_gapmer_risk_proxy", 0))
    bottom[2].metric("Renal Accumulation", mech.get("renal_accumulation_proxy", 0))
    bottom[3].metric("Length (nt)", mech.get("length_nt", 0))

    with st.expander("Intermediate calculation table"):
        st.dataframe(pd.DataFrame([mech]), use_container_width=True)

    st.caption("These are mechanistic-support proxy calculations for interpretation and triage, not final validated clinical scores.")


def render_assay_recommendations(recs: List[str]) -> None:
    st.markdown("### Assay Recommendations")
    if recs:
        for rec in recs:
            st.markdown(f"- {rec}")
    else:
        st.info("No assay recommendations generated.")


def render_last_result_snapshot(last: Dict[str, Any]) -> None:
    result = last.get("result", {}) if isinstance(last.get("result"), dict) else {}
    endpoint_df = flatten_endpoint_scores(result)
    st.markdown("### Latest Prediction Snapshot")
    st.markdown(f"**Run ID:** `{last.get('run_id', '')}`")
    st.markdown(
        f"**Sequence:** `{last.get('sequence', '')}` | "
        f"**Modality:** {last.get('modality', '')} | "
        f"**Chemistry:** {last.get('chemistry', '')} | "
        f"**Conjugate:** {last.get('conjugate', '') or 'None'}"
    )
    if result.get("dominant_endpoint"):
        st.info(f"🎯 Dominant endpoint: {result.get('dominant_endpoint')}")
    if not endpoint_df.empty:
        st.dataframe(endpoint_df, use_container_width=True)
    else:
        st.json(result)


def dataset_summary(dataset_path: Path) -> Dict[str, Any]:
    info: Dict[str, Any] = {"exists": dataset_path.exists(), "path": str(dataset_path)}
    if not dataset_path.exists():
        return info

    df = pd.read_csv(dataset_path)
    info["rows"] = len(df)
    info["columns"] = len(df.columns)
    info["column_names"] = list(df.columns)
    if "sequence" in df.columns:
        seq_series = df["sequence"].dropna().astype(str)
        info["sequence_coverage"] = int(seq_series.shape[0])
        if not seq_series.empty:
            info["median_length"] = int(seq_series.map(len).median())
            info["min_length"] = int(seq_series.map(len).min())
            info["max_length"] = int(seq_series.map(len).max())
    return info


def model_registry_df() -> pd.DataFrame:
    rows = []
    if MODEL_ARTIFACTS_DIR.exists():
        for p in sorted(list(MODEL_ARTIFACTS_DIR.glob("*.pkl")) + list(MODEL_ARTIFACTS_DIR.glob("*.joblib"))):
            rows.append(
                {
                    "Artifact": p.name,
                    "Size (KB)": round(p.stat().st_size / 1024, 2),
                    "Modified": datetime.fromtimestamp(p.stat().st_mtime).isoformat(timespec="seconds"),
                    "Path": str(p),
                }
            )
    return pd.DataFrame(rows)


def metrics_summary() -> Dict[str, Any]:
    metrics_file = METRICS_ARTIFACTS_DIR / "training_summary.json"
    out = {"metrics_file": str(metrics_file), "exists": metrics_file.exists()}
    if metrics_file.exists():
        try:
            with open(metrics_file, "r") as f:
                out["payload"] = json.load(f)
        except Exception as e:
            out["error"] = str(e)
    return out


st.set_page_config(
    page_title="OligoTox Platform",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.title(APP_TITLE)
st.caption(APP_CAPTION)

with st.sidebar:
    st.header("Platform Status")
    st.success("Backend entrypoint: run_oligotox")
    st.write(f"PK/PD bridge available: {'Yes' if PKPD_AVAILABLE else 'No'}")
    st.write(f"Processed data dir: `{PROCESSED_DATA_DIR}`")
    st.write(f"Model artifacts dir: `{MODEL_ARTIFACTS_DIR}`")
    st.write(f"Metrics dir: `{METRICS_ARTIFACTS_DIR}`")

tabs = st.tabs(TAB_NAMES)

# ----------------------------------------------------------------------------
# TAB 0: SEQUENCE ANALYSIS
# ----------------------------------------------------------------------------
with tabs[0]:
    st.subheader("Sequence Analysis")
    st.caption("Enter an oligonucleotide sequence for full toxicity risk profiling using the trained OligoTox backend.")

    col1, col2 = st.columns([2, 1])

    with col1:
        seq_input = st.text_area(
            "Input Sequence (RNA/DNA)",
            height=120,
            placeholder="e.g., AUGCUAGCUAGCUAGCUAGC",
            help="Paste ASO, siRNA, or other oligo sequence (ACGTUN alphabet; DNA T auto-normalized to U).",
        )

    with col2:
        modality = st.selectbox("Modality", ["ASO", "siRNA", "lncRNA-targeting", "Other"], index=0)
        chemistry = st.selectbox("Chemistry", ["PS", "PO", "2'-OMe", "2'-MOE", "LNA", "Mixed"], index=0)
        conjugate = st.selectbox("Conjugate", ["None", "GalNAc", "Cholesterol", "Other"], index=0)

    st.markdown("#### Optional PK/PD Context")
    pk1, pk2 = st.columns(2)
    with pk1:
        dose_mg_kg_input = st.number_input(
            "Dose (mg/kg)",
            min_value=0.0,
            value=0.0,
            step=0.1,
            help="Optional: enables PK/PD exposure display when > 0.",
        )
    with pk2:
        schedule_input = st.selectbox(
            "Schedule",
            ["", "QW", "Q2W", "Q4W", "Q8W", "Q12W"],
            index=0,
            help="Optional: choose a dosing schedule for PK/PD exposure context.",
        )

    with st.expander("Input preview"):
        seq_norm_preview = normalize_sequence(seq_input) if seq_input else ""
        st.write(
            {
                "normalized_sequence": seq_norm_preview,
                "length_nt": len(seq_norm_preview),
                "modality": modality,
                "chemistry": chemistry,
                "conjugate": conjugate,
                "dose_mg_kg": dose_mg_kg_input,
                "schedule": schedule_input,
            }
        )

    run_btn = st.button("🔬 Run Full OligoTox Prediction", type="primary", use_container_width=True)

    if run_btn:
        seq_clean = normalize_sequence(seq_input)
        is_valid, err = validate_sequence(seq_clean)

        if not is_valid:
            st.error(err)
        else:
            try:
                with st.spinner("Running OligoTox inference backend..."):
                    inp = {
                        "sequence": seq_clean,
                        "modality": modality,
                        "chemistry": chemistry,
                        "conjugate": conjugate if conjugate != "None" else "",
                        "length": len(seq_clean),
                    }

                    result = run_oligotox(inp, omics_payload=None)

                    fallback_pkpd = build_pkpd_metrics(
                        sequence=seq_clean,
                        modality=modality,
                        chemistry=chemistry,
                        dose_mg_kg=dose_mg_kg_input if dose_mg_kg_input > 0 else None,
                        schedule=schedule_input if schedule_input else None,
                    )

                modality_class = detect_modality_class(modality, chemistry, conjugate, seq_clean)
                registry = engine_dataset_registry(modality_class)
                mech = compute_intermediate_safety_metrics(seq_clean, chemistry, conjugate, modality)
                assay_recs = build_assay_recommendations(result, mech, modality_class)
                regulatory_report_md = generate_regulatory_report(
                    seq_clean,
                    modality_class,
                    result,
                    mech,
                    registry,
                    assay_recs,
                )

                st.success("✅ Prediction complete!")

                log_run(
                    st,
                    sequence=seq_clean,
                    modality=modality,
                    chemistry=chemistry,
                    conjugate=conjugate if conjugate != "None" else "",
                    inputs=inp,
                    result=result,
                    regulatory_md=regulatory_report_md,
                )

                st.markdown("### Toxicity Risk Profile")
                st.markdown(f"**Sequence:** `{seq_clean}` ({len(seq_clean)} nt)")
                st.markdown(f"**Modality:** {modality} | **Chemistry:** {chemistry} | **Conjugate:** {conjugate}")

                st.markdown("---")
                st.markdown("#### Per-Endpoint Risk Scores")

                endpoint_df = flatten_endpoint_scores(result)
                if not endpoint_df.empty:
                    st.dataframe(endpoint_df, use_container_width=True)
                else:
                    st.json(result)

                st.download_button(
                    label="📥 Download Full Result (JSON)",
                    data=json.dumps(result, indent=2, default=str),
                    file_name=f"oligotox_result_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json",
                    mime="application/json",
                )

                if "dominant_endpoint" in result:
                    st.info(f"🎯 **Dominant Risk:** {result['dominant_endpoint']}")

                st.markdown("---")
                render_modality_result_card(modality_class, registry)

                st.markdown("---")
                render_pkpd_panel(result, fallback_pkpd=fallback_pkpd)

                st.markdown("---")
                render_intermediate_outputs(mech)

                st.markdown("---")
                render_confidence_panel(result)

                st.markdown("---")
                render_assay_recommendations(assay_recs)

                st.markdown("---")
                st.markdown("### Regulatory-Style Safety Summary")
                st.markdown(regulatory_report_md)
                st.download_button(
                    label="📥 Download Regulatory Summary (MD)",
                    data=regulatory_report_md,
                    file_name=f"oligotox_regulatory_summary_{datetime.now().strftime('%Y%m%d_%H%M%S')}.md",
                    mime="text/markdown",
                )

                st.markdown("---")
                st.markdown("#### Mechanistic Attribution (D1/D2/D3)")
                st.markdown(
                    """
- **D1 (On-target):** Exaggerated pharmacology from intended target engagement  
- **D2 (Hybridization off-target):** Unintended RNA binding via sequence complementarity  
- **D3 (Chemistry/delivery):** Chemical-class or delivery-vehicle liabilities  
                    """
                )

                if "omics" in result and isinstance(result["omics"], dict):
                    if "adjustment_summary" in result["omics"]:
                        st.markdown(f"**Mechanistic Notes:** {result['omics']['adjustment_summary']}")
                    explanation = result["omics"].get("explanation", {})
                    if isinstance(explanation, dict) and explanation.get("summary"):
                        st.caption(explanation["summary"])

                with st.expander("Raw backend payload"):
                    st.json(result)

            except Exception as e:
                st.error(f"Prediction failed: {str(e)}")
                st.exception(e)

# ----------------------------------------------------------------------------
# TAB 1: OVERVIEW
# ----------------------------------------------------------------------------
with tabs[1]:
    st.subheader("Platform Overview")

    ds_path = PROCESSED_DATA_DIR / "training_dataset.csv"
    model_df = model_registry_df()
    metric_info = metrics_summary()

    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Training Dataset Present", "Yes" if ds_path.exists() else "No")
    c2.metric("Model Artifacts", len(model_df))
    c3.metric("Metrics Summary", "Yes" if metric_info.get("exists") else "No")
    c4.metric("PK/PD Bridge", "Available" if PKPD_AVAILABLE else "Unavailable")

    st.markdown("---")
    st.markdown(
        """
### OligoTox Platform

**Mission:** Mechanism-aware toxicity prediction for oligonucleotide therapeutics.

**Core Framework: D1/D2/D3**
- **D1:** On-target pharmacology
- **D2:** Hybridization-dependent off-target effects
- **D3:** Chemistry and delivery liabilities

**Current App Scope**
- Trained OligoTox backend inference via `run_oligotox`
- Session-based run logging and comparison
- Optional PK/PD exposure bridge
- Mechanistic proxy overlays for interpretation
- Regulatory-style summary generation
- Artifact-aware model and metrics inspection
        """
    )

# ----------------------------------------------------------------------------
# TAB 2: CURATION
# ----------------------------------------------------------------------------
with tabs[2]:
    st.subheader("Data Curation Status")

    try:
        dataset_path = PROCESSED_DATA_DIR / "training_dataset.csv"
        info = dataset_summary(dataset_path)

        if info["exists"]:
            st.success(f"✅ Training dataset found: {info['path']}")
            c1, c2, c3, c4 = st.columns(4)
            c1.metric("Rows", info.get("rows", 0))
            c2.metric("Columns", info.get("columns", 0))
            c3.metric("Seq Coverage", info.get("sequence_coverage", 0))
            c4.metric("Median Length", info.get("median_length", "NA"))

            preview = pd.read_csv(dataset_path, nrows=15)
            st.dataframe(preview, use_container_width=True)

            with st.expander("Dataset metadata"):
                st.json(info)
        else:
            st.warning("Training dataset not found.")
    except Exception as e:
        st.error(f"Error loading dataset: {e}")
        st.exception(e)

# ----------------------------------------------------------------------------
# TAB 3: RUNS
# ----------------------------------------------------------------------------
with tabs[3]:
    st.subheader("Pipeline Runs")
    st.caption("Session-level run history generated from Sequence Analysis activity.")

    runs = get_runs(st)
    if runs:
        run_df = summarize_runs_df(runs)
        st.success(f"Found {len(run_df)} logged run(s) in session state.")
        st.dataframe(run_df, use_container_width=True)

        latest = get_last_run(st)
        if latest:
            with st.expander("Latest run details", expanded=True):
                render_last_result_snapshot(latest)
    else:
        st.info("No logged runs yet. Execute a Sequence Analysis run to populate this tab.")

# ----------------------------------------------------------------------------
# TAB 4: MODELS
# ----------------------------------------------------------------------------
with tabs[4]:
    st.subheader("Model Registry")

    try:
        df_models = model_registry_df()
        if not df_models.empty:
            st.success(f"Found {len(df_models)} model artifact(s)")
            st.dataframe(df_models, use_container_width=True)
        else:
            st.info("No model binaries found. Train a model first.")

        with st.expander("Artifact directories"):
            st.write(
                {
                    "MODEL_ARTIFACTS_DIR": str(MODEL_ARTIFACTS_DIR),
                    "METRICS_ARTIFACTS_DIR": str(METRICS_ARTIFACTS_DIR),
                }
            )
    except Exception as e:
        st.error(f"Error reading model registry: {e}")
        st.exception(e)

# ----------------------------------------------------------------------------
# TAB 5: PREDICTIONS
# ----------------------------------------------------------------------------
with tabs[5]:
    st.subheader("Recent Predictions")
    st.caption("Detailed prediction outputs from this session.")

    runs = get_runs(st)
    if runs:
        run_options = [r.get("run_id", f"run_{i}") for i, r in enumerate(runs)]
        selected_run_id = st.selectbox("Select a run", run_options, index=len(run_options) - 1)
        selected_run = next((r for r in runs if r.get("run_id") == selected_run_id), runs[-1])

        result = selected_run.get("result", {})
        endpoint_df = flatten_endpoint_scores(result)

        st.markdown(
            f"**Sequence:** `{selected_run.get('sequence', '')}` | "
            f"**Modality:** {selected_run.get('modality', '')} | "
            f"**Chemistry:** {selected_run.get('chemistry', '')} | "
            f"**Conjugate:** {selected_run.get('conjugate', '') or 'None'}"
        )

        if result.get("dominant_endpoint"):
            st.info(f"🎯 Dominant endpoint: {result.get('dominant_endpoint')}")

        if not endpoint_df.empty:
            st.dataframe(endpoint_df, use_container_width=True)
        else:
            st.json(result)

        st.download_button(
            label="📥 Download Selected Prediction JSON",
            data=json.dumps(result, indent=2, default=str),
            file_name=f"{selected_run_id}_prediction.json",
            mime="application/json",
        )
    else:
        st.info("No prediction runs logged yet. Execute Sequence Analysis first.")

# ----------------------------------------------------------------------------
# TAB 6: EVALUATION
# ----------------------------------------------------------------------------
with tabs[6]:
    st.subheader("📊 Model Evaluation & Clinical Metrics")
    st.caption("Rigorous evaluation with leakage controls, calibration, and uncertainty quantification")

    st.markdown(
        """
### Evaluation Philosophy

1. **Leakage-Resistant Splitting**
   - Sequence similarity
   - Batch effects
   - Temporal structure
   - Study/assay families

2. **Clinically Meaningful Metrics**
   - AUROC and AUPRC
   - Brier score
   - Per-endpoint thresholds
   - Reliability-oriented tracking

3. **Uncertainty Quantification**
   - Confidence bands
   - Coverage estimation
   - OOD-aware future extensions
        """
    )

    st.markdown("---")

    try:
        metrics_file = METRICS_ARTIFACTS_DIR / "training_summary.json"
        if metrics_file.exists():
            with open(metrics_file, "r") as f:
                metrics = json.load(f)

            st.success("✅ Model evaluation metrics found")
            st.markdown("### Global Metrics")

            col1, col2, col3 = st.columns(3)
            if "composite_auroc" in metrics:
                col1.metric("Composite AUROC", f"{safe_float(metrics['composite_auroc']):.3f}")
            if "composite_brier" in metrics:
                col2.metric("Composite Brier Score", f"{safe_float(metrics['composite_brier']):.3f}")
            if "n_endpoints" in metrics:
                col3.metric("Evaluated Endpoints", metrics["n_endpoints"])

            st.markdown("---")
            st.markdown("### Per-Endpoint Performance")

            if "per_endpoint" in metrics and isinstance(metrics["per_endpoint"], dict):
                endpoint_rows = []
                for ep, vals in metrics["per_endpoint"].items():
                    endpoint_rows.append(
                        {
                            "Endpoint": ep,
                            "AUROC": vals.get("auroc", "N/A"),
                            "AUPRC": vals.get("auprc", "N/A"),
                            "Brier": vals.get("brier", "N/A"),
                            "Threshold": vals.get("threshold", "N/A"),
                        }
                    )
                df_metrics = pd.DataFrame(endpoint_rows)
                st.dataframe(df_metrics, use_container_width=True)

            st.markdown("---")
            st.markdown("### Leakage Diagnostic")

            datasail_file = METRICS_ARTIFACTS_DIR / "datasail_split_diagnostic.json"
            if datasail_file.exists():
                with open(datasail_file, "r") as f:
                    datasail_info = json.load(f)
                st.json(datasail_info)
            else:
                st.info("DataSAIL leakage diagnostic not found. Run training pipeline to generate.")
        else:
            st.warning("⚠️ No evaluation metrics found. Train the model to generate metrics.")
            st.code(
                "python train.py",
                language="bash",
            )
    except Exception as e:
        st.error(f"Error loading evaluation metrics: {e}")
        st.exception(e)

    st.markdown("---")
    st.markdown("### Calibration & Uncertainty")
    st.markdown(
        """
**Calibration:** Does a predicted 70% risk correspond to 70% observed frequency?  
**Uncertainty:** How confident is the model? Does it know when it doesn't know?

Current UI confidence labels are probability-distance proxies; formal calibrated uncertainty should come from model-level calibration or ensemble variance.
        """
    )

# ----------------------------------------------------------------------------
# TAB 7: COMPARATOR
# ----------------------------------------------------------------------------
with tabs[7]:
    st.subheader("Candidate Comparator")
    st.caption("Compare multiple oligonucleotide candidates on sequence-derived mechanistic overlays and session runs.")

    runs = get_runs(st)
    if runs:
        candidates = [
            {
                "name": r.get("run_id", "run"),
                "sequence": r.get("sequence", ""),
                "modality": r.get("modality", ""),
                "chemistry": r.get("chemistry", ""),
                "conjugate": r.get("conjugate", ""),
            }
            for r in runs
        ]
        df_cmp = compare_candidates_df(candidates)
        st.dataframe(df_cmp, use_container_width=True)
        st.caption("Comparator is populated from logged session runs.")
    else:
        sample_candidates = [
            {
                "name": "Candidate A",
                "sequence": "AUGCUAGCUAGCUAGCUAGC",
                "modality": "ASO",
                "chemistry": "PS",
                "conjugate": "None",
            },
            {
                "name": "Candidate B",
                "sequence": "AUGUUAGUUAGCUAGUUAGU",
                "modality": "siRNA",
                "chemistry": "2'-OMe",
                "conjugate": "GalNAc",
            },
            {
                "name": "Candidate C",
                "sequence": "GCUAUGCUUGGCUAUGCUUA",
                "modality": "ASO",
                "chemistry": "LNA",
                "conjugate": "None",
            },
        ]
        df_cmp = compare_candidates_df(sample_candidates)
        st.dataframe(df_cmp, use_container_width=True)
        st.caption("No logged runs yet; showing sample candidates.")

# ----------------------------------------------------------------------------
# TAB 8: ASSAYS
# ----------------------------------------------------------------------------
with tabs[8]:
    st.subheader("Assay Recommendation Center")
    st.caption("Mechanism-linked follow-up assay suggestions.")

    last = get_last_run(st)
    if last and isinstance(last.get("result"), dict):
        seq = last.get("sequence", "")
        modality = last.get("modality", "")
        chemistry = last.get("chemistry", "")
        conjugate = last.get("conjugate", "")
        result = last.get("result", {})
        modality_class = detect_modality_class(modality, chemistry, conjugate, seq)
        mech = compute_intermediate_safety_metrics(seq, chemistry, conjugate, modality)
        recs = build_assay_recommendations(result, mech, modality_class)

        render_assay_recommendations(recs)
        st.caption("Showing assay recommendations from the latest Sequence Analysis run.")
    else:
        st.markdown(
            """
- PBMC cytokine release / innate immune activation
- Platelet activation / PF4 follow-up
- Hepatocyte injury and transcriptomics follow-up
- Renal proximal tubule accumulation / injury follow-up
- Seed-mediated off-target transcriptomics
- RNase H1 / hybridization off-target confirmation
- Complement activation biomarker panel
- Coagulation and thrombocytopenia follow-up
            """
        )
        st.caption("No run-linked assay recommendations yet; showing default panel.")

# ----------------------------------------------------------------------------
# TAB 9: REGULATORY REPORT
# ----------------------------------------------------------------------------
with tabs[9]:
    last = get_last_run(st)
    if last and last.get("regulatory_md"):
        st.subheader("OSWG / FDA Run-Level Report")
        st.caption("Generated from the latest Sequence Analysis run.")
        st.markdown(last["regulatory_md"])
        st.download_button(
            label="📥 Download Latest Regulatory Summary (MD)",
            data=last["regulatory_md"],
            file_name="latest_oligotox_regulatory_summary.md",
            mime="text/markdown",
        )

        with st.expander("Latest run payload"):
            st.json(last.get("result", {}))
    else:
        st.subheader("OSWG / FDA Report Framing")
        st.caption("Run-level output should map predictions into a structured preclinical safety narrative.")
        st.markdown(
            """
### Recommended Report Sections
1. Candidate context and modality
2. D1 / D2 / D3 mechanistic interpretation
3. Endpoint risks and uncertainty
4. Intermediate calculations and exposure context
5. Assay recommendations
6. OSWG off-target workflow interpretation
7. FDA-oriented nonclinical safety planning notes
            """
        )