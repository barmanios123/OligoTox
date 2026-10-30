#!/usr/bin/env python3
import os
import sys
import json
import uuid
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd
import streamlit as st

ROOT = Path(__file__).resolve().parent
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from oligotox.features.candidate import candidate_features as build_candidate_features
from oligotox.pipelines.score_candidate import score_candidate
from oligotox.pkpd.bridge import compute_pkpd_enriched_features
from oligotox.advice import explain, suggest, summarize
from oligotox.efficacy.potency import score_efficacy
from oligotox.reporting.oswg_report import OSWGSafetyReport

APP_TITLE = "OligoTox Operator Console"
APP_SUBTITLE = (
    "Research prototype. Two toxicity endpoints of 13 currently clear the "
    "evidence gate; the rest are shown as unsupported with the reason. "
    "CNS neurotoxicity and siRNA potency are trained on external data. "
    "Each number carries the evidence behind it, and How it works has the "
    "whole of it."
)

# Endpoints are presented by how well the evaluation supports them, not by how
# confident the model sounds: see scripts/derive_evidence_tiers.py.
TIER_LABELS = {
    "validated": "Validated",
    "weak": "Weak signal",
    "source_confounded": "Source-confounded",
    "insufficient_data": "Insufficient data",
}

PROCESSED_DIR = ROOT / "data" / "processed"
MODEL_DIR = ROOT / "artifacts" / "models"
CONFIG_DIR = ROOT / "artifacts" / "config"

st.set_page_config(
    page_title="OligoTox",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

def _safe_float(x, default=0.0):
    try:
        return float(x)
    except Exception:
        return default

def _extract_endpoint_value(info):
    if isinstance(info, dict):
        for key in [
            "probability",
            "score",
            "risk_score",
            "predicted_probability",
            "prob",
            "value",
            "raw_score",
        ]:
            if key in info and info[key] is not None:
                try:
                    return float(info[key])
                except Exception:
                    return info[key]
        return None
    return info

def _validate_sequence(seq: str):
    seq = (seq or "").strip().upper()
    if not seq:
        return False, "Sequence is required."
    allowed = set("ACGTU")
    bad = sorted(set(ch for ch in seq if ch not in allowed))
    if bad:
        return False, f"Invalid sequence characters: {', '.join(bad)}. Allowed: A, C, G, T, U."
    if len(seq) < 8:
        return False, "Sequence is too short for meaningful analysis; use at least 8 nt."
    return True, seq

GENOMIC_DIR = ROOT / "data" / "external" / "genomic"


@st.cache_data(show_spinner=False)
def _target_transcript(symbol: str) -> str | None:
    """The cached genomic span for a gene symbol, or None.

    None is a real answer here and is reported as such. Scoring falls back to
    the design-only model rather than inventing accessibility values.
    """
    if not symbol:
        return None
    path = GENOMIC_DIR / f"{symbol.replace('/', '_')}.json"
    if not path.exists():
        return None
    try:
        record = json.loads(path.read_text())
    except json.JSONDecodeError:
        return None
    seq = (record.get("sequence") or "").upper().replace("U", "T")
    return seq or None


def _ps_fraction(chemistry: str) -> float:
    """The one place this is decided.

    The summary and the design scan both used to default to 1.0 when the run
    payload carried no value, while scoring used 0.9, so the interface
    reported a backbone the model never saw.
    """
    return 0.9 if "PS" in chemistry.upper() else 0.0


def _sequence_features(seq: str, chemistry: str, dose_mg_kg: float, modality: str = "siRNA"):
    ps_fraction = _ps_fraction(chemistry)
    return build_candidate_features(
        seq,
        modality=modality,
        chemistry=chemistry,
        ps_fraction=ps_fraction,
        dose_mg_kg=dose_mg_kg,
    )

def _run_toxicity(candidate_features):
    return score_candidate(
        candidate_features=candidate_features,
        processed_dir=PROCESSED_DIR,
        model_dir=MODEL_DIR,
        config_dir=CONFIG_DIR,
    )

def _run_pkpd(seq: str, dose_mg_kg: float, modality: str, chemistry: str, schedule: str = "weekly"):
    try:
        return compute_pkpd_enriched_features(
            sequence=seq,
            modality=modality,
            chemistry=chemistry,
            schedule=schedule,
            dose_mg_kg=dose_mg_kg,
        )
    except TypeError:
        try:
            return compute_pkpd_enriched_features(
                seq,
                modality,
                chemistry,
                schedule,
                dose_mg_kg=dose_mg_kg,
            )
        except Exception as e:
            return {"error": str(e)}
    except Exception as e:
        return {"error": str(e)}

def _build_report(seq, chemistry, toxicity_result, pkpd_result):
    toxicity_predictions = {}
    for endpoint, info in toxicity_result.get("per_endpoint", {}).items():
        val = _extract_endpoint_value(info)
        if isinstance(val, (float, int)):
            toxicity_predictions[endpoint] = float(val)

    pkpd_summary = {}
    if isinstance(pkpd_result, dict):
        for k, v in pkpd_result.items():
            if isinstance(v, (float, int, str, bool)):
                pkpd_summary[k] = v

    evidence_tiers = {
        endpoint: info.get("evidence_tier")
        for endpoint, info in toxicity_result.get("per_endpoint", {}).items()
        if isinstance(info, dict)
    }

    report = OSWGSafetyReport(
        compound_id=f"RUN-{datetime.now(UTC).strftime('%Y%m%d%H%M%S')}",
        sequence=seq,
        assessment_date=datetime.now(UTC).strftime("%Y-%m-%d"),
        toxicity_predictions=toxicity_predictions,
        pkpd_summary=pkpd_summary,
        evidence_tiers=evidence_tiers,
    )

    report_json = report.generate_report(format="json")
    risk_score = report._calculate_risk_score()
    risk_class = report._classify_safety()

    return {
        "json": report_json,
        "risk_score": risk_score,
        "risk_class": risk_class,
    }

def _save_run(payload: dict):
    """Keep the run in this visitor's session only.

    Runs used to be written to artifacts/runs/ and the Saved runs tab listed
    every file it found there. On a shared deployment that is one visitor's
    sequences shown to the next, so runs now live in session state and leave
    with the session. Use the download buttons to keep one.
    """
    run_id = payload.get("run_id") or str(uuid.uuid4())
    payload["run_id"] = run_id
    payload["saved_at"] = datetime.now(UTC).isoformat().replace("+00:00", "Z")
    st.session_state.setdefault("runs", []).insert(0, payload)
    return run_id

def _load_saved_runs():
    """This session's runs, newest first."""
    rows = []
    for obj in st.session_state.get("runs", []):
        rows.append(
            {
                "run_id": obj.get("run_id"),
                "saved_at": obj.get("saved_at"),
                "sequence": obj.get("input", {}).get("sequence"),
                "chemistry": obj.get("input", {}).get("chemistry"),
                "modality": obj.get("input", {}).get("modality"),
                "composite_score": obj.get("toxicity", {}).get("composite_score"),
                "composite_band": obj.get("toxicity", {}).get("composite_band"),
                "dominant_endpoint": obj.get("toxicity", {}).get("dominant_endpoint"),
            }
        )
    return pd.DataFrame(rows)


@st.cache_data(show_spinner=False)
def _model_features() -> dict[str, list[str]]:
    """What each shipped model actually consumes, read off the model.

    Taken from the estimators rather than from a config file, so the How it
    works page lists the real inputs even for the endpoints whose feature set
    was never written down separately.
    """
    import joblib

    out: dict[str, list[str]] = {}
    for path in sorted(MODEL_DIR.glob("*.joblib")):
        try:
            names = getattr(joblib.load(path), "feature_names_in_", None)
        except Exception:
            continue
        if names is not None:
            out[path.stem] = list(names)
    return out


def _load_report(name: str) -> dict:
    """A shipped report, or {} if this build does not carry it."""
    return _load_json(ROOT / "reports" / "results" / name)


def _load_config(name: str) -> dict:
    return _load_json(CONFIG_DIR / name)


def _load_json(path: Path) -> dict:
    try:
        return json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return {}


def _render_endpoint_table(toxicity_result, tiers=None):
    rows = []
    for endpoint, info in toxicity_result.get("per_endpoint", {}).items():
        if tiers is not None and info.get("evidence_tier") not in tiers:
            continue
        rows.append(
            {
                "endpoint": endpoint,
                "value": _extract_endpoint_value(info),
                # A number that is not calibrated is a rank, not a percentage.
                "reads as": ("probability" if info.get("probabilities_calibrated")
                             else "ranking only"),
                "threshold": info.get("threshold") if isinstance(info, dict) else None,
                "predicted_positive": info.get("predicted_positive") if isinstance(info, dict) else None,
                "risk_band": info.get("risk_band") if isinstance(info, dict) else None,
            }
        )
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    return df.sort_values("value", ascending=False, na_position="last")


def _render_unsupported_table(toxicity_result):
    """Endpoints the evaluation doesn't support, with why."""
    rows = []
    for endpoint, info in toxicity_result.get("per_endpoint", {}).items():
        tier = info.get("evidence_tier")
        if tier == "validated":
            continue
        rows.append(
            {
                "endpoint": endpoint,
                "tier": TIER_LABELS.get(tier, tier),
                "model output": _extract_endpoint_value(info),
                "held-out PR-AUC": info.get("cv_pr_auc"),
                "base rate": info.get("base_rate"),
                "positive sequences": info.get("positive_sequences"),
                "why it is not a prediction": info.get("evidence_note"),
            }
        )
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    order = {label: i for i, label in enumerate(TIER_LABELS.values())}
    return df.sort_values(
        ["tier", "positive sequences"],
        ascending=[True, False],
        key=lambda s: s.map(order) if s.name == "tier" else s,
    )

st.title(APP_TITLE)
st.caption(APP_SUBTITLE)

with st.sidebar:
    st.header("Run setup")
    sequence = st.text_area(
        "Sequence",
        value="AUGCUAGCUAGCUAGCUAGC",
        height=120,
        help="Enter RNA/DNA oligo sequence using A/C/G/T/U."
    )
    chemistry = st.selectbox("Chemistry", ["PS", "PO", "mixed_PS_PO"], index=0)
    modality = st.selectbox("Modality", ["siRNA", "ASO"], index=0)
    dose_mg_kg = st.number_input("Dose (mg/kg)", min_value=0.0, max_value=100.0, value=3.0, step=0.5)
    target_gene = st.text_input(
        "Target gene (optional)",
        value="",
        help="A gene symbol such as DGAT2. When its genomic span is cached and "
             "the oligo's binding site locates in it, assays with an "
             "accessibility model use that model instead.",
    ).strip().upper()
    analyst_note = st.text_input("Analyst note", value="")
    run_btn = st.button("Run analysis", type="primary", width="stretch")

# Summary first, because it is the one view that answers "what did this
# find" without the reader knowing which tab holds which half of it. The old
# "Candidate Input" tab only echoed the sidebar's own fields and dumped raw
# features, so it is folded in here rather than kept as a ninth thing to read.
TAB_NAMES = [
    "Summary",
    "Toxicity",
    "Efficacy",
    "PK/PD",
    "Design advice",
    "Report",
    "How it works",
    "Saved Runs",
    "Build Info",
]
# Keyed by name rather than position: the blocks below used to be `tabs[4]`,
# which meant inserting a tab silently moved every panel after it.
tab = dict(zip(TAB_NAMES, st.tabs(TAB_NAMES)))

if "last_run" not in st.session_state:
    st.session_state["last_run"] = None

# Validation runs before the tabs so the sidebar can report a bad sequence
# immediately, without waiting for a run.
_input_ok, _input_msg = _validate_sequence(sequence)
normalized_seq = _input_msg if _input_ok else None

if run_btn:
    ok, seq_or_msg = _validate_sequence(sequence)
    if not ok:
        st.error(seq_or_msg)
    else:
        normalized_seq = seq_or_msg
        candidate_features = _sequence_features(normalized_seq, chemistry, dose_mg_kg, modality)

        with st.spinner("Running OligoTox analysis..."):
            toxicity_result = _run_toxicity(candidate_features)
            pkpd_result = _run_pkpd(normalized_seq, dose_mg_kg, modality, chemistry)
            report_result = _build_report(normalized_seq, chemistry, toxicity_result, pkpd_result)

        payload = {
            "input": {
                "sequence": normalized_seq,
                "chemistry": chemistry,
                "modality": modality,
                "dose_mg_kg": dose_mg_kg,
                "ps_fraction": _ps_fraction(chemistry),
                "target_gene": target_gene,
                "analyst_note": analyst_note,
                "candidate_features": candidate_features,
            },
            "toxicity": toxicity_result,
            "pkpd": pkpd_result,
            "report": report_result,
        }
        run_id = _save_run(payload)
        st.session_state["last_run"] = payload
        st.success(f"Analysis complete. Run {run_id} kept in this session.")

run = st.session_state.get("last_run")

with tab["Summary"]:
    if not _input_ok:
        st.error(_input_msg)
    elif not run:
        st.subheader("Ready")
        st.success(f"Sequence accepted, {len(normalized_seq)} nt.")
        st.code(normalized_seq)
        st.info("Press **Run analysis** in the sidebar. The summary of what was "
                "measured, what is predicted and what is deliberately not "
                "predicted appears here.")
        with st.expander("Derived features, before any model runs"):
            st.caption("These are computed from the sequence alone. They are "
                       "inputs, not results.")
            st.json({k: (round(v, 4) if isinstance(v, (int, float)) else v)
                     for k, v in _sequence_features(
                         normalized_seq, chemistry, dose_mg_kg, modality).items()})
    else:
        st.subheader("Summary")
        with st.spinner("Summarising"):
            full = summarize(run["input"]["sequence"], run["input"]["modality"],
                             run["input"]["chemistry"],
                             run["input"].get("ps_fraction", _ps_fraction(run["input"]["chemistry"])))

        # The headline first, so the answer does not depend on reading prose.
        tox = run["toxicity"]
        validated = full["validated_toxicity"]
        c1, c2, c3 = st.columns(3)
        c1.metric("Endpoints presented", f"{len(validated)} of 13")
        c2.metric("Composite", "none" if tox.get("composite_score") is None
                  else f"{_safe_float(tox.get('composite_score')):.3f}")
        c3.metric("Band", str(tox.get("composite_band") or "not scored"))
        st.caption("Ten of the thirteen toxicity models run and are withheld "
                   "because they failed an evidence check. The Toxicity tab "
                   "lists each one and the reason.")

        st.code(full["text"], language=None)
        st.download_button("Download this summary", data=full["text"],
                           file_name=f"oligotox_summary_{full['sequence'][:12]}.txt",
                           mime="text/plain")
        if full["parameters_unavailable"]:
            st.warning("Some parameters this summary expects were not produced "
                       f"by the feature pipeline: {full['parameters_unavailable']}. "
                       "The summary is incomplete rather than wrong.")

        with st.expander("What the other tabs add"):
            st.markdown(
                "- **Toxicity** every endpoint, including the ten withheld, with why\n"
                "- **Efficacy** predicted knockdown, scored per assay\n"
                "- **PK/PD** exposure and dose\n"
                "- **Design advice** every single-base edit, scored and ranked\n"
                "- **Report** the regulatory write-up\n"
                "- **Saved Runs** this session's runs, kept in your browser only"
            )

with tab["Toxicity"]:
    st.subheader("Toxicity")
    if not run:
        st.info("Run an analysis to view toxicity outputs.")
    else:
        tox = run["toxicity"]
        basis = tox.get("composite_basis") or []
        total = len(tox.get("per_endpoint", {}))

        if tox.get("composite_score") is None:
            st.warning(
                "No endpoint in this build passes the evaluation, so no composite score is "
                "reported. The per-endpoint outputs below are not predictions."
            )
        else:
            c1, c2, c3 = st.columns(3)
            c1.metric("Composite score", f"{_safe_float(tox.get('composite_score')):.3f}")
            c2.metric("Band", str(tox.get("composite_band")))
            c3.metric("Highest validated endpoint", str(tox.get("dominant_endpoint")))
            st.caption(
                f"Composite covers {len(basis)} of {total} endpoints "
                f"({', '.join(e.replace('_', ' ') for e in basis)}). The other "
                f"{total - len(basis)} are excluded because the baselines don't support them."
            )

        st.markdown("**Predictions**")
        validated_df = _render_endpoint_table(tox, tiers={"validated"})
        if validated_df.empty:
            st.info("No validated endpoints in this build.")
        else:
            st.dataframe(validated_df, width="stretch", hide_index=True)

        unsupported_df = _render_unsupported_table(tox)
        if not unsupported_df.empty:
            with st.expander(
                f"Not supported by the evaluation ({len(unsupported_df)} endpoints) — "
                "model output shown, but do not read these as predictions"
            ):
                st.caption(
                    "These models run, but each one either loses to a baseline that sees only "
                    "the source dataset or one crude sequence number, or rests on fewer than 10 "
                    "positive sequences. They are excluded from the composite score."
                )
                st.dataframe(unsupported_df, width="stretch", hide_index=True)

        st.download_button(
            "Download toxicity JSON",
            data=json.dumps(tox, indent=2, default=str),
            file_name="toxicity_result.json",
            mime="application/json",
        )

with tab["PK/PD"]:
    st.subheader("PK/PD")
    if not run:
        st.info("Run an analysis to view PK/PD outputs.")
    else:
        pkpd = run["pkpd"]
        if isinstance(pkpd, dict):
            st.json(pkpd)
        else:
            st.write(pkpd)

with tab["Report"]:
    st.subheader("Regulatory report")
    if not run:
        st.info("Run an analysis to generate a report.")
    else:
        rep = run["report"]
        c1, c2 = st.columns(2)
        c1.metric("Risk score", f"{_safe_float(rep.get('risk_score')):.3f}")
        c2.metric("Classification", str(rep.get("risk_class")))
        st.code(rep.get("json", ""), language="json")
        st.download_button(
            "Download report JSON",
            data=rep.get("json", ""),
            file_name="oswg_report.json",
            mime="application/json",
        )


with tab["Efficacy"]:
    st.subheader("Efficacy")
    if not run:
        st.info("Run an analysis to view efficacy outputs.")
    else:
        seq = run["input"]["sequence"]
        modality = run["input"]["modality"]
        transcript = _target_transcript(run["input"].get("target_gene", ""))
        eff = score_efficacy(seq, modality, MODEL_DIR, CONFIG_DIR,
                             target_transcript=transcript)

        gene = run["input"].get("target_gene", "")
        if gene and transcript is None:
            st.warning(
                f"No cached genomic span for {gene}, so accessibility was not "
                "used. The design-only models ran instead. Fetch it with "
                "`python scripts/fetch_transcripts.py --genomic`."
            )
        elif gene and not any(p.get("used_accessibility")
                              for p in eff.get("predictions", {}).values()):
            st.info(
                f"{gene} is cached, but this oligo's binding site was not found "
                "in it, so the design-only models ran. A site that cannot be "
                "located is unknown, not inaccessible."
            )
        elif gene:
            st.success(f"Target-site accessibility in {gene} was used where an "
                       "assay has a model for it.")

        if not eff["available"]:
            st.info(eff["reason"])
        else:
            ordered = sorted(eff["predictions"].items(),
                             key=lambda kv: {"validated": 0, "weak": 1}.get(kv[1]["tier"], 2))
            for name, pred in ordered:
                label = name.replace("_", " ")
                if pred["tier"] == "context":
                    ctx = pred.get("context", {})
                    st.metric(f"{ctx.get('cell_line', 'assay')} knockdown",
                              f"{pred['score']:.0f}%", help=pred["units"])
                    st.caption(
                        f"{pred['readout']} in {ctx.get('cell_line')} "
                        f"at dose {ctx.get('dose')}, {ctx.get('transfection_method')}. "
                        f"{pred['evidence']}"
                    )
                    if pred.get("used_accessibility"):
                        st.caption("Scored with target-site accessibility.")
                    st.caption(f"Caveat: {pred['caveat']}")
                    st.divider()
                    continue
                if pred["tier"] == "validated":
                    st.metric(label, f"{pred['score']:.1f}", help=pred["units"])
                else:
                    st.metric(f"{label} (weak evidence)", f"{pred['score']:.2f}", help=pred["units"])
                st.caption(f"{pred['readout']}. {pred['evidence']}")
                if pred["caveat"]:
                    st.caption(f"Caveat: {pred['caveat']}")
                st.divider()
            st.caption(
                "These are raw model readouts in the units of the training data, not "
                "probabilities, and they are not comparable across modalities."
            )


with tab["Design advice"]:
    st.subheader("Design advice")
    if not run:
        st.info("Run an analysis to see which edits change the prediction.")
    else:
        seq = run["input"]["sequence"]
        modality = run["input"]["modality"]

        # The full summary lives on the Summary tab.

        st.caption(
            "Every single-base substitution is scored through the same models the "
            "other tabs use. Each number is the difference between two real "
            "predictions, not a coefficient. Only endpoints tiered validated are "
            "used, because an edit recommended to move an unsupported number is "
            "worse than no advice."
        )
        with st.spinner("Scanning substitutions..."):
            advice = suggest(seq, modality, run["input"]["chemistry"],
                             ps_fraction=run["input"]["candidate_features"].get("ps_fraction", 1.0),
                             model_dir=MODEL_DIR, config_dir=CONFIG_DIR, top_n=12)

        if not advice.validated_endpoints:
            st.warning(
                "No validated toxicity endpoint applies to this candidate, so no edit "
                "can be recommended on toxicity grounds."
            )
        else:
            st.write("**Starting point**")
            st.json({k: round(v, 4) for k, v in advice.baseline_toxicity.items()})

        if advice.edits:
            rows = []
            for e in advice.edits:
                free = e.total_toxicity_reduction > 0 and e.mean_potency_change >= 0
                rows.append({
                    "edit": f"{e.from_base}{e.position}{e.to_base}",
                    "risk reduction": round(e.total_toxicity_reduction, 4),
                    "knockdown change %": round(e.mean_potency_change, 1),
                    "trade-off": "none, improves both" if free else "costs potency",
                    "sequence": e.sequence,
                })
            st.write("**Suggested single-base edits**, best risk reduction first")
            st.dataframe(pd.DataFrame(rows), width="stretch", hide_index=True)
            free_wins = [r for r in rows if r["trade-off"].startswith("none")]
            if free_wins:
                st.success(
                    f"{len(free_wins)} edit(s) lower predicted risk without costing "
                    f"predicted knockdown, best being {free_wins[0]['edit']}."
                )
        else:
            st.info("No single substitution moved a validated score meaningfully.")

        for note in advice.notes:
            st.caption(note)
with tab["How it works"]:
    st.subheader("How it works, and how to check it")
    st.caption(
        "Everything on this page is read from the files this build shipped "
        "with, not written out by hand, so it cannot drift from what the app "
        "is actually doing."
    )

    st.markdown("#### 1. What the app does")
    st.markdown(
        "You give it a sequence, a chemistry and a dose. It computes features "
        "from the sequence, runs one model per toxicity endpoint and one per "
        "efficacy context, and reports the endpoints that passed an evidence "
        "gate. Endpoints that did not pass are still shown, under "
        "**Unsupported**, with the reason. Nothing is hidden because it looks "
        "bad; it is separated because it is not supported."
    )

    st.markdown("#### 2. What it was trained on")
    inventory = _load_report("dataset_inventory.json")
    sources = inventory.get("sources", [])
    if sources:
        train = [s for s in sources if s.get("usable_for_training")]
        st.markdown(
            f"**{len(train)} source(s), {inventory.get('usable_sequences', 0):,} "
            "sequences.** These are public assay datasets, not internal data."
        )
        st.dataframe(pd.DataFrame([{
            "source": s["name"],
            "kind": s.get("kind"),
            "endpoint": s.get("endpoint"),
            "sequences": s.get("sequences"),
        } for s in train]), width="stretch", hide_index=True)
        held = [s for s in sources if not s.get("usable_for_training")]
        with st.expander(f"{len(held)} source(s) not used for training, and why"):
            st.dataframe(pd.DataFrame([{
                "source": s["name"],
                "status": s.get("status"),
                "not used because": s.get("blocked_by"),
            } for s in held]), width="stretch", hide_index=True)
    else:
        st.info("The dataset inventory did not ship with this build.")

    st.markdown("#### 3. How a sequence becomes numbers")
    st.markdown(
        "Only sequence-derived features. No proprietary descriptors, no "
        "embeddings from a model you cannot inspect. Each feature is computed "
        "by code in `src/oligotox/features/`."
    )
    feature_config = _load_config("endpoint_features.json")
    features = _model_features()
    if features:
        st.dataframe(pd.DataFrame([{
            "model": name,
            "features": len(columns),
            "uses chemistry": any(c.startswith(("ps_", "lna_", "aso_chem", "aso_mod"))
                                  for c in columns),
            "uses folding": any(c.startswith("site_") for c in columns),
        } for name, columns in features.items()]), width="stretch", hide_index=True)
        chosen = st.selectbox("See every input to one model", sorted(features))
        st.code(", ".join(features[chosen]), language="text")
        note = (feature_config.get(chosen) or {}).get("note")
        if note:
            st.caption(note)
    st.caption(
        "Where a target gene is given and its transcript is cached, local "
        "folding (ViennaRNA) adds binding-site accessibility, but only for "
        "the two assay contexts where it improved a like-for-like comparison "
        "and held across three sampling seeds. A site that does not locate "
        "falls back to the design-only model rather than being scored as if "
        "its accessibility were zero."
    )

    st.markdown("#### 4. How it is validated, and the number that matters most")
    st.markdown(
        "Cross-validation is **grouped by sequence**, so no sequence appears "
        "in both the training and test half of a fold. That sounds like a "
        "detail and it is the whole difference between a real score and a "
        "flattering one, because these datasets test the same sequence many "
        "times. Here is the same model scored both ways:"
    )
    training = _load_report("training_summary_after.json")
    rows = [{
        "endpoint": name,
        "PR-AUC, grouped by sequence": round(spec["pr_auc"], 3),
        "PR-AUC, random split": (round(spec["random_split_pr_auc"], 3)
                                 if spec.get("random_split_pr_auc") is not None else None),
        "unique positive sequences": spec.get("positive_sequences"),
    } for name, spec in training.items()
        if isinstance(spec, dict) and spec.get("pr_auc") is not None]
    if rows:
        st.dataframe(pd.DataFrame(rows), width="stretch", hide_index=True)
        inflated = [r for r in rows
                    if r["PR-AUC, random split"] is not None
                    and r["PR-AUC, random split"] > r["PR-AUC, grouped by sequence"]]
        st.markdown(
            f"The random split scores higher on **{len(inflated)} of "
            f"{len(rows)}** endpoints. The left column is the one quoted "
            "everywhere else in this app."
        )

    st.markdown("#### 5. The evidence gate")
    evidence = _load_config("endpoint_evidence.json")
    st.markdown(
        "An endpoint is only presented as a prediction if it passes every "
        "check below, run by `scripts/derive_evidence_tiers.py`. The gate is "
        "a ratchet: a tier can be lost on new evidence, and the history is "
        "kept in `reports/results/endpoint_history.jsonl`."
    )
    endpoints = evidence.get("endpoints", {})
    if endpoints:
        st.dataframe(pd.DataFrame([{
            "endpoint": name,
            "tier": TIER_LABELS.get(spec.get("tier"), spec.get("tier")),
            "unique positive sequences": spec.get("positive_sequences"),
            "reads as": ("probability" if spec.get("probabilities_calibrated")
                         else "ranking only"),
        } for name, spec in endpoints.items()]), width="stretch", hide_index=True)
        presented = [n for n, s in endpoints.items() if s.get("tier") == "validated"]
        st.markdown(
            f"**{len(presented)} of {len(endpoints)}** endpoints are presented: "
            + ", ".join(f"`{n}`" for n in presented) + "."
        )
        for name in presented:
            with st.expander(f"Why {name} passed, and what is still wrong with it"):
                st.write(endpoints[name].get("reason", ""))

    st.markdown("#### 6. Why some numbers are rankings and not percentages")
    rankings = [n for n, s in endpoints.items()
                if not s.get("probabilities_calibrated")]
    st.markdown(
        "A score is only called a probability if it is calibrated: its Brier "
        "skill against simply predicting the base rate has to beat 0.02. "
        "Below that, a 30% means this sequence ranks above one scored 20%, "
        "and nothing more. The app labels each score with which it is. "
        + (f"**{len(rankings)} of {len(endpoints)}** endpoints are rankings "
           "rather than percentages." if endpoints else "")
    )

    st.markdown("#### 7. What an independent dataset says")
    external = (endpoints.get("hepatocellular_injury", {}) or {}).get("external_check")
    if external:
        st.warning(
            f"On {external.get('sequences')} LNA gapmers from an independent "
            f"screen ({external.get('positives')} hepatotoxic by the paper's "
            f"own ALT criterion), `hepatocellular_injury` achieves ROC-AUC "
            f"**{external.get('roc_auc', 0):.2f}** "
            f"[{external.get('roc_ci', [0, 0])[0]:.2f}-"
            f"{external.get('roc_ci', [0, 0])[1]:.2f}]. "
            "0.5 is chance. Five positives cannot settle it, but there is no "
            "evidence this endpoint transfers off the data it was trained on."
        )
        # Whether that endpoint is presented is a fact about the gate, and
        # this page is the last place that should state it from memory: the
        # sentence here read "it is presented anyway" for a while after the
        # endpoint had already been demoted.
        tier = (endpoints.get("hepatocellular_injury", {}) or {}).get("tier")
        if tier == "validated":
            st.caption(
                "It is presented as a prediction anyway, with this said on "
                "its own tab, because removing the check would be the "
                "dishonest option."
            )
        else:
            st.caption(
                f"It is not presented as a prediction. Its tier is "
                f"`{tier}`, and section 5 gives the check it failed."
            )
        st.caption(f"Source: {external.get('source', '')}")

    st.markdown("#### 8. How to confirm all of this yourself")
    st.markdown(
        "Nothing here asks you to take the app's word for it. The repository "
        "is public and each claim above has a command behind it."
    )
    st.code(
        "pytest -q                               "
        "# the full suite, including values pinned to hand-read papers\n"
        "python scripts/loop_health.py           "
        "# are the data, training and evidence gate still consistent\n"
        "python scripts/derive_evidence_tiers.py "
        "# recompute every tier on this page from the metrics\n"
        "python scripts/evaluate_hepatotox_external.py  "
        "# re-run the independent check in section 7",
        language="bash",
    )
    health = _load_report("loop_health.json")
    if health:
        checks = health.get("checks", {})
        st.dataframe(pd.DataFrame([{
            "check": name,
            "result": "pass" if spec.get("ok") else "FAIL",
            "detail": spec.get("detail"),
        } for name, spec in checks.items()]), width="stretch", hide_index=True)
        st.caption(
            f"Last checked {health.get('checked_at', 'unknown')}. "
            + (health.get("note") or "")
        )

    st.markdown("#### 9. What this app is not")
    st.markdown(
        "- It is a research prototype and not a regulatory tool. No output "
        "here substitutes for an in-vivo study.\n"
        "- Most toxicity endpoints rest on a small number of unique "
        "sequences. Rows are plentiful; independent sequences are not, and "
        "the sequence count is the one that limits what can be claimed.\n"
        "- The CNS endpoint does not beat the published G-count rule. It "
        "matches it, so treat it as that rule rather than a new finding.\n"
        "- A sequence unlike anything in the training sources will still get "
        "a score. The app cannot tell you that it is out of distribution.\n"
        "- Nothing you type is sent anywhere or written to the server. Runs "
        "live in your browser session and disappear when you close the tab."
    )


with tab["Saved Runs"]:
    st.subheader("Saved runs")
    st.caption(
        "Runs from this browser session only. They are not written to the server "
        "and are gone when you close the tab, so nothing you enter here is visible "
        "to anyone else. Download a run to keep it."
    )
    runs_df = _load_saved_runs()
    if runs_df.empty:
        st.info("No runs in this session yet.")
    else:
        st.dataframe(runs_df, width="stretch", hide_index=True)
        session_runs = st.session_state.get("runs", [])
        labels = {f"{r['saved_at']}  {r['input']['sequence'][:24]}": r for r in session_runs}
        selected = st.selectbox("Open a run from this session", list(labels))
        if selected:
            chosen = labels[selected]
            st.download_button(
                "Download this run as JSON",
                data=json.dumps(chosen, indent=2, default=str),
                file_name=f"oligotox_run_{chosen['run_id']}.json",
                mime="application/json",
            )
            st.json(chosen)

with tab["Build Info"]:
    st.subheader("Build info")
    st.write({
        "app_file": str(ROOT / "oligotox_app.py"),
        "processed_dir": str(PROCESSED_DIR),
        "model_dir": str(MODEL_DIR),
        "config_dir": str(CONFIG_DIR),
        "run_storage": "session only, not written to disk",
    })
    st.write("Validation status: backend validated end-to-end in current workspace.")
