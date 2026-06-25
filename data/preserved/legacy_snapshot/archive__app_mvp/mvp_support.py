from __future__ import annotations
from dataclasses import dataclass
from io import BytesIO
from pathlib import Path
from typing import Dict, List, Any
import json
import math
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

ENDPOINT_LABELS = [
    "hepatocellular_injury",
    "hepatic_vacuolation",
    "renal_tubular_injury",
    "renal_accumulation",
    "thrombocytopenia",
    "coagulation_effect",
    "complement_activation",
    "innate_immune_activation",
    "cytokine_release",
    "kupffer_cell_activation",
    "cns_toxicity",
    "offtarget_gene_silencing",
    "injection_site_reaction",
]

ENDPOINT_TO_BUCKET = {
    "hepatocellular_injury": "D2",
    "hepatic_vacuolation": "D2",
    "renal_tubular_injury": "D2",
    "renal_accumulation": "D3",
    "thrombocytopenia": "D3",
    "coagulation_effect": "D3",
    "complement_activation": "D3",
    "innate_immune_activation": "D2",
    "cytokine_release": "D2",
    "kupffer_cell_activation": "D2",
    "cns_toxicity": "D1",
    "offtarget_gene_silencing": "D1",
    "injection_site_reaction": "D3",
}

ASSAY_LIBRARY = [
    {
        "assay_name": "Serum chemistry panel",
        "tier": "Tier 1",
        "cost_tier": "$",
        "triggers": ["hepatocellular_injury", "hepatic_vacuolation", "renal_tubular_injury"],
        "rationale": "Screens early liver and kidney liability using routine biomarkers.",
    },
    {
        "assay_name": "Urinalysis + cystatin C",
        "tier": "Tier 1",
        "cost_tier": "$",
        "triggers": ["renal_tubular_injury", "renal_accumulation"],
        "rationale": "Useful for proximal tubular injury and renal monitoring signals.",
    },
    {
        "assay_name": "Platelet count + coagulation panel",
        "tier": "Tier 1",
        "cost_tier": "$",
        "triggers": ["thrombocytopenia", "coagulation_effect"],
        "rationale": "Captures platelet and coagulation liabilities seen in PS-linked tox.",
    },
    {
        "assay_name": "Complement split products / cytokine panel",
        "tier": "Tier 2",
        "cost_tier": "$$",
        "triggers": ["complement_activation", "innate_immune_activation", "cytokine_release"],
        "rationale": "Profiles infusion-related and innate immune activation risk.",
    },
    {
        "assay_name": "Kupffer / liver macrophage histopathology",
        "tier": "Tier 2",
        "cost_tier": "$$",
        "triggers": ["kupffer_cell_activation", "hepatocellular_injury"],
        "rationale": "Adds mechanistic resolution when liver-directed innate signals appear.",
    },
    {
        "assay_name": "Off-target transcriptomics / seed-match validation",
        "tier": "Tier 2",
        "cost_tier": "$$",
        "triggers": ["offtarget_gene_silencing"],
        "rationale": "Validates sequence-driven off-target burden and mechanistic plausibility.",
    },
    {
        "assay_name": "Repeat-dose tox package",
        "tier": "Tier 3",
        "cost_tier": "$$$",
        "triggers": ENDPOINT_LABELS,
        "rationale": "Escalation package for candidates with persistent or multi-system risk.",
    },
]

BIOMARKER_LIBRARY = {
    "liver": ["ALT", "AST", "ALP", "bilirubin"],
    "renal": ["serum creatinine", "cystatin C", "urine protein", "urine dipstick"],
    "platelet": ["platelet count", "immature platelet fraction"],
    "immune": ["C3a", "C5a", "IL-6", "TNF-alpha", "CRP"],
}

def confidence_badge(max_risk: float) -> str:
    if max_risk >= 0.66:
        return "High"
    if max_risk >= 0.33:
        return "Medium"
    return "Low"

def normalize_scores(scores: Dict[str, float]) -> Dict[str, float]:
    out = {}
    for k, v in scores.items():
        try:
            out[k] = max(0.0, min(1.0, float(v)))
        except Exception:
            out[k] = 0.0
    return out

def mechanism_bucket_scores(scores: Dict[str, float]) -> Dict[str, float]:
    scores = normalize_scores(scores)
    buckets = {"D1": 0.0, "D2": 0.0, "D3": 0.0}
    for endpoint, score in scores.items():
        bucket = ENDPOINT_TO_BUCKET.get(endpoint)
        if bucket:
            buckets[bucket] = max(buckets[bucket], score)
    return buckets

def top_endpoint(scores: Dict[str, float]) -> Dict[str, Any]:
    scores = normalize_scores(scores)
    if not scores:
        return {"endpoint": None, "score": 0.0}
    endpoint = max(scores, key=scores.get)
    return {"endpoint": endpoint, "score": scores[endpoint]}

def build_assay_plan(candidate: Dict[str, Any], scores: Dict[str, float], stage: str = "pre-IND") -> pd.DataFrame:
    scores = normalize_scores(scores)
    active = {k for k, v in scores.items() if v >= 0.25}
    stage = str(stage).strip()

    rows = []
    for assay in ASSAY_LIBRARY:
        overlap = sorted(active.intersection(set(assay["triggers"])))
        include = len(overlap) > 0

        if stage.upper() == "CN" and assay["tier"] == "Tier 3":
            include = False
        if stage.upper() == "GLP" and assay["tier"] in {"Tier 1", "Tier 2", "Tier 3"} and len(active) > 0:
            include = True

        if include:
            rows.append({
                "assay_name": assay["assay_name"],
                "tier": assay["tier"],
                "cost_tier": assay["cost_tier"],
                "triggered_by": ", ".join(overlap) if overlap else "stage-based escalation",
                "rationale": assay["rationale"],
                "stage": stage,
                "modality": candidate.get("modality", ""),
            })

    if not rows:
        rows.append({
            "assay_name": "Baseline serum chemistry + urinalysis",
            "tier": "Tier 1",
            "cost_tier": "$",
            "triggered_by": "baseline",
            "rationale": "Default low-risk screen for early candidate profiling.",
            "stage": stage,
            "modality": candidate.get("modality", ""),
        })

    return pd.DataFrame(rows)

def simple_seed_matches(sequence: str, target: str) -> pd.DataFrame:
    seq = (sequence or "").upper().replace("T", "U")
    target = (target or "").upper().replace("T", "U")
    seed = seq[1:8] if len(seq) >= 8 else seq
    hits = []

    if not seed or not target:
        return pd.DataFrame(columns=["seed", "position", "window", "match_type", "score"])

    for i in range(0, max(0, len(target) - len(seed) + 1)):
        window = target[i:i + len(seed)]
        matches = sum(1 for a, b in zip(seed, window) if a == b)
        score = matches / len(seed) if len(seed) else 0.0
        if score >= 0.57:
            hits.append({
                "seed": seed,
                "position": i + 1,
                "window": window,
                "match_type": "seed_like" if score < 1.0 else "exact_seed",
                "score": round(score, 3),
            })
    return pd.DataFrame(hits).sort_values(["score", "position"], ascending=[False, True])

def summarize_offtarget(sequence: str, target: str) -> Dict[str, Any]:
    hits = simple_seed_matches(sequence, target)
    if hits.empty:
        return {
            "summary_risk": 0.05,
            "top_hit_count": 0,
            "best_seed_score": 0.0,
            "hits": hits,
        }
    best = float(hits["score"].max())
    count = int((hits["score"] >= 0.71).sum())
    summary = min(0.95, 0.15 + best * 0.5 + min(count, 10) * 0.03)
    return {
        "summary_risk": round(summary, 3),
        "top_hit_count": count,
        "best_seed_score": round(best, 3),
        "hits": hits,
    }

def build_translational_bridge(preclinical: Dict[str, Any]) -> pd.DataFrame:
    species = str(preclinical.get("species", "rat")).lower()
    alt_fold = float(preclinical.get("alt_fold", 1.0) or 1.0)
    platelet_nadir = float(preclinical.get("platelet_nadir", 200.0) or 200.0)
    creatinine_fold = float(preclinical.get("creatinine_fold", 1.0) or 1.0)

    recommendations = []

    if alt_fold >= 2.0:
        recommendations.append({
            "signal_domain": "liver",
            "preclinical_signal": f"ALT fold={alt_fold}",
            "clinical_panel": ", ".join(BIOMARKER_LIBRARY["liver"]),
            "species_correction": "Rat liver signals may overcall mild human risk; interpret with exposure context.",
        })

    if creatinine_fold >= 1.3:
        correction = "Rat renal findings can overestimate human renal liability for some oligos."
        if species == "monkey":
            correction = "Cynomolgus renal findings may be more translational for some delivery systems."
        recommendations.append({
            "signal_domain": "renal",
            "preclinical_signal": f"Creatinine fold={creatinine_fold}",
            "clinical_panel": ", ".join(BIOMARKER_LIBRARY["renal"]),
            "species_correction": correction,
        })

    if platelet_nadir < 150:
        recommendations.append({
            "signal_domain": "platelet",
            "preclinical_signal": f"Platelet nadir={platelet_nadir}",
            "clinical_panel": ", ".join(BIOMARKER_LIBRARY["platelet"]),
            "species_correction": "Monitor for exaggerated PS-linked platelet effects and compare across species.",
        })

    if not recommendations:
        recommendations.append({
            "signal_domain": "general",
            "preclinical_signal": "No major bridge trigger",
            "clinical_panel": "ALT, AST, creatinine, cystatin C, platelets, CRP",
            "species_correction": "Use routine biomarker panel and exposure-normalized interpretation.",
        })

    return pd.DataFrame(recommendations)

def candidate_brief_html(candidate: Dict[str, Any], scores: Dict[str, float]) -> str:
    scores = normalize_scores(scores)
    top = top_endpoint(scores)
    buckets = mechanism_bucket_scores(scores)
    badge = confidence_badge(top["score"])

    score_rows = "".join(
        f"<tr><td>{k}</td><td>{v:.3f}</td><td>{ENDPOINT_TO_BUCKET.get(k,'')}</td></tr>"
        for k, v in sorted(scores.items(), key=lambda kv: kv[1], reverse=True)
    )

    payload = json.dumps(candidate, indent=2)

    return f"""
    <html>
    <head>
      <meta charset="utf-8">
      <style>
        body {{ font-family: Arial, sans-serif; margin: 28px; color: #1f2937; }}
        h1,h2 {{ margin-bottom: 8px; }}
        .badge {{ display:inline-block; padding:6px 10px; border-radius:999px; background:#e5e7eb; }}
        .grid {{ display:grid; grid-template-columns: 1fr 1fr; gap: 16px; }}
        .card {{ border:1px solid #d1d5db; border-radius:12px; padding:14px; margin-bottom:16px; }}
        table {{ border-collapse: collapse; width:100%; }}
        th,td {{ border:1px solid #d1d5db; padding:8px; text-align:left; font-size: 12px; }}
        pre {{ background:#f3f4f6; padding:12px; border-radius:8px; overflow:auto; }}
      </style>
    </head>
    <body>
      <h1>OligoTox Candidate Safety Brief</h1>
      <div class="card">
        <strong>Top endpoint:</strong> {top["endpoint"]} ({top["score"]:.3f})<br/>
        <strong>Confidence:</strong> <span class="badge">{badge}</span><br/>
        <strong>D1:</strong> {buckets["D1"]:.3f} |
        <strong>D2:</strong> {buckets["D2"]:.3f} |
        <strong>D3:</strong> {buckets["D3"]:.3f}
      </div>

      <div class="grid">
        <div class="card">
          <h2>Candidate</h2>
          <pre>{payload}</pre>
        </div>
        <div class="card">
          <h2>Interpretation</h2>
          <p>This brief is intended for exploratory safety triage and mechanism-aware candidate comparison.</p>
          <p>Highest current liability is <strong>{top["endpoint"]}</strong>, with mechanism signal distributed across D1/D2/D3 as above.</p>
        </div>
      </div>

      <div class="card">
        <h2>Endpoint scores</h2>
        <table>
          <thead><tr><th>Endpoint</th><th>Score</th><th>Bucket</th></tr></thead>
          <tbody>{score_rows}</tbody>
        </table>
      </div>
    </body>
    </html>
    """

def html_to_pdf_bytes(html: str) -> bytes:
    try:
        from weasyprint import HTML
        return HTML(string=html).write_pdf()
    except Exception:
        return html.encode("utf-8")

def assay_plan_xlsx_bytes(df: pd.DataFrame) -> bytes:
    output = BytesIO()
    with pd.ExcelWriter(output, engine="openpyxl") as writer:
        df.to_excel(writer, index=False, sheet_name="Assay Plan")
        ws = writer.sheets["Assay Plan"]
        for col_cells in ws.columns:
            max_len = 0
            col_letter = col_cells[0].column_letter
            for cell in col_cells:
                val = "" if cell.value is None else str(cell.value)
                max_len = max(max_len, len(val))
            ws.column_dimensions[col_letter].width = min(max_len + 2, 48)
    return output.getvalue()
