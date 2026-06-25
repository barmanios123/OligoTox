import pytest
pytest.importorskip("app.legacy_mvp.mvp_support", reason="legacy MVP archived")
from app.legacy_mvp.mvp_support import (
    build_assay_plan,
    summarize_offtarget,
    build_translational_bridge,
    candidate_brief_html,
    mechanism_bucket_scores,
)

def test_assay_plan():
    candidate = {"modality": "siRNA"}
    scores = {"offtarget_gene_silencing": 0.51, "renal_tubular_injury": 0.30}
    df = build_assay_plan(candidate, scores, stage="pre-IND")
    assert len(df) >= 2
    assert "assay_name" in df.columns

def test_offtarget():
    out = summarize_offtarget("AUGGCUACGUAGCUAGCUAG", "CCAUGGCUACGUAGCUAGCUAGGGA")
    assert "summary_risk" in out
    assert out["summary_risk"] >= 0.0

def test_translation():
    df = build_translational_bridge({"species": "rat", "alt_fold": 2.5, "platelet_nadir": 120})
    assert len(df) >= 1

def test_html():
    html = candidate_brief_html({"sequence": "AUGG"}, {"offtarget_gene_silencing": 0.5})
    assert "OligoTox Candidate Safety Brief" in html

def test_buckets():
    buckets = mechanism_bucket_scores({"offtarget_gene_silencing": 0.5, "renal_accumulation": 0.2})
    assert buckets["D1"] == 0.5
