"""smoke test for run_oligotox inference API"""
from oligotox.inference.api import run_oligotox

def test_run_oligotox_smoke():
    inp = {
        "sequence": "AUGCUACGUACGUACGUACG",
        "modality": "siRNA_LNP",
        "target_gene": "TTR",
        "is_lnp": True,
        "is_galnac": False,
        "is_iv": True,
        "is_it": False,
        "dose_mg_kg": 0.3,
        "repeat_dose": True,
        "ps_fraction": 0.7,
        "lna_fraction": 0.0,
        "particle_size": 80.0,
        "pdi": 0.15,
        "peg_pct": 1.5,
        "c3a": 0.1,
        "c5a": 0.05,
        "sc5b9": 0.2,
        "alt_fold": 1.0,
    }
    result = run_oligotox(inp)
    assert isinstance(result, dict), "run_oligotox should return a dict"
    assert "per_endpoint" in result or "composite_score" in result, f"Missing expected keys: {result.keys()}"
