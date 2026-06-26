# tests/smoke_test_api.py

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from oligotox.inference.api import run_oligotox


def test_smoke_candidate():
    candidate = {
        "sequence": "GCUAUGCUAGCUAGCUAGCUAGCUAGC",
        "is_lnp": True,
        "is_galnac": False,
        "is_iv": True,
        "is_it": False,
        "dose_mg_kg": 1.5,
        "repeat_dose": False,
        "ps_fraction": 0.7,
        "lna_fraction": 0.2,
        "particle_size": 85,
        "pdi": 0.18,
        "peg_pct": 1.5,
        "c3a": 0.2,
        "c5a": 0.1,
        "sc5b9": 0.25,
        "alt_fold": 1.1,
        "ast_fold": 1.0,
        "bun_fold": 1.0,
        "cr_fold": 1.0,
        "il6": 0.1,
        "tnfa": 0.05,
        "ifna": 0.05,
        "pw_nfkb": 0.2,
        "pw_ifn": 0.1,
        "pw_apop": 0.1,
        "pw_comp": 0.25,
        "pw_plt": 0.15,
        "pw_kyn": 0.05,
        "micropeptide_risk": 0.1,
        "human_specific": True,
    }

    result = run_oligotox(candidate)

    assert "composite_score" in result
    assert "composite_band" in result
    assert "dominant_endpoint" in result
    assert "per_endpoint" in result
    assert isinstance(result["composite_score"], float)
    assert isinstance(result["per_endpoint"], dict)
    assert len(result["per_endpoint"]) > 0

    print("=" * 80)
    print("SMOKE TEST: OligoTox scoring")
    print("=" * 80)
    print(f"Composite score : {result['composite_score']:.3f}")
    print(f"Composite band  : {result['composite_band']}")
    print(f"Dominant endpoint: {result['dominant_endpoint']}")
    print()

    for endpoint, info in sorted(result["per_endpoint"].items()):
        print(
            f"{endpoint:28s} "
            f"prob={info['prob']:.3f} "
            f"threshold={info['threshold']} "
            f"band={info['risk_band']:8s} "
            f"positive={info['predicted_positive']}"
        )


if __name__ == "__main__":
    test_smoke_candidate()