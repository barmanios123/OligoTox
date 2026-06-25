import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from api import run_oligotox_unified


class Candidate:
    sequence = "GCAUCAAGUGCCUACAGUA"
    is_lnp = True
    is_galnac = False
    is_iv = True
    is_it = False
    dose_mg_kg = 0.3
    repeat_dose = True
    ps_fraction = 0.4
    lna_fraction = 0.05
    particle_size = 80
    pdi = 0.2
    c3a = 0.2
    c5a = 0.15
    sc5b9 = 0.25
    alt_fold = 1.0
    ast_fold = 1.0
    bun_fold = 1.0
    cr_fold = 1.0
    il6 = 0.1
    tnfa = 0.1
    ifna = 0.1
    pw_nfkb = 0.2
    pw_ifn = 0.2
    pw_apop = 0.1
    pw_comp = 0.2
    pw_plt = 0.3
    pw_kyn = 0.1
    target_gene = "PCSK9"
    disease_context = "hypercholesterolemia"
    timeline_weeks = 16
    species = "rat+NHP"
    has_surrogate = True
    dart_relevant = False


def test_unified_api_smoke():
    result = run_oligotox_unified(Candidate())
    assert "v03_model" in result
    assert "v04_decision_support" in result
    assert isinstance(result["v04_decision_support"], dict)  # v04 deprecated, returns empty dict
    print(result["v04_decision_support"])  # v04 deprecated


if __name__ == "__main__":
    test_unified_api_smoke()
