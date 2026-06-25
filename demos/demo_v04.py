"""
OligoTox v0.4 — Demo: three candidates with full module output.
"""
from v04.api import run_oligotox_v04

CANDIDATES = [
    {
        "name": "OTX-001 — IV LNP-siRNA (PCSK9 liver)",
        "sequence": "GCAUCAAGUGCCUACAGUA",
        "features": {
            "ps_fraction": 0.40, "lna_fraction": 0.05, "gc_content": 0.52,
            "seed_match_score": 0.58, "cpg_density": 0.02, "tlr7_score": 0.30,
            "tlr8_score": 0.25, "tlr9_score": 0.15, "self_complementarity": 0.18,
            "is_lnp": 1.0, "is_iv": 1.0, "is_galnacc": 0.0, "is_intrathecal": 0.0,
            "dose_mg_kg": 0.3, "repeat_dose": 1.0, "conservation_score": 0.87,
            "is_aso": 0.0,
        },
        "n_of_1_config": {
            "target_gene": "PCSK9",
            "disease_context": "hypercholesterolemia",
            "timeline_weeks": 16,
        },
        "translational_config": {
            "features": {
                "ps_fraction": 0.40, "is_lnp": 1.0, "is_iv": 1.0,
                "tlr7_score": 0.30, "tlr8_score": 0.25, "complement_activation_preclinical": 0.45,
                "conservation_score": 0.87,
            },
            "species": "rat+NHP",
            "disease_context": "metabolic_liver",
            "has_surrogate": True,
            "dart_relevant": False,
        },
    },
    {
        "name": "OTX-002 — SC GalNAc-siRNA (TTR liver, rare disease)",
        "sequence": "CUGUUCAAGGAGACCAUUC",
        "features": {
            "ps_fraction": 0.25, "lna_fraction": 0.0, "gc_content": 0.47,
            "seed_match_score": 0.32, "cpg_density": 0.005, "tlr7_score": 0.12,
            "tlr8_score": 0.10, "tlr9_score": 0.05, "self_complementarity": 0.10,
            "is_lnp": 0.0, "is_iv": 0.0, "is_galnacc": 1.0, "is_intrathecal": 0.0,
            "dose_mg_kg": 0.1, "repeat_dose": 0.0, "conservation_score": 0.91,
            "is_aso": 0.0,
        },
        "n_of_1_config": {
            "target_gene": "TTR",
            "disease_context": "rare hereditary transthyretin amyloidosis",
            "timeline_weeks": 12,
        },
        "translational_config": {
            "features": {
                "ps_fraction": 0.25, "is_lnp": 0.0, "is_galnacc": 1.0,
                "conservation_score": 0.91,
            },
            "species": "rat+NHP",
            "disease_context": "rare_hereditary_amyloidosis",
            "has_surrogate": True,
            "dart_relevant": False,
        },
    },
    {
        "name": "OTX-003 — IT ASO-gapmer (HTT CNS, rare disease n-of-1)",
        "sequence": "CTGTCCGTGTGAGCGAAGG",
        "features": {
            "ps_fraction": 0.80, "lna_fraction": 0.30, "gc_content": 0.65,
            "seed_match_score": 0.70, "cpg_density": 0.055, "tlr7_score": 0.35,
            "tlr8_score": 0.28, "tlr9_score": 0.45, "self_complementarity": 0.38,
            "is_lnp": 0.0, "is_iv": 0.0, "is_galnacc": 0.0, "is_intrathecal": 1.0,
            "dose_mg_kg": 2.0, "repeat_dose": 1.0, "conservation_score": 0.55,
            "is_aso": 1.0,
        },
        "n_of_1_config": {
            "target_gene": "HTT",
            "disease_context": "rare CNS huntington disease n-of-1",
            "timeline_weeks": 20,
        },
        "translational_config": {
            "features": {
                "ps_fraction": 0.80, "is_intrathecal": 1.0, "repeat_dose": 1.0,
                "tlr9_score": 0.45, "conservation_score": 0.55,
            },
            "species": "rat+NHP",
            "disease_context": "rare CNS huntington n-of-1",
            "has_surrogate": False,
            "dart_relevant": True,
        },
    },
]


def section(title: str):
    print("\n" + "=" * 72)
    print(f"  {title}")
    print("=" * 72)


def run_demo():
    for cand in CANDIDATES:
        section(cand["name"])
        result = run_oligotox_v04(
            sequence=cand["sequence"],
            features=cand["features"],
            n_of_1_config=cand.get("n_of_1_config"),
            translational_config=cand.get("translational_config"),
        )

        ov = result["overall_v04_risk"]
        print(f"\nOVERALL RISK: {ov['composite_score']:.3f}  [{ov['band']}]")
        for k, v in ov["drivers"].items():
            print(f"  {k:<32} {v:.3f}")

        qa = result["qAOP"]
        print(f"\nqAOP TOP: {qa['top_aop']}  (P={qa['top_ao_probability']:.3f})")
        print(f"  {qa['narratives'][qa['top_aop']][:120]}...")

        rd = result["redesign"]
        if rd["top_recommendations"]:
            top = rd["top_recommendations"][0]
            print(f"\nBEST REDESIGN: {top['rule_name']}")
            print(f"  composite delta = {top['composite_delta']:.3f}  |  tradeoff: {top['tradeoff'][:80]}")

        imp = result["impurity_tox"]
        print(f"\nIMPURITY: {imp['total_variants_assessed']} variants assessed, {imp['high_risk_variants']} high-risk, index={imp['impurity_liability_index']}")
        for v in imp["top_flagged"][:2]:
            print(f"  [{v['type']}] pos={v['position']} delta={v['delta_risk_score']:.3f} flags={v['flags']}")

        n1 = result["n_of_1"]["human_offtarget_screen"]
        print(f"\nN-OF-1 OT SCREEN: {n1['total_hits']} hits, {n1['critical_hits']} CRITICAL")
        if n1["critical_genes"]:
            print(f"  Critical genes: {', '.join(n1['critical_genes'][:5])}")
        if result["n_of_1"]["minimum_safety_package"]:
            pkg = result["n_of_1"]["minimum_safety_package"]
            print(f"  ICH-S12 applies: {pkg['ics12_applies']}")
            print(f"  Min assay package ({len(pkg['minimum_assay_package'])} assays):")
            for a in pkg["minimum_assay_package"][:3]:
                print(f"    P{a['priority']} | {a['assay'][:70]}")

        tf = result["translational_failure"]
        if "translational_failure_risk" in tf:
            print(f"\nTRANSLATIONAL FAILURE RISK: {tf['translational_failure_risk']:.3f}  [{tf['failure_risk_level']}]")
            print(f"  Most likely failure mode: {tf['most_likely_failure_mode'][:80]}")
            if tf["clinical_monitoring_triggers"]:
                print("  Clinical monitoring triggers:")
                for t in tf["clinical_monitoring_triggers"][:3]:
                    print(f"    - {t[:90]}")

    print("\n" + "=" * 72)
    print("  OligoTox v0.4 demo complete.")
    print("=" * 72)


if __name__ == "__main__":
    run_demo()
