"""Three-Domain Mechanistic Rules Engine v1"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from typing import Dict, Tuple
from dataclasses import dataclass
from oligotox.schemas.data_models import OligoToxInput, Modality, DeliveryVehicle
from oligotox.schemas.tox_endpoints import ToxEndpoint, ToxPrediction, RiskClass, RiskDomain, ENDPOINT_DOMAIN

@dataclass
class RuleResult:
    score: float; confidence: float; features: list; note: str

def classify_risk(s: float) -> RiskClass:
    return RiskClass.HIGH if s >= 0.60 else (RiskClass.MODERATE if s >= 0.30 else RiskClass.LOW)

def score_hepatocellular_injury(f: dict) -> RuleResult:
    s = 0.0
    s += min(f['lna_5flank'] * 1.8, 0.55)  # Dieckmann 2017: LNA 5'-flank caspase-3
    s += min(f['ps_fraction'] * 0.30, 0.25)
    s += min(f['alt_fold'] * 0.20 + f['ast_fold'] * 0.15, 0.40)
    s += f['pw_apop'] * 0.25
    note = 'LNA 5-flank caspase-3 (Dieckmann 2017); PS class hepatotox; ALT/AST biomarkers.'
    return RuleResult(min(s,1.0), 0.80, ['lna_5flank','ps_fraction','alt_fold'], note)

def score_renal_tubular_injury(f: dict) -> RuleResult:
    s = min(f['ps_fraction'] * 0.55, 0.50)
    s += min(f['dose_mg_kg'] / 200, 0.25)
    s += min(f['bun_fold'] * 0.20 + f['cr_fold'] * 0.20, 0.30)
    s += f['repeat_dose'] * 0.10
    note = 'PS-ASO renal accumulation; KIM-1/NGAL biomarkers. Rat overestimates human (ICH S7A).'
    return RuleResult(min(s,1.0), 0.75, ['ps_fraction','dose_mg_kg','bun_fold'], note)

def score_thrombocytopenia(f: dict) -> RuleResult:
    s = min(f['ps_fraction'] * 0.60, 0.55)
    s += min(f['platelet_drop_fraction'] * 0.60, 0.45)
    s += min(f['pf4_od'] * 0.20, 0.25)
    s += f['c3d_positive'] * 0.15
    note = 'PF4-PS complex mechanism; platelet-bound C3d/C4d (CARPA). Inotersen black box ref.'
    return RuleResult(min(s,1.0), 0.78, ['ps_fraction','platelet_drop_fraction','pf4_od'], note)

def score_complement_activation(f: dict) -> RuleResult:
    s = 0.0
    s += f['is_lnp'] * 0.65  # LNP CARPA (Patisiran NDA)
    s += min((f['c3a']/500 + f['c5a']/300 + f['sc5b9']/30) * 0.35, 0.45)
    s += f['carpa'] * 0.25
    s += f['pw_comp'] * 0.20
    s += f['c3d_positive'] * 0.10
    note = 'LNP CARPA (Patisiran IND); C3a/C5a/SC5b-9 analytes; CH50/AH50 ex vivo.'
    return RuleResult(min(s,1.0), 0.79, ['is_lnp','c3a','c5a','carpa'], note)

def score_innate_immune(f: dict) -> RuleResult:
    s = min(f['cpg_density'] * 3.0, 0.40)  # TLR9
    s += min(f['poly_u_density'] * 4.0, 0.30)  # TLR7/8
    s += f['self_complementarity'] * 0.15  # RIG-I/MDA-5
    s += f['pw_ifn'] * 0.25
    s += min(f['ifna'] / 20, 0.30)
    note = 'TLR7/8/9 (CpG, poly-U); RIG-I/MDA-5 (dsRNA). IFN-a Phase 1 monitoring (ICH S7B).'
    return RuleResult(min(s,1.0), 0.76, ['cpg_density','poly_u_density','self_complementarity'], note)

def score_cytokine_release(f: dict) -> RuleResult:
    s = min(f['il6'] / 10, 0.40)
    s += min(f['tnfa'] / 8, 0.30)
    s += f['pw_nfkb'] * 0.25
    s += f['pw_kyn'] * 0.20  # IDO1/Trp:Kyn pathway
    note = 'IDO1/Trp:Kyn (IL4I1-linked); MSD/Luminex; ASTCT CRS grading. IL-6 primary analyte.'
    return RuleResult(min(s,1.0), 0.72, ['il6','tnfa','pw_nfkb','pw_kyn'], note)

def score_hepatic_vacuolation(f: dict) -> RuleResult:
    s = f['is_lnp'] * 0.55  # LNP lysosomal vacuolation
    s += f['is_galnac'] * 0.25
    s += min(f['dose_mg_kg'] / 300, 0.20)
    note = 'LNP/GalNAc lysosomal vacuolation (Oil Red O, LAMP1 IHC). GLP reversibility satellite.'
    return RuleResult(min(s,1.0), 0.74, ['is_lnp','is_galnac','dose_mg_kg'], note)

def score_renal_accumulation(f: dict) -> RuleResult:
    s = min(f['ps_fraction'] * 0.45, 0.40)
    s += min(f['dose_mg_kg'] / 250, 0.20)
    note = 'PS-ASO lysosomal renal accumulation (H&E+PAS+LAMP1). Rat caveat: ICH S7A.'
    return RuleResult(min(s,1.0), 0.67, ['ps_fraction','dose_mg_kg'], note)

def score_coagulation(f: dict) -> RuleResult:
    s = min(f['ps_fraction'] * 0.50, 0.45)
    s += min(f['platelet_drop_fraction'] * 0.20, 0.15)
    note = 'PS APTT prolongation (nusinersen label); ATII binding. Monitor APTT/PT/fibrinogen.'
    return RuleResult(min(s,1.0), 0.74, ['ps_fraction','platelet_drop_fraction'], note)

def score_kupffer_activation(f: dict) -> RuleResult:
    s = f['is_lnp'] * 0.35
    s += min(f['ps_fraction'] * 0.20, 0.20)
    s += f['pw_nfkb'] * 0.15
    note = 'CD68/F4/80 IHC; Kupffer activation correlates with LNP dose and PS class.'
    return RuleResult(min(s,1.0), 0.70, ['is_lnp','ps_fraction','pw_nfkb'], note)

def score_cns_toxicity(f: dict) -> RuleResult:
    s = f['is_it'] * 0.35
    s += f['is_lncrna_hub'] * 0.40  # lncRNA hub/TDP-43 risk
    s += f['human_specific'] * 0.15
    note = 'IT delivery TDP-43 neuroinflam; lncRNA hub D1 risk. NfL/GFAP CSF+serum (FDA 2022).'
    return RuleResult(min(s,1.0), 0.77, ['is_it','is_lncrna_hub','human_specific'], note)

def score_offtarget_silencing(f: dict) -> RuleResult:
    s = min(f['cpg_density'] * 1.5, 0.25)
    s += f['is_lncrna_hub'] * 0.35
    s += f['micropeptide_risk'] * 0.20
    s += f['self_complementarity'] * 0.10
    note = 'Seed-region OT (pos 2-8); lncRNA hub D2 network; micropeptide co-depletion.'
    return RuleResult(min(s,1.0), 0.72, ['is_lncrna_hub','cpg_density','micropeptide_risk'], note)

def score_injection_site(f: dict) -> RuleResult:
    s = min(f['ps_fraction'] * 0.25, 0.20)
    s += f['is_galnac'] * 0.15
    s += min(f['dose_mg_kg'] / 200, 0.15)
    note = 'SC ISR per inclisiran/nusinersen label; H&E+IHC grading.'
    return RuleResult(min(s,1.0), 0.62, ['ps_fraction','is_galnac','dose_mg_kg'], note)

RULE_MAP = {
    ToxEndpoint.HEPATOCELLULAR_INJURY:    score_hepatocellular_injury,
    ToxEndpoint.RENAL_TUBULAR_INJURY:     score_renal_tubular_injury,
    ToxEndpoint.THROMBOCYTOPENIA:         score_thrombocytopenia,
    ToxEndpoint.COMPLEMENT_ACTIVATION:    score_complement_activation,
    ToxEndpoint.INNATE_IMMUNE_ACTIVATION: score_innate_immune,
    ToxEndpoint.CYTOKINE_RELEASE:         score_cytokine_release,
    ToxEndpoint.HEPATIC_VACUOLATION:      score_hepatic_vacuolation,
    ToxEndpoint.RENAL_ACCUMULATION:       score_renal_accumulation,
    ToxEndpoint.COAGULATION_EFFECT:       score_coagulation,
    ToxEndpoint.KUPFFER_CELL_ACTIVATION:  score_kupffer_activation,
    ToxEndpoint.CNS_TOXICITY:             score_cns_toxicity,
    ToxEndpoint.OFFTARGET_GENE_SILENCING: score_offtarget_silencing,
    ToxEndpoint.INJECTION_SITE_REACTION:  score_injection_site,
}

class MechanisticRulesEngine:
    def score_all(self, features: dict, inp: OligoToxInput) -> Dict[str, ToxPrediction]:
        results = {}
        for ep, rule_fn in RULE_MAP.items():
            rr = rule_fn(features)
            results[ep.value] = ToxPrediction(
                endpoint=ep, risk_score=round(rr.score, 4),
                risk_class=classify_risk(rr.score),
                confidence=rr.confidence,
                primary_domain=ENDPOINT_DOMAIN.get(ep, RiskDomain.D3),
                top_features=rr.features, mechanistic_note=rr.note)
        return results
