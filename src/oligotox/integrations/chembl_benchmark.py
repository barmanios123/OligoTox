"""ChEMBL/SIDER benchmark against approved oligonucleotide drugs"""
from dataclasses import dataclass, field
from typing import Dict, List, Optional

APPROVED_DRUG_PROFILES = {
    'nusinersen':   dict(chembl='CHEMBL3545177', modality='ASO_SSO', backbone='full_PS',
        delivery='intrathecal', approved=2016,
        known_tox={'coagulation_effect':0.75,'renal_tubular_injury':0.40,'thrombocytopenia':0.35,
                   'injection_site_reaction':0.20,'cns_toxicity':0.15,'hepatocellular_injury':0.10},
        note='Black box: thrombocytopenia, coagulation. Monitor APTT+CBC.'),
    'patisiran':    dict(chembl='CHEMBL4297610', modality='siRNA_LNP', backbone='partial_PS',
        delivery='LNP_IV', approved=2018,
        known_tox={'complement_activation':0.85,'cytokine_release':0.60,'innate_immune_activation':0.55,
                   'hepatic_vacuolation':0.40,'injection_site_reaction':0.30,'thrombocytopenia':0.20},
        note='CARPA: mandatory premedication. Patisiran NDA CARPA framework.'),
    'givosiran':    dict(chembl='CHEMBL4523760', modality='siRNA_GalNAc', backbone='partial_PS',
        delivery='GalNAc_SC', approved=2019,
        known_tox={'renal_tubular_injury':0.65,'hepatocellular_injury':0.55,'injection_site_reaction':0.40,
                   'cytokine_release':0.25,'innate_immune_activation':0.20},
        note='Hepatotoxicity warning. Renal monitoring mandatory.'),
    'inclisiran':   dict(chembl='CHEMBL4297613', modality='siRNA_GalNAc', backbone='partial_PS',
        delivery='GalNAc_SC', approved=2021,
        known_tox={'injection_site_reaction':0.85,'renal_tubular_injury':0.25,'hepatocellular_injury':0.15,'thrombocytopenia':0.10},
        note='ISR most common AE (~8%). Favorable hepatic/renal profile.'),
    'mipomersen':   dict(chembl='CHEMBL2103872', modality='ASO_gapmer', backbone='full_PS',
        delivery='naked_SC', approved=2013,
        known_tox={'hepatocellular_injury':0.80,'hepatic_vacuolation':0.70,'injection_site_reaction':0.90,
                   'renal_tubular_injury':0.35,'thrombocytopenia':0.30,'coagulation_effect':0.40},
        note='REMS: hepatotoxicity. High ISR. Full-PS gapmer class ref.'),
    'inotersen':    dict(chembl='CHEMBL4297599', modality='ASO_gapmer', backbone='full_PS',
        delivery='naked_SC', approved=2018,
        known_tox={'thrombocytopenia':0.90,'renal_tubular_injury':0.75,'hepatocellular_injury':0.50,
                   'coagulation_effect':0.45,'injection_site_reaction':0.60,'complement_activation':0.35},
        note='BLACK BOX: thrombocytopenia+glomerulonephritis. REMS. PS-ASO benchmark.'),
    'eplontersen':  dict(chembl='CHEMBL4650417', modality='ASO_gapmer_GalNAc',
        delivery='GalNAc_SC', approved=2023,
        known_tox={'injection_site_reaction':0.45,'renal_tubular_injury':0.30,'hepatocellular_injury':0.20,'thrombocytopenia':0.15},
        note='GalNAc reduces systemic PS; improved safety vs inotersen.'),
    'viltolarsen':  dict(chembl='CHEMBL4297656', modality='ASO_SSO_PMO',
        delivery='naked_IV', approved=2020,
        known_tox={'renal_tubular_injury':0.30,'injection_site_reaction':0.25,'coagulation_effect':0.05,'thrombocytopenia':0.05},
        note='PMO backbone: minimal PS effects. Lower tox profile benchmark.'),
}

@dataclass
class BenchmarkComparison:
    drug_name: str; predicted_scores: Dict; known_scores: Dict
    per_endpoint_error: Dict; mae: float; rmse: float
    correct_risk_class: Dict; accuracy: float; regulatory_note: str = ''

def benchmark_against_drug(predicted: Dict[str,float], drug_name: str) -> Optional[BenchmarkComparison]:
    drug = APPROVED_DRUG_PROFILES.get(drug_name.lower())
    if drug is None: return None
    known = drug['known_tox']; endpoints = set(predicted.keys()) & set(known.keys())
    if not endpoints: return None
    errors = {}; correct = {}
    def rc(s): return 'high' if s>=0.60 else ('moderate' if s>=0.30 else 'low')
    for ep in endpoints:
        errors[ep] = abs(predicted[ep]-known[ep])
        correct[ep] = rc(predicted[ep]) == rc(known[ep])
    mae = sum(errors.values())/len(errors)
    rmse = (sum(v**2 for v in errors.values())/len(errors))**0.5
    acc = sum(correct.values())/max(len(correct),1)
    return BenchmarkComparison(drug_name=drug_name, predicted_scores=predicted, known_scores=known,
        per_endpoint_error=errors, mae=round(mae,4), rmse=round(rmse,4),
        correct_risk_class=correct, accuracy=round(acc,4),
        regulatory_note=drug.get('note',''))

def summary(r: BenchmarkComparison) -> str:
    lines = [f'  [ChEMBL] Benchmark: {r.drug_name.upper()}',
             f'  MAE:{r.mae:.4f} RMSE:{r.rmse:.4f} Risk-class accuracy:{r.accuracy:.1%}',
             f'  {"Endpoint":<38}{"Pred":>8}{"Known":>8}{"Err":>8}  Match']
    lines.append('  ' + '-'*68)
    for ep in sorted(r.per_endpoint_error,key=r.per_endpoint_error.get,reverse=True):
        lines.append(f'  {ep:<38}{r.predicted_scores.get(ep,0):>8.3f}'
                     f'{r.known_scores.get(ep,0):>8.3f}'
                     f'{r.per_endpoint_error[ep]:>8.3f}  '
                     f'{"✓" if r.correct_risk_class.get(ep) else "✗"}')
    lines.append(f'  Regulatory note: {r.regulatory_note}')
    return '\n'.join(lines)
