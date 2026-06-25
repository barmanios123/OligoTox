"""Translational Bridge — preclinical to clinical biomarker mapping"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from dataclasses import dataclass, field
from typing import Dict, List
from oligotox.schemas.data_models import Species
from oligotox.schemas.tox_endpoints import ToxProfile

BRIDGE_MAP = {
    'hepatocellular_injury':    dict(preclinical='ALT/AST/GDH/H&E/casp-3', clinical='ALT/AST; Hy Law',
        note='ICH S2R1/S7A. NOAEL-based safety margin. Hy Law per FDA DILI 2009.'),
    'hepatic_vacuolation':      dict(preclinical='H&E; Oil Red O; LAMP1', clinical='ALT; liver imaging',
        note='GLP histopathology; reversibility satellite group.'),
    'renal_tubular_injury':     dict(preclinical='KIM-1/NGAL/clusterin; H&E', clinical='eGFR; ACR; urine KIM-1',
        note='Do not use rat renal for human extrapolation (ICH S7A).'),
    'renal_accumulation':       dict(preclinical='H&E+PAS+LAMP1; tissue TK', clinical='eGFR; urine protein',
        note='Same caveat as renal tubular injury.'),
    'thrombocytopenia':         dict(preclinical='CBC; PF4-PS; C3d/C4d flow', clinical='CBC; anti-PF4 Ab; C3d/C4d',
        note='Phase 1 CBC weekly; anti-PF4 if PLT drops >30%.'),
    'coagulation_effect':       dict(preclinical='APTT/PT/fibrinogen/ATII', clinical='APTT; PT; fibrinogen',
        note='Nusinersen APTT monitoring precedent.'),
    'complement_activation':    dict(preclinical='C3a/C5a/SC5b-9/CH50/AH50', clinical='C3/C4/CH50; infusion reaction',
        note='Patisiran IND CARPA framework; FDA IV oligo guidance.'),
    'innate_immune_activation': dict(preclinical='TLR7/8/9; RIG-I/MDA-5; IFN-a', clinical='IFN-a; PBMC; cytokines',
        note='IFN-a monitoring Phase 1 per ICH S7B guidance.'),
    'cytokine_release':         dict(preclinical='MSD/Luminex; Trp:Kyn; IDO1', clinical='IL-6/TNF-a; Trp:Kyn ratio',
        note='ASTCT CRS grading; IL-6 primary monitoring analyte.'),
    'kupffer_cell_activation':  dict(preclinical='CD68/F4/80 IHC; cytokines', clinical='CRP; ferritin; cytokines',
        note='Correlated with histopathology and biomarker data.'),
    'cns_toxicity':             dict(preclinical='NfL/GFAP CSF+serum; TDP-43 IHC', clinical='Serum NfL+GFAP; CSF',
        note='Nusinersen NDA NfL/GFAP precedent; FDA 2022 IT oligo guidance.'),
    'offtarget_gene_silencing': dict(preclinical='RNA-seq; seed alignment; OT RT-PCR', clinical='Target PD; transcriptomics',
        note='Human-specific sequences require human transcriptomics (ICH S6R1).'),
    'injection_site_reaction':  dict(preclinical='SC ISR H&E+IHC; cytokines', clinical='ISR scoring; biopsy if severe',
        note='ISR per inclisiran/nusinersen label precedent.'),
}

@dataclass
class BridgeEntry:
    endpoint: str; preclinical_marker: str; clinical_marker: str; risk: float; note: str

@dataclass
class TranslationalBridge:
    sample_id: str; species: str; entries: List[BridgeEntry] = field(default_factory=list)
    species_caveats: List[str] = field(default_factory=list)

def generate_translational_bridge(profile: ToxProfile, species: Species = Species.CYNOMOLGUS) -> TranslationalBridge:
    bridge = TranslationalBridge(sample_id=profile.sample_id, species=species.value)
    for ep, pred in sorted(profile.predictions.items(), key=lambda x:-x[1].risk_score):
        bm = BRIDGE_MAP.get(ep.value, {})
        bridge.entries.append(BridgeEntry(
            endpoint=ep.value, preclinical_marker=bm.get('preclinical','See label'),
            clinical_marker=bm.get('clinical','See label'),
            risk=pred.risk_score, note=bm.get('note','')))
    caveats = []
    if species == Species.CYNOMOLGUS:
        caveats += ['Cynomolgus ASGPR and complement closely model human.',
                    'Preferred species for complement, platelet, and hepatic endpoints.',
                    'DART: EFD MPS or surrogate sequence; pre-IND FDA/EMA consult.',
                    'lncRNA: verify cynomolgus ortholog; hub network analysis required.']
    elif species == Species.RAT:
        caveats += ['Rat renal tubular injury overestimates human (ICH S7A).',
                    'Rat TLR pattern differs from human; PBMC assay required.',
                    'Rat complement less predictive than cynomolgus.']
    bridge.species_caveats = caveats
    return bridge

def print_bridge(bridge: TranslationalBridge):
    print(f'  Translational Bridge | {bridge.sample_id} | Species: {bridge.species}')
    print(f'  {"Endpoint":<30}{"Preclinical":<24}{"Clinical":<18}{"Risk"}')
    print('  ' + '-' * 78)
    for e in bridge.entries:
        print(f'  {e.endpoint:<30}{e.preclinical_marker[:23]:<24}{e.clinical_marker[:17]:<18}{e.risk:.3f}')
        if e.note: print(f'    -> {e.note}')
    if bridge.species_caveats:
        print('\n  Species Caveats:')
        for c in bridge.species_caveats: print(f'    * {c}')
