"""Assay Planner — stage/modality/risk-adapted assay packages"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from dataclasses import dataclass, field
from typing import List, Dict
from oligotox.schemas.tox_endpoints import ToxProfile, RiskClass, ToxEndpoint
from oligotox.schemas.data_models import OligoInput, DeliveryInput, DeliveryVehicle, Modality

@dataclass
class AssayItem:
    assay: str; endpoint: str; species: str; priority: int; rationale: str

@dataclass
class AssayPlan:
    sample_id: str; stage: str; modality: str; delivery: str
    p1_assays: List[AssayItem] = field(default_factory=list)
    p2_assays: List[AssayItem] = field(default_factory=list)
    p3_assays: List[AssayItem] = field(default_factory=list)

def generate_assay_plan(profile: ToxProfile, oligo: OligoInput, delivery: DeliveryInput) -> AssayPlan:
    plan = AssayPlan(sample_id=profile.sample_id, stage='pre_IND',
                     modality=oligo.modality.value, delivery=delivery.vehicle.value)
    scores = {ep.value: pred.risk_score for ep, pred in profile.predictions.items()}
    mod = oligo.modality; veh = delivery.vehicle

    def add(priority, assay, endpoint, species, rationale):
        item = AssayItem(assay=assay, endpoint=endpoint, species=species, priority=priority, rationale=rationale)
        if priority == 1: plan.p1_assays.append(item)
        elif priority == 2: plan.p2_assays.append(item)
        else: plan.p3_assays.append(item)

    # Universal P1
    add(1,'28-day repeat-dose tox (rat+cynomolgus)','hepatocellular_injury','rat+cyno',
        'GLP histopathology: H&E+Masson. ICH S7A core package.')

    # PS/ASO class
    if scores.get('thrombocytopenia',0) >= 0.25 or mod in (Modality.ASO_GAPMER, Modality.ASO_SSO):
        add(1,'PF4-PS ELISA','thrombocytopenia','cynomolgus','PF4-PS complex; anti-PF4 Ab.')
        add(1,'Platelet-bound C3d/C4d flow cytometry','complement_activation','human ex vivo',
            'Complement-opsonized platelet clearance; CARPA mechanism.')
        add(1,'Trp:Kyn HPLC-MS','cytokine_release','rat/plasma','IDO1/kynurenine pathway.')

    # LNP
    if veh == DeliveryVehicle.LNP:
        add(1,'CH50/AH50 ex vivo human serum','complement_activation','human serum',
            'CARPA classical+alternative pathway. Patisiran NDA framework.')
        add(2,'RIG-I/MDA-5 reporter assay','innate_immune_activation','human PBMC',
            'dsRNA innate sensing. ICH S7B IFN-a monitoring.')
        add(2,'Liver H&E + Oil Red O','hepatic_vacuolation','rat',
            'LNP lysosomal vacuolation assessment.')

    # lncRNA targeting
    if mod == Modality.LNCRNA_TARGETING or oligo.target_is_lncrna_hub:
        add(1,'lncRNA hub connectivity + network RNA-seq','offtarget_gene_silencing','human cell line',
            'lncRNA hub disruption; D2 ncRNA off-target network (ICH S6R1).')
        add(1,'Micropeptide ELISA / ribosome profiling','offtarget_gene_silencing','human cell line',
            'Micropeptide co-depletion risk.')
        add(1,'Cynomolgus ortholog verification','offtarget_gene_silencing','cynomolgus',
            'Human-specific lncRNA bridging strategy.')

    # Intrathecal
    if 'intrathecal' in delivery.route.lower() or oligo.is_human_specific_sequence:
        add(1,'EFD microphysiological system (ICH S5.3)','cns_toxicity','MPS/human',
            'DART human-specific sequence alternative.')
        add(1,'Surrogate sequence FDA pre-IND consultation','cns_toxicity','rat/cyno',
            'DART bridging via surrogate; pre-IND meeting.')

    # General CNS
    if scores.get('cns_toxicity',0) >= 0.20:
        add(2,'NfL/GFAP CSF+serum ELISA','cns_toxicity','cynomolgus',
            'Nusinersen NDA precedent; FDA 2022 IT oligo guidance.')

    # Renal
    if scores.get('renal_tubular_injury',0) >= 0.30:
        add(2,'Renal H&E + PAS + LAMP1 IHC','renal_accumulation','rat',
            'PS-ASO lysosomal renal accumulation. Rat caveat.')

    # Recovery
    add(2,'Recovery satellite group (28d)','hepatocellular_injury','rat',
        'Reversibility of liver/kidney findings.')

    # P3
    add(3,'TLR7/8/9 reporter assay (HEK-Blue)','innate_immune_activation','human cell line',
        'Sequence-based TLRpred upgrade; upgrades regex to probability scores.')
    add(3,'immuneSIM repertoire simulation','innate_immune_activation','in silico',
        'Validate innate immune module against simulated immune repertoire.')
    add(3,'ChEMBL/SIDER DILI flag check','hepatocellular_injury','database',
        'Programmatic DILI flag from approved drug comparator database.')

    return plan

def print_assay_plan(plan: AssayPlan):
    D = '=' * 78
    print(f'  Assay Planner | {plan.sample_id} | {plan.modality} | {plan.delivery}')
    for priority, items, label in [(1,plan.p1_assays,'Essential'),(2,plan.p2_assays,'Recommended'),(3,plan.p3_assays,'Exploratory')]:
        if not items: continue
        print(f'\n  P{priority} -- {label} ({len(items)} assays)')
        print(f'  {"Assay":<48}{"Endpoint":<24}{"Species"}')
        print('  ' + '-' * 78)
        for a in items:
            print(f'  {a.assay:<48}{a.endpoint[:23]:<24}{a.species}')
            print(f'    >> {a.rationale}')
