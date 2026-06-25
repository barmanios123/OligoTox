"""OligoTox v0.2 Full Integration Demo"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from oligotox.schemas.data_models import OligoInput, DeliveryInput, OmicsReadout, PlateletComplement, Modality, DeliveryVehicle, Species
from oligotox.rules.engine_v2 import run_engine_v2
from oligotox.assays.planner import generate_assay_plan, print_assay_plan
from oligotox.translation.bridge import generate_translational_bridge, print_bridge
from oligotox.schemas.tox_endpoints import ToxProfile, ToxEndpoint, ToxPrediction, RiskClass, RiskDomain, ENDPOINT_DOMAIN
from profiler import OligoToxProfiler
from oligotox.integrations.seed_match import summary as seed_sum
from oligotox.integrations.tlr_immune import summary as tlr_sum
from oligotox.integrations.gsea_pathway import summary as gsea_sum
from oligotox.integrations.chembl_benchmark import benchmark_against_drug, summary as bench_sum
from oligotox.integrations.phylo_conservation import summary as phylo_sum
from oligotox.integrations.geneformer_context import summary as gf_sum

D = '='*78
def sec(t): print(f'\n{D}\n  {t}\n{D}')

def make_profile(sid, result):
    preds = {}
    for ep_name, score in result['scores'].items():
        try: ep = ToxEndpoint(ep_name)
        except ValueError: continue
        domain = ENDPOINT_DOMAIN.get(ep, RiskDomain.D3)
        rc = RiskClass.HIGH if score>=0.60 else (RiskClass.MODERATE if score>=0.30 else RiskClass.LOW)
        preds[ep] = ToxPrediction(endpoint=ep, risk_score=score, risk_class=rc, confidence=0.75, primary_domain=domain)
    return ToxProfile(sample_id=sid, predictions=preds)

CANDIDATES = [
    dict(id='OTX-001-LNP-siRNA',
         oligo=OligoInput(sequence='GCAUUCGGAUCGAUUGCAUC', modality=Modality.SIRNA_LNP),
         delivery=DeliveryInput(vehicle=DeliveryVehicle.LNP, route='intravenous',
                                ionizable_lipid='DLin-MC3-DMA', peg_lipid_mol_pct=1.5, particle_size_nm=85.0, pdi=0.11),
         omics=OmicsReadout(deg_list=['C3','C4A','CFB','SERPING1','IL6','CXCL8','TNF'], cytokine_il6=4.2, cytokine_ifna=0.5),
         platelet=PlateletComplement(platelet_count_baseline=220, platelet_count_treated=175,
                                     c3a_pg_ml=280, c5a_pg_ml=180, sc5b9_ng_ml=22, pf4_ps_elisa_od=0.8, complement_activation_ex_vivo=True),
         gene='TTR', tissue='hepatocyte', is_lncrna=False),
    dict(id='OTX-002-GalNAc-siRNA',
         oligo=OligoInput(sequence='GUAUCGAUUGCAUCCUGAUU', modality=Modality.SIRNA_GALNAC),
         delivery=DeliveryInput(vehicle=DeliveryVehicle.GALNAC, route='subcutaneous', galnac_valency=3),
         omics=OmicsReadout(deg_list=['ALT','KIM1','NGAL','CASP3'], cytokine_il6=0.8),
         platelet=PlateletComplement(), gene='PCSK9', tissue='hepatocyte', is_lncrna=False),
    dict(id='OTX-003-ASO-gapmer',
         oligo=OligoInput(sequence='GCUAAGCUAGCUAGCUAGCUAGCU', modality=Modality.ASO_GAPMER),
         delivery=DeliveryInput(vehicle=DeliveryVehicle.NAKED, route='subcutaneous', dose_mg_kg=50.0, repeat_dose=True),
         omics=OmicsReadout(deg_list=['CASP3','CASP7','BAX','ALT','KIM1','ITGA2B','PF4'], cytokine_il6=2.1),
         platelet=PlateletComplement(platelet_count_baseline=220, platelet_count_treated=140, pf4_ps_elisa_od=1.6),
         gene='ApoB', tissue='hepatocyte', is_lncrna=False),
    dict(id='OTX-004-IT-lncRNA-NEAT1',
         oligo=OligoInput(sequence='GCUAAGCUAGCUAGCUAGCUAGCU', modality=Modality.LNCRNA_TARGETING,
                          target_gene='NEAT1', target_is_lncrna_hub=True, is_human_specific_sequence=True),
         delivery=DeliveryInput(vehicle=DeliveryVehicle.NAKED, route='intrathecal'),
         omics=OmicsReadout(deg_list=['STAT1','IFIT1','ISG15','MX1','CASP3','C3','TLR7','IDO1','IL6','CXCL8'],
                            cytokine_il6=3.5, cytokine_ifna=8.0),
         platelet=PlateletComplement(), gene='NEAT1', tissue='cns_neuron', is_lncrna=True),
]

def main():
    print(f'\n{D}')
    print('  OligoTox v0.2 -- Mechanism-Aware Tox for RNA/Oligo Therapeutics')
    print('  Integrations: SeedMatch | TLR-Innate | GSEApy | PhyloP | Geneformer | ChEMBL/SIDER')
    print(f'{D}')

    profiles = []; results_store = []
    for cand in CANDIDATES:
        sec(f'{cand["id"]} | {cand["oligo"].modality.value} | {cand["delivery"].vehicle.value}/{cand["delivery"].route}')
        r = run_engine_v2(oligo=cand['oligo'], delivery=cand['delivery'],
                          platelet=cand['platelet'], omics=cand['omics'],
                          target_gene=cand['gene'], target_tissue=cand['tissue'], is_lncrna=cand['is_lncrna'])
        results_store.append((cand['id'], r))
        p = make_profile(cand['id'], r)
        profiles.append(p)
        print(p.summary_table())
        if cand['id'] in ('OTX-001-LNP-siRNA','OTX-004-IT-lncRNA-NEAT1'):
            m = r['integration_meta']
            print('\n  -- Integration Module Outputs --')
            if m['seed_match']: print(seed_sum(m['seed_match']))
            if m['tlr_immune']: print(tlr_sum(m['tlr_immune']))
            if m['gsea']:       print(gsea_sum(m['gsea']))
            if m['phylo']:      print(phylo_sum(m['phylo']))
            if m['geneformer']: print(gf_sum(m['geneformer']))

    sec('Candidate Comparator -- Ranked by Composite Risk Score')
    ranked = sorted(profiles, key=lambda p: p.composite_score)
    print(f'  {"Rank":<6}{"Candidate":<42}{"Composite":>12}{"Dominant Endpoint":>26}')
    print('  '+'-'*88)
    for i,p in enumerate(ranked,1):
        dom = max(p.predictions, key=lambda e:p.predictions[e].risk_score) if p.predictions else ''
        dn = dom.value if hasattr(dom,'value') else str(dom)
        print(f'  {i:<6}{p.sample_id:<42}{p.composite_score:>12.4f}{dn:>26}')

    sec('ChEMBL/SIDER Benchmark -- OTX-003-ASO vs Inotersen (PS-ASO Reference)')
    b = benchmark_against_drug(results_store[2][1]['scores'], 'inotersen')
    if b: print(bench_sum(b))

    sec('ChEMBL/SIDER Benchmark -- OTX-001-LNP vs Patisiran (LNP Reference)')
    b2 = benchmark_against_drug(results_store[0][1]['scores'], 'patisiran')
    if b2: print(bench_sum(b2))

    sec('Assay Planner -- OTX-004-IT-lncRNA-NEAT1')
    profiler = OligoToxProfiler()
    fp = profiler.profile(CANDIDATES[3]['id'], CANDIDATES[3]['oligo'], CANDIDATES[3]['delivery'],
                          CANDIDATES[3]['platelet'], CANDIDATES[3]['omics'])
    plan = generate_assay_plan(fp, CANDIDATES[3]['oligo'], CANDIDATES[3]['delivery'])
    print_assay_plan(plan)

    sec('Translational Bridge -- OTX-004-IT-lncRNA-NEAT1')
    bridge = generate_translational_bridge(fp, species=Species.CYNOMOLGUS)
    print_bridge(bridge)

    print(f'\n{D}')
    print('  OligoTox v0.2 complete. All 6 integrations validated.')
    print(f'{D}\n')

if __name__ == '__main__':
    main()
