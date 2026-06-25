"""Rules Engine v2 — wires all integration modules"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from typing import Dict, Optional, Any
from oligotox.rules.engine import MechanisticRulesEngine
from oligotox.features.engineer import extract_features
from oligotox.schemas.data_models import OligoInput, DeliveryInput, OmicsReadout, OligoToxInput, PlateletComplement, OffTargetContext

try:
    from integrations.seed_match import score_seed_matches
    SEED_MATCH_AVAIL = True
except Exception: SEED_MATCH_AVAIL = False

try:
    from integrations.tlr_immune import score_tlr_innate
    TLR_AVAIL = True
except Exception: TLR_AVAIL = False

try:
    from integrations.gsea_pathway import score_pathways
    GSEA_AVAIL = True
except Exception: GSEA_AVAIL = False

try:
    from integrations.phylo_conservation import get_conservation_profile
    PHYLO_AVAIL = True
except Exception: PHYLO_AVAIL = False

try:
    from integrations.geneformer_context import get_cell_context
    GENEFORMER_AVAIL = True
except Exception: GENEFORMER_AVAIL = False

_eng = MechanisticRulesEngine()

def run_engine_v2(oligo: OligoInput, delivery: DeliveryInput,
                  platelet: Optional[PlateletComplement] = None,
                  omics: Optional[OmicsReadout] = None,
                  target_gene: str = '', target_tissue: str = 'hepatocyte',
                  is_lncrna: bool = False) -> Dict[str, Any]:
    if platelet is None: platelet = PlateletComplement()
    if omics is None: omics = OmicsReadout()
    inp = OligoToxInput(oligo=oligo, delivery=delivery, platelet_complement=platelet, omics=omics)
    feats = extract_features(inp)
    base_preds = _eng.score_all(feats, inp)
    scores = {ep_name: pred.risk_score for ep_name, pred in base_preds.items()}

    meta: Dict[str, Any] = {'seed_match':None,'tlr_immune':None,'gsea':None,'phylo':None,'geneformer':None}

    if SEED_MATCH_AVAIL and oligo.sequence:
        sm = score_seed_matches(oligo.sequence)
        meta['seed_match'] = sm
        old = scores.get('offtarget_gene_silencing', 0.0)
        scores['offtarget_gene_silencing'] = round(min(old*0.50 + sm.d2_risk_contribution*0.50, 1.0), 4)
        if sm.dsrna_probability > 0.20:
            old_ii = scores.get('innate_immune_activation', 0.0)
            scores['innate_immune_activation'] = round(min(old_ii + sm.dsrna_probability*0.12, 1.0), 4)

    if TLR_AVAIL and oligo.sequence:
        ps_frac = feats.get('ps_fraction', 0.0)
        is_ds = oligo.modality.value.lower() in ('sirna_lnp','sirna_galnac','mirna_mimic')
        tlr = score_tlr_innate(oligo.sequence, ps_fraction=ps_frac, is_dsrna=is_ds)
        meta['tlr_immune'] = tlr
        scores['innate_immune_activation'] = round(min(scores.get('innate_immune_activation',0)*0.55 + tlr.composite_innate_score*0.45, 1.0), 4)
        scores['cytokine_release'] = round(min(scores.get('cytokine_release',0)*0.65 + (tlr.tlr7_score+tlr.tlr8_score)*0.175, 1.0), 4)

    if GSEA_AVAIL and omics.deg_list:
        gsea = score_pathways(omics.deg_list)
        meta['gsea'] = gsea
        pw = gsea.pathway_scores
        scores['innate_immune_activation'] = round(min(scores.get('innate_immune_activation',0)*0.65 + pw.get('HALLMARK_IFN_ALPHA',0)*0.35, 1.0), 4)
        scores['cytokine_release'] = round(min(scores.get('cytokine_release',0)*0.60 + pw.get('HALLMARK_TNFA_NFKB',0)*0.30 + pw.get('REACTOME_KYNURENINE',0)*0.10, 1.0), 4)
        scores['hepatocellular_injury'] = round(min(scores.get('hepatocellular_injury',0)*0.70 + pw.get('HALLMARK_APOPTOSIS',0)*0.30, 1.0), 4)
        scores['complement_activation'] = round(min(scores.get('complement_activation',0)*0.70 + pw.get('REACTOME_COMPLEMENT',pw.get('HALLMARK_COMPLEMENT',0))*0.30, 1.0), 4)
        scores['thrombocytopenia'] = round(min(scores.get('thrombocytopenia',0)*0.70 + pw.get('REACTOME_PLATELET',0)*0.30, 1.0), 4)

    if PHYLO_AVAIL and target_gene:
        phylo = get_conservation_profile(target_gene, base_corrections=scores, is_lncrna=is_lncrna)
        meta['phylo'] = phylo
        for ep, vals in phylo.modulated_factors.items():
            if ep in scores: scores[ep] = vals['conservation_modulated']

    if GENEFORMER_AVAIL and target_gene and is_lncrna:
        hub_score = scores.get('offtarget_gene_silencing', 0.5)
        gf = get_cell_context(target_gene, target_tissue, lncrna_network_score=hub_score)
        meta['geneformer'] = gf
        for ep in gf.tissue_relevant_endpoints:
            if ep in scores:
                scores[ep] = round(min(scores[ep]*0.70 + gf.hub_connectivity_score*0.30, 1.0), 4)

    return {'scores': {k: round(float(v),4) for k,v in scores.items()},
            'base_predictions': base_preds, 'integration_meta': meta,
            'integrations_active': {k: v is not None for k, v in meta.items()}}
