"""GSEApy pathway enrichment integration"""
from dataclasses import dataclass, field
from typing import List, Dict

HALLMARK = {
    'HALLMARK_TNFA_NFKB':       ['TNFAIP3','NFKBIA','ICAM1','IL6','CXCL8','CCL2','PTGS2','BIRC3'],
    'HALLMARK_IFN_ALPHA':        ['IFIT1','IFIT2','IFIT3','ISG15','MX1','OAS1','OAS2','IRF7','RSAD2'],
    'HALLMARK_IFN_GAMMA':        ['STAT1','STAT2','IRF1','CXCL10','HLA-A','GBP1','GBP2'],
    'HALLMARK_IL6_JAK_STAT3':    ['IL6','IL6R','JAK1','STAT3','SOCS3','IL10'],
    'HALLMARK_APOPTOSIS':        ['CASP3','CASP7','CASP8','CASP9','BAX','BCL2L11','PARP1','CYCS'],
    'HALLMARK_P53':              ['TP53','CDKN1A','MDM2','BAX','GADD45A','FAS'],
    'HALLMARK_COMPLEMENT':       ['C1QA','C1QB','C3','C4A','C5','CFB','CFD','CD55','CD59'],
    'HALLMARK_INFLAMMATORY':     ['IL1A','IL1B','IL6','IL8','CXCL1','PTGS2','NOS2','CCL5'],
    'HALLMARK_UPR':              ['XBP1','ATF6','DDIT3','EIF2AK3','HSPA5','HERPUD1'],
    'REACTOME_INNATE_IMMUNE':    ['TLR7','TLR8','TLR9','RIGI','MDA5','MAVS','IRF3','IRF7'],
    'REACTOME_COMPLEMENT':       ['C1QA','C3','C5','CFB','CFD','C4BPA','MASP1','MASP2'],
    'REACTOME_PLATELET':         ['GP6','ITGA2B','ITGB3','VWF','TBXA2R','P2RY12','PF4'],
    'REACTOME_KYNURENINE':       ['IDO1','IDO2','TDO2','KMO','AFMID','HAAO','IL4I1'],
    'REACTOME_FATTY_ACID':       ['ACSL1','FASN','HADHA','PPARA','CPT1A'],
    'HALLMARK_DNA_REPAIR':       ['BRCA1','BRCA2','RAD51','MSH2','ATM','CHEK2'],
    'HALLMARK_OXIDATIVE_PHOS':   ['NDUFA1','UQCRC1','COX5A','ATP5A1','SDHA'],
}

@dataclass
class PathwayResult:
    pathway: str; n_overlap: int; overlap_genes: List[str]; enrichment_score: float

@dataclass
class GSEAResult:
    deg_list: List[str]; pathway_scores: Dict[str,float]
    top_pathways: List[PathwayResult]; gsea_available: bool = False
    total_pathways_tested: int = 0

def score_pathways(deg_list: List[str], gene_sets: Dict = None) -> GSEAResult:
    if gene_sets is None: gene_sets = HALLMARK
    deg_upper = set(g.upper() for g in deg_list)
    pathway_scores = {}; top_pathways = []
    gsea_avail = False
    try:
        import gseapy as gp
        enr = gp.enrichr(gene_list=list(deg_upper), gene_sets=['MSigDB_Hallmark_2020'],
                         organism='Human', outdir=None, verbose=False)
        if enr is not None and hasattr(enr,'results'):
            for _,row in enr.results.iterrows():
                pw = row.get('Term',''); pval = float(row.get('P-value',1.0))
                pathway_scores[pw[:40]] = round(max(0.0,1.0-pval*10),4)
            gsea_avail = True
    except Exception:
        pass
    if not gsea_avail:
        for pw_name, pw_genes in gene_sets.items():
            pw_upper = set(g.upper() for g in pw_genes)
            overlap = deg_upper & pw_upper
            if not overlap: score = 0.0
            else:
                jaccard = len(overlap)/(len(pw_upper|deg_upper))
                recall = len(overlap)/max(len(pw_upper),1)
                score = min(jaccard*2.0 + recall*1.5, 1.0)
            pathway_scores[pw_name] = round(score,4)
            if overlap:
                top_pathways.append(PathwayResult(pw_name, len(overlap), sorted(overlap)[:5], score))
    top_pathways.sort(key=lambda x:-x.enrichment_score)
    return GSEAResult(deg_list=deg_list, pathway_scores=pathway_scores,
                      top_pathways=top_pathways[:10], gsea_available=gsea_avail,
                      total_pathways_tested=len(gene_sets))

def summary(r: GSEAResult) -> str:
    lines = [f'  [GSEApy] Pathway enrichment ({"GSEApy/Enrichr" if r.gsea_available else "keyword overlap"}) [{r.total_pathways_tested} sets]',
             f'  DEGs: {len(r.deg_list)}',
             f'  {"Pathway":<45}{"Score":>8}  Overlap']
    lines.append('  ' + '-'*72)
    for p in r.top_pathways[:8]:
        lines.append(f'  {p.pathway:<45}{p.enrichment_score:>8.4f}  {",".join(p.overlap_genes[:3])}')
    return '\n'.join(lines)
