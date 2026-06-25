"""Geneformer cell-type context (curated proxy; replace with HF model for v1.0)"""
from dataclasses import dataclass, field
from typing import Dict, List

TISSUE_EXPRESSION = {
    'hepatocyte':           {'PCSK9':1,'APOB':1,'TTR':1,'ALAS1':1,'ALB':1,'NEAT1':3,'MALAT1':2,'TLR7':4},
    'renal_proximal_tubule':{'AQP1':1,'SLC22A6':1,'CUBN':1,'NEAT1':3,'MALAT1':2,'TLR7':4},
    'cns_neuron':           {'SMN1':1,'HUNTINGTIN':2,'SOD1':2,'NEAT1':2,'MALAT1':1,'TDP43':1,'NfL':1,'GFAP':2},
    'pbmc':                 {'TLR7':1,'TLR8':1,'TLR9':2,'MYD88':1,'IRF7':1,'ISG15':2,'STAT1':1},
    'platelet':             {'GP6':1,'ITGA2B':1,'ITGB3':1,'PF4':1,'C3':3},
}
TISSUE_ENDPOINT_MAP = {
    'hepatocyte':            ['hepatocellular_injury','hepatic_vacuolation','kupffer_cell_activation'],
    'renal_proximal_tubule': ['renal_tubular_injury','renal_accumulation'],
    'cns_neuron':            ['cns_toxicity','offtarget_gene_silencing'],
    'pbmc':                  ['innate_immune_activation','cytokine_release'],
    'platelet':              ['thrombocytopenia','complement_activation'],
}

@dataclass
class GeneformerContext:
    target_gene: str; target_tissue: str; expression_rank: int
    hub_connectivity_score: float; context_note: str
    tissue_relevant_endpoints: List[str] = field(default_factory=list)

def get_cell_context(gene: str, tissue: str, lncrna_network_score: float = 0.0) -> GeneformerContext:
    td = TISSUE_EXPRESSION.get(tissue.lower(), {})
    rank = td.get(gene.upper(), 3)
    hub = lncrna_network_score
    if rank == 1: hub = min(hub*1.3, 1.0)
    elif rank >= 4: hub = hub*0.6
    eps = TISSUE_ENDPOINT_MAP.get(tissue.lower(), [])
    note = (f'{gene} in {tissue}: expr rank {rank}/5. Hub score: {hub:.3f}. '
            'Replace with ctheodoris/Geneformer (HuggingFace) for full embeddings.')
    return GeneformerContext(target_gene=gene, target_tissue=tissue, expression_rank=rank,
        hub_connectivity_score=round(hub,3), context_note=note,
        tissue_relevant_endpoints=eps)

def summary(ctx: GeneformerContext) -> str:
    return '\n'.join(['  [Geneformer] Cell-type context (curated proxy)',
        f'  Gene:{ctx.target_gene}  Tissue:{ctx.target_tissue}  '
        f'ExprRank:{ctx.expression_rank}/5  HubScore:{ctx.hub_connectivity_score:.3f}',
        f'  Relevant endpoints: {", ".join(ctx.tissue_relevant_endpoints)}',
        f'  Note: {ctx.context_note}'])
