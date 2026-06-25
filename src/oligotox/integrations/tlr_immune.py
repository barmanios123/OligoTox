"""TLR/Innate Immune Sequence Scorer — upgrades D3 innate module"""
import re
from dataclasses import dataclass, field
from typing import Dict

TLR7_MOTIFS  = [r'U{2,}G', r'UGA', r'UCA', r'UGUGU', r'U{3,}']
TLR8_MOTIFS  = [r'UUGU', r'UUAU', r'UGAU', r'UU[AGC]U']
TLR9_MOTIFS  = [r'GTCGTT', r'GACGTC', r'AACGTT', r'AGCGTT']
PS_TLR9_BOOST = 1.8

@dataclass
class TLRScoreResult:
    sequence: str; tlr7_score: float; tlr8_score: float; tlr9_score: float
    rigi_score: float; mda5_score: float; composite_innate_score: float
    motif_hits: Dict[str,int] = field(default_factory=dict)
    dominant_receptor: str = ''; mechanistic_note: str = ''

def score_tlr_innate(seq: str, ps_fraction: float = 0.0, is_dsrna: bool = False) -> TLRScoreResult:
    s = seq.upper().replace('T','U'); L = max(len(s),1); hits = {}
    tlr7_h = sum(len(re.findall(p,s)) for p in TLR7_MOTIFS); hits['tlr7'] = tlr7_h
    tlr7 = min(tlr7_h / max(L*0.15,1), 1.0)
    tlr8_h = sum(len(re.findall(p,s)) for p in TLR8_MOTIFS); hits['tlr8'] = tlr8_h
    tlr8 = min(tlr8_h / max(L*0.10,1), 1.0)
    s_dna = s.replace('U','T')
    tlr9_h = sum(len(re.findall(p,s_dna)) for p in TLR9_MOTIFS)
    cpg = s_dna.count('CG') / max(L-1,1); hits['tlr9'] = tlr9_h; hits['cpg'] = round(cpg,4)
    tlr9 = min((tlr9_h*0.25 + cpg*2.0) * (1 + ps_fraction*(PS_TLR9_BOOST-1)), 1.0)
    poly_u = len(re.findall(r'U{4,}',s)); hits['poly_u4'] = poly_u
    rigi = min((float(is_dsrna)*0.50 + poly_u*0.10 + ps_fraction*0.10), 1.0)
    gc = (s.count('G')+s.count('C'))/L; gc_runs = len(re.findall(r'[GC]{6,}',s))
    mda5 = min(gc*0.40 + gc_runs*0.15, 1.0) if (is_dsrna or L>100) else 0.0
    composite = tlr7*0.30 + tlr8*0.20 + tlr9*0.25 + rigi*0.15 + mda5*0.10
    dominant = max({'TLR7':tlr7,'TLR8':tlr8,'TLR9':tlr9,'RIG-I':rigi,'MDA-5':mda5}, key=lambda k: {'TLR7':tlr7,'TLR8':tlr8,'TLR9':tlr9,'RIG-I':rigi,'MDA-5':mda5}[k])
    return TLRScoreResult(sequence=seq, tlr7_score=round(tlr7,3), tlr8_score=round(tlr8,3),
        tlr9_score=round(tlr9,3), rigi_score=round(rigi,3), mda5_score=round(mda5,3),
        composite_innate_score=round(composite,3), motif_hits=hits, dominant_receptor=dominant,
        mechanistic_note=f'Dominant: {dominant}. PS amplifies TLR9 {PS_TLR9_BOOST}x (Vollmer 2004).')

def summary(r: TLRScoreResult) -> str:
    return '\n'.join(['  [TLR-Innate] Sequence-based innate immune scorer',
        f'  TLR7:{r.tlr7_score:.3f} TLR8:{r.tlr8_score:.3f} TLR9:{r.tlr9_score:.3f} '
        f'RIG-I:{r.rigi_score:.3f} MDA-5:{r.mda5_score:.3f}',
        f'  Composite: {r.composite_innate_score:.3f}  Dominant: {r.dominant_receptor}',
        f'  Note: {r.mechanistic_note}'])
