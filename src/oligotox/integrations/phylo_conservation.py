"""PhyloP Conservation Layer — species correction modulation"""
from dataclasses import dataclass, field
from typing import Dict, Optional

GENE_CONSERVATION = {
    'TP53':0.95,'PCSK9':0.88,'TTR':0.90,'ALAS1':0.85,'APOB':0.82,'DMD':0.78,
    'SMN2':0.70,'SMN1':0.72,'NEAT1':0.40,'MALAT1':0.45,'HOTAIR':0.30,'XIST':0.35,
    'APOC3':0.80,'APOA1':0.83,'HBB':0.95,'CFTR':0.88,'HUNTINGTIN':0.92,'SOD1':0.96,
}

@dataclass
class ConservationProfile:
    gene: str; conservation_score: float; is_lncrna: bool
    cross_species_validity: str; correction_note: str
    modulated_factors: Dict[str,Dict[str,float]] = field(default_factory=dict)

def get_conservation_profile(gene: str, base_corrections: Optional[Dict[str,float]] = None,
                              is_lncrna: bool = False) -> ConservationProfile:
    cons = GENE_CONSERVATION.get(gene.upper(), 0.40 if is_lncrna else 0.70)
    validity = 'high' if cons>=0.80 else ('moderate' if cons>=0.55 else 'low')
    if cons < 0.50 or is_lncrna:
        note = (f'{gene}: LOW conservation (PhyloP proxy={cons:.2f}). '
                'Verify cynomolgus ortholog before NHP study (ICH S6R1).')
    elif cons < 0.75:
        note = f'{gene}: MODERATE conservation ({cons:.2f}). Human transcriptomics essential.'
    else:
        note = f'{gene}: HIGH conservation ({cons:.2f}). Standard IND bridging applies.'
    modulated = {}
    if base_corrections:
        for ep, factor in base_corrections.items():
            uncertainty = 1.0 - cons
            mod = factor*(1.0-uncertainty*0.4) + 0.5*(uncertainty*0.4)
            modulated[ep] = {'base': factor, 'conservation_modulated': round(min(max(mod,0.20),1.0),4),
                             'delta': round(min(max(mod,0.20),1.0)-factor,4)}
    return ConservationProfile(gene=gene, conservation_score=cons, is_lncrna=is_lncrna,
        cross_species_validity=validity, correction_note=note, modulated_factors=modulated)

def summary(cp: ConservationProfile) -> str:
    lines = ['  [PhyloP] Conservation-based species correction',
             f'  Gene:{cp.gene}  Conservation(proxy):{cp.conservation_score:.2f}  '
             f'Cross-species validity:{cp.cross_species_validity.upper()}',
             f'  Note: {cp.correction_note}']
    if cp.modulated_factors:
        lines.append(f'  {"Endpoint":<38}{"Base":>8}{"Modulated":>12}{"Delta":>8}')
        lines.append('  '+'-'*68)
        for ep,v in list(cp.modulated_factors.items())[:6]:
            lines.append(f'  {ep:<38}{v["base"]:>8.3f}{v["conservation_modulated"]:>12.3f}{v["delta"]:>8.3f}')
    return '\n'.join(lines)
