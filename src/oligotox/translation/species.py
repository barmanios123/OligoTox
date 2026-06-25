"""Species correction factors per endpoint"""
from typing import Dict
SPECIES_CORRECTION = {
    'rat': {
        'renal_tubular_injury':0.35,'renal_accumulation':0.40,'hepatocellular_injury':0.75,
        'thrombocytopenia':0.60,'complement_activation':0.50,'innate_immune_activation':0.65,
        'coagulation_effect':0.70,'cns_toxicity':0.80,'offtarget_gene_silencing':0.55,
    },
    'mouse': {
        'renal_tubular_injury':0.30,'complement_activation':0.35,'thrombocytopenia':0.45,
        'hepatocellular_injury':0.65,'offtarget_gene_silencing':0.45,
    },
    'cynomolgus': {
        'complement_activation':0.90,'thrombocytopenia':0.85,'hepatocellular_injury':0.88,
        'renal_tubular_injury':0.80,'coagulation_effect':0.85,'innate_immune_activation':0.82,
        'offtarget_gene_silencing':0.75,'cns_toxicity':0.78,
    },
}

def apply_species_correction(scores: Dict[str,float], species: str) -> Dict[str,float]:
    corrections = SPECIES_CORRECTION.get(species.lower(), {})
    corrected = {}
    for ep, score in scores.items():
        cf = corrections.get(ep, 1.0)
        corrected[ep] = round(min(score * cf, 1.0), 4)
    return corrected
