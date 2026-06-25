from __future__ import annotations

from typing import Dict, List, Set


def _s(*genes: str) -> Set[str]:
    return {g.upper() for g in genes}


GENE_SETS: Dict[str, Set[str]] = {
    "interferon": _s("ISG15", "OAS1", "OAS2", "OAS3", "MX1", "MX2", "IFI44", "IFI44L", "IFI6", "IRF7", "STAT1", "CXCL10"),
    "nfkb": _s("NFKB1", "NFKB2", "RELA", "RELB", "TNFA", "IL6", "IL1B", "CXCL8", "CCL2"),
    "tlr_innate": _s("TLR3", "TLR7", "TLR8", "TLR9", "MYD88", "TRAF6", "IRF5", "IRF7"),
    "complement": _s("C1QA", "C1QB", "C1QC", "C2", "C3", "C4A", "C4B", "CFB", "CFD", "C5", "C6", "C7", "C8A", "C8B", "C9"),
    "coagulation": _s("F2", "F3", "F5", "F7", "F8", "F9", "F10", "F11", "F12", "SERPINE1", "THBD", "TFPI"),
    "platelet": _s("PF4", "PPBP", "ITGA2B", "ITGB3", "GP1BB", "VWF", "GP9", "THBS1"),
    "apoptosis": _s("BAX", "BAK1", "CASP3", "CASP8", "CASP9", "FAS", "FADD", "TNFSF10", "BCL2L11", "PMAIP1"),
    "p53_damage": _s("TP53", "CDKN1A", "GADD45A", "MDM2", "BBC3"),
    "mitochondrial_stress": _s("HMOX1", "SOD2", "UCP2", "PPARGC1A"),
    "oxidative_stress": _s("NQO1", "HMOX1", "SOD1", "SOD2", "TXNRD1", "GCLC", "GCLM"),
    "er_stress": _s("DDIT3", "ATF4", "XBP1", "HSPA5"),
    "lysosome": _s("LAMP1", "LAMP2", "CTSB", "CTSD", "CTSL", "GAA", "GBA"),
    "autophagy": _s("MAP1LC3B", "BECN1", "ATG5", "SQSTM1"),
    "kupffer_activation": _s("CD68", "CSF1R", "MARCO", "EMR1", "CLEC4F", "LYZ"),
    "inflammasome": _s("NLRP3", "PYCARD", "CASP1", "IL1B", "IL18"),
    "cytokine_signaling": _s("IL6", "IL1B", "TNFA", "CXCL10", "CCL2", "CCL5"),
    "chemokine_signaling": _s("CCL2", "CCL5", "CCL20", "CXCL8", "CXCL10", "CXCL9"),
    "jnk_p38_stress": _s("MAPK8", "MAPK9", "MAPK14", "FOS", "JUN"),
    "hepatocyte_injury": _s("ALDH1A1", "CYP3A4", "CYP2E1", "CYP2C9", "SLC10A1", "KRT8", "KRT18", "GLUL"),
    "cholestasis_bile_acid": _s("BSEP", "ABCB11", "SLC10A1", "ABCC2", "NR1H4"),
    "renal_tubular_stress": _s("HAVCR1", "LCN2", "SLC22A2", "SLC22A8"),
    "glomerular_injury": _s("NPHS1", "NPHS2", "PODXL", "COL4A3"),
    "endosomal_trafficking": _s("RAB5A", "RAB7A", "EEA1"),
    "rnaseh1_long_transcript_risk": _s("RNASEH1", "FOXP1", "DMD", "TTN"),
    "hybridization_offtarget": _s("DDIT3", "GADD45A", "ATF3"),
    "cell_cycle_arrest": _s("CDKN1A", "CDKN1B", "CDKN2A"),
    "dna_damage": _s("GADD45A", "GADD45B", "GADD45G", "RAD51"),
    "necrosis_damage": _s("HMGB1", "S100A8", "S100A9"),
    "fibrosis_remodeling": _s("COL1A1", "COL1A2", "COL3A1", "TGFB1", "ACTA2"),
    "kynurenine_inflammation": _s("IDO1", "TDO2"),
}
