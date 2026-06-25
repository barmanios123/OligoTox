"""OligoTox Data Models — typed input schemas"""
from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Optional


class Modality(str, Enum):
    SIRNA_LNP        = "siRNA_LNP"
    SIRNA_GALNAC     = "siRNA_GalNAc"
    SIRNA_GENERIC    = "siRNA"
    ASO_GAPMER       = "ASO_gapmer"
    ASO_SSO          = "ASO_splice"
    PMO              = "PMO"
    APTAMER          = "aptamer"
    ANTI_MIR         = "anti-miR"
    MIRNA_MIMIC      = "miRNA_mimic"
    MIRNA_INHIBITOR  = "miRNA_inhibitor"
    LNCRNA_TARGETING = "lncRNA_targeting"
    CIRCRNA          = "circRNA"
    OTHER            = "other"


class DeliveryVehicle(str, Enum):
    NAKED     = "none"
    GALNAC    = "GalNAc"
    LNP       = "LNP"
    CPP       = "CPP"
    SALINE    = "saline"
    LIPOSOMAL = "liposomal"
    OTHER     = "other"


class Route(str, Enum):
    SC    = "SC"
    IV    = "IV"
    IT    = "IT"
    PO    = "PO"
    IM    = "IM"
    ICV   = "ICV"
    OTHER = "other"


class Species(str, Enum):
    RAT        = "rat"
    MOUSE      = "mouse"
    CYNOMOLGUS = "cynomolgus"
    HUMAN      = "human"


class TargetClass(str, Enum):
    CODING_RNA    = "coding_RNA"
    LNCRNA        = "lncRNA"
    MIRNA_NETWORK = "miRNA_network"
    CIRCRNA       = "circRNA"


class DevelopmentStage(str, Enum):
    CANDIDATE_NOM = "candidate_nomination"
    PRE_IND       = "pre_IND"
    GLP_REPEAT    = "GLP_repeat_dose"
    RECOVERY      = "recovery"


class BackboneModification(str, Enum):
    PS     = "PS"
    PO     = "PO"
    PMO_BB = "PMO"
    PNA    = "PNA"


class SugarModification(str, Enum):
    NONE  = "none"
    MOE_2 = "2-MOE"
    OME_2 = "2-OMe"
    F_2   = "2-F"
    LNA   = "LNA"
    DNA   = "DNA"


@dataclass
class NucleotideUnit:
    position:     int
    base:         str
    sugar_mod:    str = "none"
    backbone_mod: str = "PO"


@dataclass
class OligoInput:
    sequence:                       str
    modality:                       Modality
    nucleotides:                    List[NucleotideUnit] = field(default_factory=list)
    target_gene:                    Optional[str]        = None
    target_class:                   TargetClass          = TargetClass.CODING_RNA
    conjugate_5prime:               Optional[str]        = None
    conjugate_3prime:               Optional[str]        = None
    is_human_specific_sequence:     bool                 = False
    target_is_lncrna_hub:           bool                 = False
    micropeptide_co_depletion_risk: bool                 = False


@dataclass
class DeliveryInput:
    vehicle:           DeliveryVehicle
    route:             Route           = Route.SC
    dose_mg_kg:        float           = 1.0
    repeat_dose:       bool            = False
    ionizable_lipid:   Optional[str]   = None
    peg_lipid_mol_pct: Optional[float] = None
    particle_size_nm:  Optional[float] = None
    pdi:               Optional[float] = None
    lipid_np_ratio:    Optional[float] = None
    galnac_valency:    Optional[int]   = None


@dataclass
class OffTargetContext:
    mirna_seed_overlap:          bool  = False
    n_predicted_offtargets:      int   = 0
    transcriptome_perturb_score: float = 0.0
    ncrna_complexity_flag:       bool  = False


@dataclass
class PlateletComplement:
    platelet_count_baseline:       float = 220.0
    platelet_count_treated:        float = 220.0
    pf4_ps_elisa_od:               float = 0.0
    c3a_pg_ml:                     float = 0.0
    c5a_pg_ml:                     float = 0.0
    sc5b9_ng_ml:                   float = 0.0
    complement_activation_ex_vivo: bool  = False
    platelet_c3d_c4d_positive:     bool  = False


@dataclass
class OmicsReadout:
    alt_iu_l:         float            = 20.0
    ast_iu_l:         float            = 20.0
    bun_mg_dl:        float            = 15.0
    creatinine_mg_dl: float            = 0.7
    cytokine_il6:     float            = 0.0
    cytokine_tnfa:    float            = 0.0
    cytokine_ifna:    float            = 0.0
    pathway_scores:   Dict[str, float] = field(default_factory=dict)
    deg_list:         List[str]        = field(default_factory=list)


@dataclass
class OligoToxInput:
    oligo:               OligoInput
    delivery:            DeliveryInput
    off_target:          OffTargetContext    = field(default_factory=OffTargetContext)
    platelet_complement: PlateletComplement  = field(default_factory=PlateletComplement)
    omics:               OmicsReadout        = field(default_factory=OmicsReadout)
    species:             Species             = Species.CYNOMOLGUS
    dev_stage:           DevelopmentStage    = DevelopmentStage.CANDIDATE_NOM
    sample_id:           Optional[str]       = None
    program_type:        Optional[str]       = None
