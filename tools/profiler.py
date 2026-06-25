"""OligoTox Candidate Profiler"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from typing import Optional
from oligotox.schemas.data_models import OligoInput, DeliveryInput, PlateletComplement, OmicsReadout, OligoToxInput
from oligotox.features.engineer import extract_features
from oligotox.rules.engine import MechanisticRulesEngine
from oligotox.schemas.tox_endpoints import ToxProfile, ToxEndpoint, ToxPrediction, RiskClass, ENDPOINT_DOMAIN, RiskDomain

_eng = MechanisticRulesEngine()

class OligoToxProfiler:
    def profile(self, sample_id: str, oligo: OligoInput, delivery: DeliveryInput,
                platelet: Optional[PlateletComplement] = None,
                omics: Optional[OmicsReadout] = None) -> ToxProfile:
        if platelet is None: platelet = PlateletComplement()
        if omics is None: omics = OmicsReadout()
        inp = OligoToxInput(oligo=oligo, delivery=delivery, platelet_complement=platelet, omics=omics)
        feats = extract_features(inp)
        preds_raw = _eng.score_all(feats, inp)
        preds = {ToxEndpoint(ep): pred for ep, pred in preds_raw.items()}
        flags = []
        if not oligo.nucleotides:
            flags.append('No per-nucleotide chemistry; LNA/PS estimated from modality defaults.')
        if not omics.deg_list:
            flags.append('No omics/DEG list; confidence on cytokine/offtarget endpoints reduced.')
        return ToxProfile(sample_id=sample_id, predictions=preds, flags=flags)
