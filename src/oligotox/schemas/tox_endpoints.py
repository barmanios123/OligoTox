"""OligoTox Endpoint Registry"""
from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, Optional, List

class RiskDomain(str, Enum):
    D1 = 'D1'; D2 = 'D2'; D3 = 'D3'

class RiskClass(str, Enum):
    LOW = 'low'; MODERATE = 'moderate'; HIGH = 'high'

class ToxEndpoint(str, Enum):
    HEPATOCELLULAR_INJURY    = 'hepatocellular_injury'
    HEPATIC_VACUOLATION      = 'hepatic_vacuolation'
    RENAL_TUBULAR_INJURY     = 'renal_tubular_injury'
    RENAL_ACCUMULATION       = 'renal_accumulation'
    THROMBOCYTOPENIA         = 'thrombocytopenia'
    COAGULATION_EFFECT       = 'coagulation_effect'
    COMPLEMENT_ACTIVATION    = 'complement_activation'
    INNATE_IMMUNE_ACTIVATION = 'innate_immune_activation'
    CYTOKINE_RELEASE         = 'cytokine_release'
    KUPFFER_CELL_ACTIVATION  = 'kupffer_cell_activation'
    CNS_TOXICITY             = 'cns_toxicity'
    OFFTARGET_GENE_SILENCING = 'offtarget_gene_silencing'
    INJECTION_SITE_REACTION  = 'injection_site_reaction'

ENDPOINT_DOMAIN: Dict['ToxEndpoint', RiskDomain] = {
    ToxEndpoint.CNS_TOXICITY:             RiskDomain.D1,
    ToxEndpoint.OFFTARGET_GENE_SILENCING: RiskDomain.D2,
    ToxEndpoint.HEPATOCELLULAR_INJURY:    RiskDomain.D3,
    ToxEndpoint.HEPATIC_VACUOLATION:      RiskDomain.D3,
    ToxEndpoint.RENAL_TUBULAR_INJURY:     RiskDomain.D3,
    ToxEndpoint.RENAL_ACCUMULATION:       RiskDomain.D3,
    ToxEndpoint.THROMBOCYTOPENIA:         RiskDomain.D3,
    ToxEndpoint.COAGULATION_EFFECT:       RiskDomain.D3,
    ToxEndpoint.COMPLEMENT_ACTIVATION:    RiskDomain.D3,
    ToxEndpoint.INNATE_IMMUNE_ACTIVATION: RiskDomain.D3,
    ToxEndpoint.CYTOKINE_RELEASE:         RiskDomain.D3,
    ToxEndpoint.KUPFFER_CELL_ACTIVATION:  RiskDomain.D3,
    ToxEndpoint.INJECTION_SITE_REACTION:  RiskDomain.D3,
}

@dataclass
class ToxPrediction:
    endpoint: ToxEndpoint
    risk_score: float
    risk_class: RiskClass
    confidence: float
    primary_domain: RiskDomain
    top_features: List[str]          = field(default_factory=list)
    mechanistic_note: Optional[str]  = None
    assay_recommendation: Optional[str] = None

@dataclass
class ToxProfile:
    sample_id: str
    predictions: Dict[ToxEndpoint, ToxPrediction] = field(default_factory=dict)
    flags: List[str]                               = field(default_factory=list)

    @property
    def composite_score(self) -> float:
        if not self.predictions: return 0.0
        return sum(p.risk_score for p in self.predictions.values()) / len(self.predictions)

    def summary_table(self) -> str:
        if not self.predictions: return '  No predictions.'
        sorted_preds = sorted(self.predictions.values(), key=lambda p: -p.risk_score)
        composite = self.composite_score
        dom_counts = {d: 0 for d in RiskDomain}
        for p in sorted_preds: dom_counts[p.primary_domain] += 1
        dominant_ep = sorted_preds[0].endpoint.value if sorted_preds else 'n/a'
        d1_risk = max((p.risk_score for p in sorted_preds if p.primary_domain == RiskDomain.D1), default=0.0)
        d2_risk = max((p.risk_score for p in sorted_preds if p.primary_domain == RiskDomain.D2), default=0.0)
        d3_risk = max((p.risk_score for p in sorted_preds if p.primary_domain == RiskDomain.D3), default=0.0)
        def rc(s): return 'high' if s >= 0.60 else ('moderate' if s >= 0.30 else 'low')
        lines = [
            f'  Composite : {composite:.3f}  (low<0.30 / moderate 0.30-0.60 / high>=0.60)',
            f'  Dominant  : {dominant_ep}',
            f'  D1 : risk={d1_risk:.2f} ({rc(d1_risk)})  '
            f'D2 : risk={d2_risk:.2f} ({rc(d2_risk)})  '
            f'D3 : risk={d3_risk:.2f} ({rc(d3_risk)})',
            f'  {"Endpoint":<42}{"Dom":>5}{"Score":>8}  {"Class":<12}{"Conf":>6}',
            '  ' + '-' * 76,
        ]
        for p in sorted_preds:
            lines.append(f'  {p.endpoint.value:<42}{p.primary_domain.value:>5}{p.risk_score:>8.3f}'
                         f'  {p.risk_class.value:<12}{p.confidence:>6.2f}')
        for flag in self.flags:
            lines.append(f'\n  FLAG: {flag}')
        return '\n'.join(lines)
