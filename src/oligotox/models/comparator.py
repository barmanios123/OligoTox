"""Candidate Comparator"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional
from oligotox.schemas.tox_endpoints import ToxProfile, ToxEndpoint

@dataclass
class CandidateResult:
    sample_id: str; composite_score: float; rank: int
    endpoint_scores: Dict[str,float] = field(default_factory=dict)
    dominant_endpoint: str = ''

def compare_candidates(profiles: List[ToxProfile]) -> List[CandidateResult]:
    results = []
    for p in profiles:
        if not p.predictions: continue
        scores = {ep.value: pred.risk_score for ep, pred in p.predictions.items()}
        dominant = max(scores, key=scores.get, default='')
        results.append(CandidateResult(sample_id=p.sample_id, composite_score=p.composite_score,
                                        rank=0, endpoint_scores=scores, dominant_endpoint=dominant))
    results.sort(key=lambda x: x.composite_score)
    for i, r in enumerate(results, 1): r.rank = i
    return results
