"""D2 Seed-Match Scorer (TargetScan-style, positions 2-8)"""
import re
from dataclasses import dataclass, field
from typing import List, Optional

COMPLEMENT = str.maketrans('AUGC','UACG')

def reverse_complement(seq: str) -> str:
    return seq.upper().translate(COMPLEMENT)[::-1]

@dataclass
class SeedMatchResult:
    query_seq: str; seed_region: str; seed_rc: str
    n_perfect_seed_hits: int; n_partial_seed_hits: int
    seed_match_score: float; top_hits: List[str] = field(default_factory=list)
    rnafold_available: bool = False; mfe_kcal: Optional[float] = None
    dsrna_probability: float = 0.0; d2_risk_contribution: float = 0.0

def score_seed_matches(query_seq: str, transcriptome_sample: Optional[List[str]] = None,
                       n_transcripts: int = 20000) -> SeedMatchResult:
    seq = query_seq.upper().replace('T','U')
    seed = seq[1:8]; seed_rc = reverse_complement(seed)
    seed_rc_6 = seed_rc[1:]
    perfect_hits = partial_hits = 0; top_hits = []
    if transcriptome_sample:
        for t in transcriptome_sample:
            t_up = t.upper().replace('T','U')
            if seed_rc in t_up: perfect_hits += 1; top_hits.append(t[:40])
            elif seed_rc_6 in t_up: partial_hits += 1
    else:
        gc = (seed_rc.count('G') + seed_rc.count('C')) / max(len(seed_rc),1)
        est = 0.01 + gc * 0.02
        perfect_hits = int(n_transcripts * est); partial_hits = int(n_transcripts * est * 2.5)
    seed_score = min((perfect_hits + partial_hits * 0.3) / (n_transcripts * 0.05), 1.0)
    # dsRNA via palindrome proxy (RNAfold if installed)
    mfe = None; dsrna_prob = 0.0; rnafold_avail = False
    try:
        import RNA
        fc = RNA.fold_compound(seq)
        structure, mfe = fc.mfe()
        paired = structure.count('(') + structure.count(')')
        dsrna_prob = paired / max(len(seq),1); rnafold_avail = True
    except Exception:
        L = len(seq)
        windows = [seq[i:i+6] for i in range(L-5)]
        pals = sum(1 for w in windows if len(w)==6 and w == reverse_complement(w))
        dsrna_prob = min(pals / max(len(windows),1) * 4.0, 1.0)
    d2_risk = min(seed_score * 0.70 + dsrna_prob * 0.30, 1.0)
    return SeedMatchResult(query_seq=query_seq, seed_region=seed, seed_rc=seed_rc,
        n_perfect_seed_hits=perfect_hits, n_partial_seed_hits=partial_hits,
        seed_match_score=seed_score, top_hits=top_hits[:5],
        rnafold_available=rnafold_avail, mfe_kcal=mfe,
        dsrna_probability=dsrna_prob, d2_risk_contribution=round(d2_risk,4))

def summary(r: SeedMatchResult) -> str:
    lines = ['  [SeedMatch] TargetScan-style D2 scorer',
             f'  Seed(2-8): {r.seed_region}  RC: {r.seed_rc}',
             f'  7-mer hits: {r.n_perfect_seed_hits}  6-mer: {r.n_partial_seed_hits}  '
             f'Seed score: {r.seed_match_score:.3f}',
             f'  dsRNA prob ({"RNAfold" if r.rnafold_available else "palindrome proxy"}): '
             f'{r.dsrna_probability:.3f}  D2 risk: {r.d2_risk_contribution:.3f}']
    if r.mfe_kcal: lines.append(f'  MFE: {r.mfe_kcal:.2f} kcal/mol')
    return '\n'.join(lines)
