"""
Integration 2: ASO Chemistry Features (ASOptimizer-inspired)
=============================================================
Extracts 14 chemistry features from sequence + modification annotation.

Local feature extraction: always available.
ASOptimizer API: optional (set use_api=True + provide API key/URL).

Ref: Kang et al., NAR 2025; ASOptimizer: DOI:10.1093/nar/gkaf349
     OligoAI: Liang et al., 2022
"""
from __future__ import annotations
import ast
import math
import re
import math
import logging
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)


def _coerce_numeric_positions(x):
    if x is None:
        return []
    try:
        if isinstance(x, float) and math.isnan(x):
            return []
    except Exception:
        pass

    if isinstance(x, list):
        vals = x
    elif isinstance(x, tuple):
        vals = list(x)
    elif isinstance(x, str):
        s = x.strip()
        if not s:
            return []
        try:
            parsed = ast.literal_eval(s)
            if isinstance(parsed, list):
                vals = parsed
            elif isinstance(parsed, tuple):
                vals = list(parsed)
            else:
                vals = [parsed]
        except Exception:
            if "," in s:
                vals = [p.strip() for p in s.split(",") if p.strip()]
            else:
                vals = [s]
    else:
        vals = [x]

    out = []
    for v in vals:
        if v is None:
            continue
        try:
            if isinstance(v, float) and math.isnan(v):
                continue
        except Exception:
            pass
        if isinstance(v, (int, float)):
            out.append(int(v))
            continue
        sv = str(v).strip()
        if not sv:
            continue
        try:
            out.append(int(float(sv)))
        except Exception:
            continue
    return out



def extract_chemistry_features(
    sequence: str,
    modification_string: Optional[str] = None,
    ps_fraction: Optional[float] = None,
    lna_positions: Optional[List[int]] = None,
    sugar_mods: Optional[Dict[str, float]] = None,
) -> Dict[str, float]:
    """
    Extract ASO chemistry features from sequence + modification metadata.

    modification_string format (Exiqon/IDT style):
      +N = LNA nucleotide
      mN = 2'-OMe
      dN = DNA (unmodified)
      fN = 2'-F
      * = phosphorothioate linkage (between nucleotides)
    """
    seq = sequence.upper().replace("U", "T")
    n = len(seq)
    feats: Dict[str, float] = {}

    # --- PS fraction ---
    if ps_fraction is not None:
        feats["chem_ps_fraction"] = float(ps_fraction)
    else:
        # Infer from modification string
        if modification_string:
            ps_count = modification_string.count("*")
            feats["chem_ps_fraction"] = ps_count / max(n - 1, 1)
        else:
            feats["chem_ps_fraction"] = 1.0  # conservative default for ASOs

    # PS clustering: consecutive PS linkages (PS patches amplify protein binding)
    ps_frac = feats["chem_ps_fraction"]
    feats["chem_ps_cluster_score"] = ps_frac ** 2  # non-linear risk
    feats["chem_full_ps"] = float(ps_frac >= 0.95)

    # --- LNA features ---
    if lna_positions is not None:
        lna_pos = _coerce_numeric_positions(lna_positions)
    elif modification_string:
        lna_pos = [i for i, c in enumerate(modification_string) if c == "+"]
    else:
        lna_pos = []

    feats["chem_lna_count"] = float(len(lna_pos))
    feats["chem_lna_fraction"] = len(lna_pos) / max(n, 1)

    # LNA flank asymmetry: 5' vs 3' wing ratio (Papargyri finding)
    if lna_pos:
        mid = n // 2
        lna_5p = sum(1 for p in lna_pos if p < mid)
        lna_3p = sum(1 for p in lna_pos if p >= mid)
        feats["chem_lna_flank_asym"] = abs(lna_5p - lna_3p) / max(len(lna_pos), 1)
    else:
        feats["chem_lna_flank_asym"] = 0.0

    # High-LNA gapmer flag (≥3 LNA flanks)
    feats["chem_high_lna_gapmer_flag"] = float(len(lna_pos) >= 6)

    # --- Sugar mods ---
    sugar_mods = sugar_mods or {}
    feats["chem_moe_fraction"]  = float(sugar_mods.get("2-MOE", 0.0))
    feats["chem_ome_fraction"]  = float(sugar_mods.get("2-OMe", 0.0))
    feats["chem_2f_fraction"]   = float(sugar_mods.get("2-F",   0.0))

    # --- Gapmer structure (from modification string) ---
    gap_len, wing_len = _parse_gapmer_structure(modification_string, n)
    feats["chem_gap_length"]  = float(gap_len)
    feats["chem_wing_length"] = float(wing_len)

    # --- Activity proxy (inverted = off-target burden) ---
    # Empirical: high-LNA + full PS + short gap → high activity → higher off-target risk
    activity_proxy = (
        feats["chem_lna_fraction"] * 0.4 +
        feats["chem_ps_fraction"]  * 0.3 +
        (1 / max(gap_len, 1)) * 0.3
    )
    feats["chem_activity_proxy"] = min(activity_proxy, 1.0)
    feats["chem_offtarget_activity_proxy"] = feats["chem_activity_proxy"] * feats["chem_ps_fraction"]

    return feats


def _parse_gapmer_structure(mod_string: Optional[str], seq_len: int):
    """Parse gap and wing lengths from modification string."""
    if not mod_string:
        return max(seq_len - 6, 5), 3  # typical 3-10-3 defaults
    # Count runs of unmodified DNA bases (gap) vs modified bases (wings)
    has_lna = "+" in mod_string
    has_moe = "m" in mod_string.lower()
    if not has_lna and not has_moe:
        return seq_len, 0
    # Simple heuristic: gap is the central unmodified run
    positions = list(mod_string)
    unmod = sum(1 for c in positions if c in "ACGT")
    mod   = sum(1 for c in positions if c in "+m")
    gap_len  = unmod
    wing_len = mod // 2 if mod > 0 else 0
    return max(gap_len, 1), max(wing_len, 0)


def batch_aso_features(
    sequences: List[str],
    mods: Optional[List[Optional[str]]] = None,
    ps_fractions: Optional[List[Optional[float]]] = None,
) -> List[Dict[str, float]]:
    """Batch extract ASO chemistry features for a list of sequences."""
    mods = mods or [None] * len(sequences)
    ps_fractions = ps_fractions or [None] * len(sequences)
    return [
        extract_chemistry_features(seq, mod, psf)
        for seq, mod, psf in zip(sequences, mods, ps_fractions)
    ]


if __name__ == "__main__":
    test_cases = [
        ("ACGTACGTACGTACGTACGT", "+A+C+GGTTACATGmAmAmAmA", 1.0, [0,1,2,17,18,19]),
        ("GCGCGCGCGCGCGCGCGCGC", None, 0.2, None),
    ]
    for seq, mod, ps, lna in test_cases:
        feats = extract_chemistry_features(seq, mod, ps, lna)
        print(f"\nSequence: {seq[:15]}...")
        for k, v in feats.items():
            print(f"  {k}: {v:.3f}")
