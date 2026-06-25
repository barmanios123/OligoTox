"""Feature matrix builder — assembles all 5 integrations into one row."""
from __future__ import annotations
import numpy as np
import logging
from typing import Dict, List, Optional
import math, ast

logger = logging.getLogger(__name__)


def _is_missing(x):
    try:
        return x is None or (isinstance(x, float) and math.isnan(x))
    except Exception:
        return False

def _normalize_listlike(x):
    if _is_missing(x) or x == "":
        return []
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
        return []

    out = []
    for v in vals:
        if _is_missing(v) or v == "":
            continue
        if isinstance(v, (int, float)):
            out.append(int(v))
            continue
        sv = str(v).strip()
        if not sv:
            continue
        try:
            if "." in sv:
                fv = float(sv)
                if fv.is_integer():
                    out.append(int(fv))
                else:
                    out.append(fv)
            else:
                out.append(int(sv))
        except Exception:
            out.append(sv)
    return out
    if isinstance(x, list):
        return x
    if isinstance(x, tuple):
        return list(x)
    if isinstance(x, str):
        s = x.strip()
        if not s:
            return []
        try:
            parsed = ast.literal_eval(s)
            if isinstance(parsed, list):
                return parsed
            if isinstance(parsed, tuple):
                return list(parsed)
        except Exception:
            pass
        if "," in s:
            return [p.strip() for p in s.split(",") if p.strip()]
        return [s]
    return []


ENDPOINTS = [
    "hepatocellular_injury","renal_tubular_injury","thrombocytopenia",
    "coagulation_effect","complement_activation","innate_immune_activation",
    "cytokine_release","cns_neurotoxicity","reproductive_developmental",
    "offtarget_gene_silencing","injection_site_reaction",
    "lnp_immune_effect","proinflammatory_vasculitis",
]


def build_feature_row(
    candidate: dict,
    use_embeddings: bool = False,
    use_viennarna: bool = True,
    embedding_model: str = "nucleotide_transformer",
    embedding_device: str = "cpu",
) -> Dict[str, float]:
    from oligotox.integrations_v03.aso_optimizer import extract_chemistry_features
    from oligotox.integrations_v03.offtarget_burden import compute_transcriptome_burden
    from oligotox.integrations_v03.structure_signals import compute_structure_features

    seq = candidate.get("sequence", "")
    feats: Dict[str, float] = {}

    # Integration 2: chemistry
    feats.update(extract_chemistry_features(
        seq,
        modification_string=None if _is_missing(candidate.get("modification_string")) else candidate.get("modification_string"),
        ps_fraction=0.0 if _is_missing(candidate.get("ps_fraction")) else candidate.get("ps_fraction"),
        lna_positions=_normalize_listlike(candidate.get("lna_positions")),
        sugar_mods=_normalize_listlike(candidate.get("sugar_mods")),
    ))

    # Integration 3: off-target burden
    ot = compute_transcriptome_burden(seq, use_mock_transcriptome=True)
    ot.pop("d2_seed_seq", None)
    feats.update(ot)

    # Integration 4: structure/innate
    feats.update(compute_structure_features(seq, use_viennarna=use_viennarna))

    # Integration 1: embeddings (optional — GPU recommended)
    if use_embeddings:
        try:
            from oligotox.integrations_v03.sequence_embeddings import SequenceEmbedder
            emb = SequenceEmbedder(model=embedding_model, device=embedding_device)
            cols = emb.embed_sequences([seq])
            feats.update({k: v[0] for k, v in cols.items()})
        except Exception as e:
            logger.warning(f"Embedding skipped: {e}")

    # Delivery/route context flags
    feats["ctx_delivery_lnp"]    = float(candidate.get("delivery","") == "LNP")
    feats["ctx_delivery_galnac"] = float(candidate.get("delivery","") == "GalNAc")
    feats["ctx_delivery_naked"]  = float(candidate.get("delivery","") == "naked")
    feats["ctx_route_iv"]  = float(candidate.get("route","") == "IV")
    feats["ctx_route_sc"]  = float(candidate.get("route","") == "SC")
    feats["ctx_route_it"]  = float(candidate.get("route","") == "IT")

    return feats


def build_feature_matrix(candidates: List[dict], **kwargs):
    import pandas as pd
    rows = [build_feature_row(c, **kwargs) for c in candidates]
    return pd.DataFrame(rows).fillna(0.0)


if __name__ == "__main__":
    import pandas as pd
    cs = [
        {"sequence":"ACGTACGTACGTACGTACGT","modality":"gapmer_ASO",
         "ps_fraction":1.0,"lna_positions":[0,1,2,17,18,19],"delivery":"naked","route":"SC"},
        {"sequence":"GCGCGCGCGCGCGCGCGCGC","modality":"siRNA",
         "ps_fraction":0.2,"delivery":"LNP","route":"IV"},
    ]
    df = build_feature_matrix(cs, use_viennarna=True, use_embeddings=False)
    print(f"Feature matrix: {df.shape[0]} rows × {df.shape[1]} features")
    print(df.T.to_string())
