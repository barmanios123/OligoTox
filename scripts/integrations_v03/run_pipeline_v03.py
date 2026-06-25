#!/usr/bin/env python3
"""
OligoTox v0.3 Integration Pipeline
Usage: PYTHONPATH=. python scripts/run_pipeline.py [--train] [--split scaffold|group|datasail]
"""
import argparse, sys, os, logging
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

from oligotox.integrations_v03.feature_matrix import build_feature_matrix, ENDPOINTS
from oligotox.integrations_v03.datasail_splits import scaffold_sequence_split, datasail_split, group_shuffle_split, check_leakage, evaluate_endpoint

DEMO_CANDIDATES = [
    {"sequence":"ACGTACGTACGTACGTACGT","modality":"gapmer_ASO","ps_fraction":1.0,
     "lna_positions":[0,1,2,17,18,19],"delivery":"naked","route":"SC",
     "hepatocellular_injury":1,"offtarget_gene_silencing":1},
    {"sequence":"GCGCGCGCGCGCGCGCGCGC","modality":"siRNA","ps_fraction":0.2,
     "delivery":"LNP","route":"IV","complement_activation":1,"lnp_immune_effect":1},
    {"sequence":"UUUUACGUACGUACGUACGU","modality":"siRNA","ps_fraction":0.3,
     "delivery":"GalNAc","route":"SC","innate_immune_activation":1},
    {"sequence":"TTTTCGTTCGTTCGTTCGTT","modality":"gapmer_ASO","ps_fraction":0.9,
     "delivery":"naked","route":"SC","thrombocytopenia":1},
    {"sequence":"AAAAACGUACGUACGUACGU","modality":"siRNA","ps_fraction":0.2,
     "delivery":"LNP","route":"IV","hepatocellular_injury":0},
    {"sequence":"CCCCGCGCGCGCGCGCGCGC","modality":"siRNA","ps_fraction":0.3,
     "delivery":"GalNAc","route":"SC","offtarget_gene_silencing":0},
]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input",  default=None)
    ap.add_argument("--output", default="data/features/features_v03.csv")
    ap.add_argument("--split",  default="scaffold", choices=["scaffold","group","datasail"])
    ap.add_argument("--train",  action="store_true")
    ap.add_argument("--no-viennarna", action="store_true")
    ap.add_argument("--embeddings",   action="store_true")
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    # Load candidates
    if args.input and os.path.exists(args.input):
        import pandas as pd
        df_raw = pd.read_csv(args.input)
        candidates = df_raw.to_dict(orient="records")
    else:
        logger.info("No input CSV — using demo candidates.")
        import pandas as pd
        candidates = DEMO_CANDIDATES
        df_raw = pd.DataFrame(candidates)

    # Build features
    logger.info("Building integrated feature matrix...")
    df_feats = build_feature_matrix(
        candidates, use_viennarna=not args.no_viennarna,
        use_embeddings=args.embeddings)
    label_cols = [ep for ep in ENDPOINTS if ep in df_raw.columns]
    for col in label_cols:
        df_feats[col] = df_raw[col].fillna(0).values

    df_feats.to_csv(args.output, index=False)
    feat_cols = [c for c in df_feats.columns if c not in ENDPOINTS]
    logger.info(f"Saved {df_feats.shape[0]} rows × {len(feat_cols)} features → {args.output}")

    # Leakage-free splits + training
    if args.train and label_cols:
        seqs = [c.get("sequence","") for c in candidates]
        X = df_feats[feat_cols].fillna(0).values
        y_all = {ep: df_feats[ep].values for ep in label_cols}
        y_dummy = y_all[label_cols[0]]

        if args.split == "datasail":
            tr, va, te = datasail_split(seqs, y_dummy)
        elif args.split == "group":
            groups = df_raw.get("target_gene", df_raw.get("drug_name",
                     [f"g{i}" for i in range(len(candidates))])).tolist()
            tr, va, te = group_shuffle_split(X, y_dummy, groups)
        else:
            tr, va, te = scaffold_sequence_split(seqs, y_dummy)

        report = check_leakage([seqs[i] for i in tr], [seqs[i] for i in te])
        logger.info(f"Leakage report: {report}")

        from sklearn.ensemble import RandomForestClassifier
        from sklearn.calibration import CalibratedClassifierCV
        for ep in label_cols:
            y_tr = y_all[ep][tr]
            if len(np.unique(y_tr)) < 2:
                logger.warning(f"  {ep}: single class — skip"); continue
            clf = CalibratedClassifierCV(
                RandomForestClassifier(100, class_weight="balanced", random_state=42),
                cv=min(3, int(y_tr.sum()) or 3), method="isotonic")
            clf.fit(X[tr], y_tr)
            y_va = y_all[ep][va]
            if len(np.unique(y_va)) >= 2:
                r = evaluate_endpoint(y_va, clf.predict_proba(X[va])[:,1], ep)
                logger.info(f"  {ep}: PR-AUC={r['pr_auc']:.3f} ROC-AUC={r['roc_auc']:.3f} Brier={r['brier']:.3f}")
    logger.info("Pipeline complete.")

if __name__ == "__main__":
    main()
