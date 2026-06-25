import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score

DATA_PATH = "data/processed/training_dataset.csv"
TARGET = "hepatocellular_injury"  # as an example

def main():
    df = pd.read_csv(DATA_PATH)
    df = df.copy()

    features = [c for c in df.columns
                if c not in ["candidate_id","source","notes","series_id",
                             "hepatocellular_injury","hepatic_vacuolation",
                             "renal_tubular_injury","renal_accumulation",
                             "thrombocytopenia","coagulation_effect",
                             "complement_activation","innate_immune_activation",
                             "cytokine_release","kupffer_cell_activation",
                             "cns_toxicity","offtarget_gene_silencing",
                             "injection_site_reaction"]]

    series = df["series_id"].unique()
    print("Series:", series)

    for s in series:
        train = df[df["series_id"] != s]
        test  = df[df["series_id"] == s]

        if train[TARGET].sum() == 0 or test[TARGET].sum() == 0:
            print(f"Skipping series {s}: no positives in train or test for {TARGET}")
            continue

        X_train = train[features]
        y_train = train[TARGET]
        X_test  = test[features]
        y_test  = test[TARGET]

        clf = RandomForestClassifier(
            n_estimators=200,
            random_state=42,
            class_weight="balanced"
        )
        clf.fit(X_train, y_train)
        probs = clf.predict_proba(X_test)[:,1]
        auc = roc_auc_score(y_test, probs)
        print(f"Series {s}: n_train={len(train)}, n_test={len(test)}, AUC={auc:.3f}")

if __name__ == "__main__":
    main()
