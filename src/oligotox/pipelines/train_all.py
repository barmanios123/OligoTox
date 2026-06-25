from oligotox.models.trainer import train_all


def run(feature_csv, label_csv, model_dir, metrics_dir):
    return train_all(
        feature_csv=feature_csv,
        label_csv=label_csv,
        model_dir=model_dir,
        metrics_dir=metrics_dir,
    )