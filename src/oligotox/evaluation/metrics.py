"""Clinical prediction evaluation metrics.

Implements AUROC, AUPRC, Brier score, calibration curves, and
leakage-aware evaluation protocols for oligonucleotide toxicity prediction.
"""
import numpy as np
import pandas as pd
from typing import Dict, Tuple, Optional
from pathlib import Path

try:
    from sklearn.metrics import roc_auc_score, average_precision_score, brier_score_loss
    from sklearn.calibration import calibration_curve
except ImportError:
    roc_auc_score = average_precision_score = brier_score_loss = calibration_curve = None


def compute_auroc(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    **kwargs
) -> float:
    """Compute Area Under ROC curve."""
    if roc_auc_score is None:
        raise ImportError("sklearn required for AUROC calculation")
    try:
        return roc_auc_score(y_true, y_pred, **kwargs)
    except ValueError:
        return np.nan


def compute_auprc(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    **kwargs
) -> float:
    """Compute Area Under Precision-Recall curve."""
    if average_precision_score is None:
        raise ImportError("sklearn required for AUPRC calculation")
    try:
        return average_precision_score(y_true, y_pred, **kwargs)
    except ValueError:
        return np.nan


def compute_brier_score(
    y_true: np.ndarray,
    y_pred: np.ndarray
) -> float:
    """Compute Brier score (calibration metric)."""
    if brier_score_loss is None:
        raise ImportError("sklearn required for Brier score calculation")
    try:
        return brier_score_loss(y_true, y_pred)
    except ValueError:
        return np.nan


def compute_calibration_curve(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    n_bins: int = 10,
    strategy: str = "uniform"
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute calibration curve (reliability diagram).
    
    Returns:
        prob_true: fraction of positives in each bin
        prob_pred: mean predicted probability in each bin
    """
    if calibration_curve is None:
        raise ImportError("sklearn required for calibration curve")
    try:
        return calibration_curve(y_true, y_pred, n_bins=n_bins, strategy=strategy)
    except ValueError:
        return np.array([]), np.array([])


def evaluate_endpoint(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    endpoint_name: str,
    threshold: Optional[float] = None
) -> Dict[str, float]:
    """Comprehensive evaluation for a single endpoint.
    
    Args:
        y_true: Binary labels
        y_pred: Predicted probabilities
        endpoint_name: Name of the toxicity endpoint
        threshold: Classification threshold (default 0.5)
    
    Returns:
        Dictionary with AUROC, AUPRC, Brier score, and threshold info
    """
    if threshold is None:
        threshold = 0.5
    
    result = {
        "endpoint": endpoint_name,
        "n_samples": len(y_true),
        "n_positive": int(y_true.sum()),
        "prevalence": float(y_true.mean()),
        "threshold": threshold,
    }
    
    try:
        result["auroc"] = compute_auroc(y_true, y_pred)
        result["auprc"] = compute_auprc(y_true, y_pred)
        result["brier"] = compute_brier_score(y_true, y_pred)
        
        # Binary classification metrics at threshold
        y_pred_binary = (y_pred >= threshold).astype(int)
        tp = ((y_true == 1) & (y_pred_binary == 1)).sum()
        tn = ((y_true == 0) & (y_pred_binary == 0)).sum()
        fp = ((y_true == 0) & (y_pred_binary == 1)).sum()
        fn = ((y_true == 1) & (y_pred_binary == 0)).sum()
        
        result["accuracy"] = float((tp + tn) / len(y_true))
        result["precision"] = float(tp / (tp + fp)) if (tp + fp) > 0 else np.nan
        result["recall"] = float(tp / (tp + fn)) if (tp + fn) > 0 else np.nan
        result["specificity"] = float(tn / (tn + fp)) if (tn + fp) > 0 else np.nan
        
    except Exception as e:
        result["error"] = str(e)
    
    return result


def generate_evaluation_report(
    y_true_dict: Dict[str, np.ndarray],
    y_pred_dict: Dict[str, np.ndarray],
    output_path: Optional[Path] = None
) -> Dict:
    """Generate comprehensive evaluation report across all endpoints.
    
    Args:
        y_true_dict: {endpoint: binary labels}
        y_pred_dict: {endpoint: predicted probabilities}
        output_path: Optional path to save JSON report
    
    Returns:
        Complete evaluation report dictionary
    """
    report = {
        "n_endpoints": len(y_true_dict),
        "per_endpoint": {},
    }
    
    auroc_list = []
    brier_list = []
    
    for endpoint in y_true_dict.keys():
        if endpoint not in y_pred_dict:
            continue
        
        y_t = y_true_dict[endpoint]
        y_p = y_pred_dict[endpoint]
        
        # Filter out NaNs
        valid_mask = ~(np.isnan(y_t) | np.isnan(y_p))
        y_t = y_t[valid_mask]
        y_p = y_p[valid_mask]
        
        if len(y_t) < 10:  # Min samples for meaningful metrics
            continue
        
        endpoint_result = evaluate_endpoint(y_t, y_p, endpoint)
        report["per_endpoint"][endpoint] = endpoint_result
        
        if not np.isnan(endpoint_result.get("auroc", np.nan)):
            auroc_list.append(endpoint_result["auroc"])
        if not np.isnan(endpoint_result.get("brier", np.nan)):
            brier_list.append(endpoint_result["brier"])
    
    # Composite metrics
    if auroc_list:
        report["composite_auroc"] = float(np.mean(auroc_list))
    if brier_list:
        report["composite_brier"] = float(np.mean(brier_list))
    
    if output_path:
        import json
        with open(output_path, "w") as f:
            json.dump(report, f, indent=2)
    
    return report
