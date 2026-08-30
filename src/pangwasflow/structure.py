from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import chi2


def pca_scores(features: pd.DataFrame, n_components: int = 2) -> tuple[pd.DataFrame, np.ndarray]:
    """Estimate population-structure covariates from a binary feature matrix."""
    if "sample_id" not in features.columns:
        raise ValueError("features must contain sample_id")
    matrix = features.drop(columns=["sample_id"]).to_numpy(dtype=float)
    scale = matrix.std(axis=0, ddof=0)
    keep = scale > 0
    if not np.any(keep):
        raise ValueError("no variable features available for PCA")
    matrix = matrix[:, keep]
    matrix = (matrix - matrix.mean(axis=0)) / matrix.std(axis=0, ddof=0)
    u, s, vt = np.linalg.svd(matrix, full_matrices=False)
    n_components = min(int(n_components), u.shape[1])
    for index in range(n_components):
        anchor = int(np.argmax(np.abs(vt[index])))
        if vt[index, anchor] < 0:
            u[:, index] *= -1
            vt[index] *= -1
    scores = u[:, :n_components] * s[:n_components]
    total = float(np.sum(s**2))
    variance = (s[:n_components] ** 2) / total if total else np.zeros(n_components)
    result = pd.DataFrame({"sample_id": features["sample_id"].astype(str)})
    for index in range(n_components):
        result[f"PC{index + 1}"] = scores[:, index]
    return result, variance


def genomic_inflation(pvalues: np.ndarray) -> float:
    """Return genomic-inflation lambda from one-degree-of-freedom p-values."""
    p = np.asarray(pvalues, dtype=float)
    p = p[np.isfinite(p)]
    if p.size == 0:
        return float("nan")
    p = np.clip(p, np.finfo(float).tiny, 1.0)
    statistics = chi2.isf(p, df=1)
    return float(np.median(statistics) / chi2.ppf(0.5, df=1))
