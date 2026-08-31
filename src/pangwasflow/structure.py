from __future__ import annotations

from pathlib import Path

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


def read_distance_matrix(path: str | Path) -> pd.DataFrame:
    """Read and validate a square sample-by-sample distance matrix."""
    matrix = pd.read_csv(path, sep="\t", index_col=0)
    matrix.index = matrix.index.astype(str)
    matrix.columns = matrix.columns.astype(str)
    if matrix.index.duplicated().any() or matrix.columns.duplicated().any():
        raise ValueError("distance matrix contains duplicate sample identifiers")
    if matrix.shape[0] != matrix.shape[1]:
        raise ValueError("distance matrix must be square")
    if set(matrix.index) != set(matrix.columns):
        raise ValueError("distance matrix row and column sample identifiers must match")
    matrix = matrix.loc[matrix.index, matrix.index]
    try:
        values = matrix.to_numpy(dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError("distance matrix must contain numeric values") from exc
    if not np.isfinite(values).all():
        raise ValueError("distance matrix contains missing or non-finite values")
    if np.any(values < -1e-12):
        raise ValueError("distance matrix cannot contain negative distances")
    if not np.allclose(values, values.T, atol=1e-8, rtol=1e-6):
        raise ValueError("distance matrix must be symmetric")
    if not np.allclose(np.diag(values), 0.0, atol=1e-8):
        raise ValueError("distance matrix diagonal must be zero")
    return matrix.astype(float)


def classical_mds_scores(
    distances: pd.DataFrame,
    sample_ids: list[str] | pd.Series | np.ndarray | None = None,
    n_components: int = 2,
) -> tuple[pd.DataFrame, np.ndarray]:
    """Compute classical multidimensional scaling coordinates from genomic distances."""
    if distances.shape[0] != distances.shape[1]:
        raise ValueError("distances must be a square matrix")
    if sample_ids is None:
        ordered_samples = distances.index.astype(str).tolist()
    else:
        ordered_samples = [str(sample) for sample in sample_ids]
        if len(set(ordered_samples)) != len(ordered_samples):
            raise ValueError("requested MDS sample identifiers must be unique")
        missing = sorted(set(ordered_samples).difference(distances.index.astype(str)))
        if missing:
            preview = ", ".join(missing[:5])
            raise ValueError(f"distance matrix is missing {len(missing)} required sample(s): {preview}")
    matrix = distances.loc[ordered_samples, ordered_samples].to_numpy(dtype=float)
    n_samples = matrix.shape[0]
    if n_samples < 2:
        raise ValueError("at least two samples are required for MDS")
    centering = np.eye(n_samples) - np.full((n_samples, n_samples), 1.0 / n_samples)
    gram = -0.5 * centering @ (matrix**2) @ centering
    eigenvalues, eigenvectors = np.linalg.eigh(gram)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]
    tolerance = max(float(np.max(np.abs(eigenvalues))), 1.0) * 1e-12
    positive = eigenvalues > tolerance
    if not np.any(positive):
        raise ValueError("distance matrix has no positive MDS dimensions")
    positive_values = eigenvalues[positive]
    positive_vectors = eigenvectors[:, positive]
    n_components = min(int(n_components), positive_values.size)
    if n_components < 1:
        raise ValueError("n_components must be at least 1")
    coordinates = positive_vectors[:, :n_components] * np.sqrt(positive_values[:n_components])
    for index in range(n_components):
        anchor = int(np.argmax(np.abs(coordinates[:, index])))
        if coordinates[anchor, index] < 0:
            coordinates[:, index] *= -1
    total_positive = float(positive_values.sum())
    variance = positive_values[:n_components] / total_positive if total_positive else np.zeros(n_components)
    result = pd.DataFrame({"sample_id": ordered_samples})
    for index in range(n_components):
        result[f"MDS{index + 1}"] = coordinates[:, index]
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
