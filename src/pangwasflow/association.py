from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import expit
from scipy.stats import fisher_exact, norm

from .multiple_testing import benjamini_hochberg


def fisher_scan(features: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    """Run an unadjusted 2x2 association scan for binary genomic features."""
    merged = metadata.merge(features, on="sample_id", how="inner")
    rows: list[dict[str, float | str]] = []
    for feature in features.columns:
        if feature == "sample_id":
            continue
        values = merged[feature].to_numpy(dtype=float)
        if np.unique(values).size < 2:
            rows.append({"feature": feature, "odds_ratio": float("nan"), "pvalue": 1.0, "prevalence": float(values.mean())})
            continue
        table = pd.crosstab(merged[feature], merged["label"]).reindex(index=[0, 1], columns=[0, 1], fill_value=0)
        odds_ratio, pvalue = fisher_exact(table.to_numpy())
        rows.append({"feature": feature, "odds_ratio": float(odds_ratio), "pvalue": float(pvalue), "prevalence": float(values.mean())})
    result = pd.DataFrame(rows)
    result["qvalue"] = benjamini_hochberg(result["pvalue"].to_numpy())
    return result.sort_values("pvalue").reset_index(drop=True)


def _logistic_wald(y: np.ndarray, design: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    y = np.asarray(y, dtype=float)
    design = np.asarray(design, dtype=float)
    x = np.column_stack([np.ones(y.size), design])

    def objective(beta: np.ndarray) -> float:
        eta = x @ beta
        return float(np.sum(np.logaddexp(0.0, eta) - y * eta))

    fit = minimize(objective, np.zeros(x.shape[1]), method="BFGS")
    beta = fit.x
    fitted = expit(x @ beta)
    weights = fitted * (1.0 - fitted)
    information = x.T @ (x * weights[:, None])
    covariance = np.linalg.pinv(information)
    standard_error = np.sqrt(np.maximum(np.diag(covariance), 0.0))
    with np.errstate(divide="ignore", invalid="ignore"):
        zscore = beta / standard_error
    pvalue = 2.0 * norm.sf(np.abs(zscore))
    pvalue[~np.isfinite(pvalue)] = 1.0
    return beta, standard_error, pvalue


def adjusted_logistic_scan(
    features: pd.DataFrame,
    metadata: pd.DataFrame,
    covariates: pd.DataFrame,
) -> pd.DataFrame:
    """Scan binary features with logistic regression adjusted for structure covariates."""
    merged = metadata.merge(covariates, on="sample_id", how="inner").merge(features, on="sample_id", how="inner")
    covariate_columns = [column for column in covariates.columns if column != "sample_id"]
    rows: list[dict[str, float | str]] = []
    for feature in features.columns:
        if feature == "sample_id":
            continue
        values = merged[feature].to_numpy(dtype=float)
        if np.unique(values).size < 2:
            rows.append({"feature": feature, "log_odds": float("nan"), "odds_ratio": float("nan"), "standard_error": float("nan"), "pvalue": 1.0, "prevalence": float(values.mean())})
            continue
        design = np.column_stack([values, merged[covariate_columns].to_numpy(dtype=float)])
        beta, standard_error, pvalue = _logistic_wald(merged["label"].to_numpy(dtype=float), design)
        rows.append({"feature": feature, "log_odds": float(beta[1]), "odds_ratio": float(np.exp(beta[1])), "standard_error": float(standard_error[1]), "pvalue": float(pvalue[1]), "prevalence": float(values.mean())})
    result = pd.DataFrame(rows)
    result["qvalue"] = benjamini_hochberg(result["pvalue"].to_numpy())
    return result.sort_values("pvalue").reset_index(drop=True)
