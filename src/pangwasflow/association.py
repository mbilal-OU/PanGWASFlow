from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import expit
from scipy.stats import chi2, fisher_exact, norm

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


def _fit_logistic(
    y: np.ndarray,
    design: np.ndarray,
    initial: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, float, bool, np.ndarray]:
    """Fit an unpenalized logistic model and return coefficients, SEs and NLL."""
    y = np.asarray(y, dtype=float)
    design = np.asarray(design, dtype=float)
    if design.ndim == 1:
        design = design[:, None]
    x = np.column_stack([np.ones(y.size), design])
    if initial is None:
        initial = np.zeros(x.shape[1], dtype=float)
    else:
        initial = np.asarray(initial, dtype=float)
        if initial.shape != (x.shape[1],):
            raise ValueError("initial logistic coefficients have the wrong shape")

    def objective(beta: np.ndarray) -> float:
        eta = x @ beta
        return float(np.sum(np.logaddexp(0.0, eta) - y * eta))

    def gradient(beta: np.ndarray) -> np.ndarray:
        eta = x @ beta
        return x.T @ (expit(eta) - y)

    fit = minimize(
        objective,
        initial,
        jac=gradient,
        method="BFGS",
        options={"gtol": 1e-7, "maxiter": 250},
    )
    beta = np.asarray(fit.x, dtype=float)
    fitted = expit(x @ beta)
    weights = fitted * (1.0 - fitted)
    information = x.T @ (x * weights[:, None])
    covariance = np.linalg.pinv(information)
    standard_error = np.sqrt(np.maximum(np.diag(covariance), 0.0))
    with np.errstate(divide="ignore", invalid="ignore"):
        zscore = beta / standard_error
    wald_pvalue = 2.0 * norm.sf(np.abs(zscore))
    wald_pvalue[~np.isfinite(wald_pvalue)] = 1.0
    nll = objective(beta)
    gradient_ok = np.linalg.norm(gradient(beta), ord=np.inf) < 1e-4
    success = bool(fit.success or gradient_ok)
    return beta, standard_error, nll, success, wald_pvalue


def adjusted_logistic_scan(
    features: pd.DataFrame,
    metadata: pd.DataFrame,
    covariates: pd.DataFrame,
) -> pd.DataFrame:
    """Scan binary features with a structure-adjusted logistic LRT.

    A structure-only logistic model is fitted once. Each genomic feature is then
    added to that model and assessed with a one-degree-of-freedom likelihood-ratio
    test. Wald quantities are retained as effect diagnostics, but LRT p-values are
    used for multiple-testing correction because they are more stable under strong
    or near-separated binary associations.
    """
    merged = metadata.merge(covariates, on="sample_id", how="inner").merge(features, on="sample_id", how="inner")
    covariate_columns = [column for column in covariates.columns if column != "sample_id"]
    y = merged["label"].to_numpy(dtype=float)
    covariate_matrix = merged[covariate_columns].to_numpy(dtype=float)

    null_beta, _, null_nll, null_success, _ = _fit_logistic(y, covariate_matrix)
    if not null_success:
        raise RuntimeError("structure-only logistic null model did not converge")

    rows: list[dict[str, float | str | bool]] = []
    min_pvalue = float(np.finfo(float).tiny)
    for feature in features.columns:
        if feature == "sample_id":
            continue
        values = merged[feature].to_numpy(dtype=float)
        prevalence = float(values.mean())
        if np.unique(values).size < 2:
            rows.append(
                {
                    "feature": feature,
                    "log_odds": float("nan"),
                    "odds_ratio": float("nan"),
                    "standard_error": float("nan"),
                    "wald_pvalue": 1.0,
                    "pvalue": 1.0,
                    "prevalence": prevalence,
                    "fit_success": True,
                    "effect_stable": False,
                }
            )
            continue

        design = np.column_stack([values, covariate_matrix])
        initial = np.concatenate(([null_beta[0], 0.0], null_beta[1:]))
        beta, standard_error, alt_nll, fit_success, wald_pvalue = _fit_logistic(y, design, initial=initial)
        likelihood_ratio = max(0.0, 2.0 * (null_nll - alt_nll))
        lrt_pvalue = max(float(chi2.sf(likelihood_ratio, df=1)), min_pvalue)
        feature_beta = float(beta[1])
        feature_se = float(standard_error[1])
        effect_stable = bool(
            fit_success
            and np.isfinite(feature_beta)
            and np.isfinite(feature_se)
            and feature_se > 1e-8
            and abs(feature_beta) < 20.0
        )
        rows.append(
            {
                "feature": feature,
                "log_odds": feature_beta,
                "odds_ratio": float(np.exp(np.clip(feature_beta, -700.0, 700.0))),
                "standard_error": feature_se,
                "wald_pvalue": float(wald_pvalue[1]),
                "pvalue": lrt_pvalue,
                "prevalence": prevalence,
                "fit_success": bool(fit_success),
                "effect_stable": effect_stable,
            }
        )

    result = pd.DataFrame(rows)
    result["qvalue"] = benjamini_hochberg(result["pvalue"].to_numpy())
    return result.sort_values("pvalue").reset_index(drop=True)
