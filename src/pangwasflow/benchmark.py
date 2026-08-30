from __future__ import annotations

from pathlib import Path
import json

import numpy as np
import pandas as pd
from scipy.special import expit
from scipy.stats import fisher_exact

from .association import adjusted_logistic_scan
from .multiple_testing import benjamini_hochberg
from .plots import write_demo_figures
from .structure import genomic_inflation, pca_scores


def make_demo(seed: int = 42, samples: int = 400, features: int = 120) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Create a deterministic population-structure stress test."""
    if features < 40:
        raise ValueError("the structure benchmark requires at least 40 features")
    rng = np.random.default_rng(seed)
    group = rng.binomial(1, 0.5, samples)
    matrix = rng.binomial(1, 0.20, size=(samples, features))
    matrix[:, 0] = rng.binomial(1, 0.30, samples)
    matrix[:, 1] = rng.binomial(1, np.where(group == 1, 0.65, 0.15))
    for index in range(2, 36):
        matrix[:, index] = rng.binomial(1, np.where(group == 1, 0.75, 0.08))
    probability = expit(-2.4 + 1.7 * matrix[:, 0] + 1.8 * group)
    label = rng.binomial(1, probability)
    feature_names = [f"feature_{index:03d}" for index in range(features)]
    feature_table = pd.DataFrame(matrix, columns=feature_names)
    feature_table.insert(0, "sample_id", [f"sample_{index:03d}" for index in range(samples)])
    metadata = pd.DataFrame({"sample_id": feature_table["sample_id"], "label": label, "latent_group": group})
    roles = ["causal"] + ["lineage_marker"] * 35 + ["null"] * (features - 36)
    truth = pd.DataFrame({"feature": feature_names, "role": roles})
    return feature_table, metadata, truth


def fisher_scan(features: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    """Run an unadjusted 2x2 association scan."""
    merged = metadata.merge(features, on="sample_id", how="inner")
    rows = []
    for column in features.columns:
        if column == "sample_id":
            continue
        table = pd.crosstab(merged[column], merged["label"]).reindex(index=[0, 1], columns=[0, 1], fill_value=0)
        odds_ratio, pvalue = fisher_exact(table.to_numpy())
        rows.append({"feature": column, "odds_ratio": float(odds_ratio), "pvalue": float(pvalue), "prevalence": float(merged[column].mean())})
    result = pd.DataFrame(rows)
    result["qvalue"] = benjamini_hochberg(result["pvalue"].to_numpy())
    return result.sort_values("pvalue").reset_index(drop=True)


def run_demo(outdir: str | Path, seed: int = 42, samples: int = 400, features: int = 120, pcs: int = 2) -> None:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    feature_table, metadata, truth = make_demo(seed=seed, samples=samples, features=features)
    scores, variance = pca_scores(feature_table, n_components=pcs)
    naive = fisher_scan(feature_table, metadata).merge(truth, on="feature", how="left")
    adjusted = adjusted_logistic_scan(feature_table, metadata, scores).merge(truth, on="feature", how="left")
    lambda_naive = genomic_inflation(naive["pvalue"].to_numpy())
    lambda_adjusted = genomic_inflation(adjusted["pvalue"].to_numpy())
    feature_table.to_csv(outdir / "features.tsv", sep="\t", index=False)
    metadata.to_csv(outdir / "metadata.tsv", sep="\t", index=False)
    truth.to_csv(outdir / "feature_truth.tsv", sep="\t", index=False)
    scores.to_csv(outdir / "pca_scores.tsv", sep="\t", index=False)
    naive.to_csv(outdir / "association_naive.tsv", sep="\t", index=False)
    adjusted.to_csv(outdir / "association_adjusted.tsv", sep="\t", index=False)
    naive_lineage = int(((naive["role"] == "lineage_marker") & (naive["qvalue"] < 0.05)).sum())
    adjusted_lineage = int(((adjusted["role"] == "lineage_marker") & (adjusted["qvalue"] < 0.05)).sum())
    causal_naive = naive.loc[naive["role"] == "causal"].iloc[0]
    causal_adjusted = adjusted.loc[adjusted["role"] == "causal"].iloc[0]
    summary = {"seed": seed, "samples": samples, "features": features, "pcs": pcs, "cases": int(metadata["label"].sum()), "controls": int((metadata["label"] == 0).sum()), "pc_variance_explained": [float(value) for value in variance], "lambda_naive": lambda_naive, "lambda_adjusted": lambda_adjusted, "significant_lineage_markers_naive": naive_lineage, "significant_lineage_markers_adjusted": adjusted_lineage, "causal_feature": "feature_000", "causal_qvalue_naive": float(causal_naive["qvalue"]), "causal_qvalue_adjusted": float(causal_adjusted["qvalue"])}
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    write_demo_figures(metadata, scores, variance, naive, adjusted, truth, lambda_naive, lambda_adjusted, outdir / "figures")
