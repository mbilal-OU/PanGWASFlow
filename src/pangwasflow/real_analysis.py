from __future__ import annotations

from pathlib import Path
import json

import pandas as pd

from .association import adjusted_logistic_scan, fisher_scan
from .structure import classical_mds_scores, genomic_inflation, pca_scores, read_distance_matrix


def run_binary_gwas(
    features: pd.DataFrame,
    metadata: pd.DataFrame,
    outdir: str | Path,
    pcs: int = 2,
    covariates: pd.DataFrame | None = None,
    structure_method: str = "feature_pca",
    structure_variance: list[float] | None = None,
) -> dict[str, object]:
    """Run baseline and structure-adjusted GWAS models on a prepared binary feature matrix."""
    feature_samples = set(features["sample_id"].astype(str))
    metadata_samples = set(metadata["sample_id"].astype(str))
    if feature_samples != metadata_samples:
        raise ValueError("features and metadata must contain the same prepared sample set")
    if metadata["label"].nunique() != 2:
        raise ValueError("both phenotype classes must be present")

    if covariates is None:
        scores, variance = pca_scores(features, n_components=pcs)
        structure_method = "feature_pca"
    else:
        if "sample_id" not in covariates.columns:
            raise ValueError("structure covariates must contain sample_id")
        if covariates["sample_id"].astype(str).duplicated().any():
            raise ValueError("structure covariates contain duplicate sample identifiers")
        covariate_samples = set(covariates["sample_id"].astype(str))
        missing = metadata_samples.difference(covariate_samples)
        if missing:
            preview = ", ".join(sorted(missing)[:5])
            raise ValueError(f"structure covariates are missing {len(missing)} prepared sample(s): {preview}")
        scores = covariates[covariates["sample_id"].astype(str).isin(metadata_samples)].copy()
        score_columns = [column for column in scores.columns if column != "sample_id"]
        if not score_columns:
            raise ValueError("structure covariates must contain at least one numeric covariate")
        scores[score_columns] = scores[score_columns].apply(pd.to_numeric, errors="raise")
        variance = pd.Series(structure_variance if structure_variance is not None else [], dtype=float).to_numpy()

    naive = fisher_scan(features, metadata)
    adjusted = adjusted_logistic_scan(features, metadata, scores)
    score_columns = [column for column in scores.columns if column != "sample_id"]
    summary: dict[str, object] = {
        "samples": int(metadata.shape[0]),
        "features": int(features.shape[1] - 1),
        "cases": int(metadata["label"].sum()),
        "controls": int((metadata["label"] == 0).sum()),
        "structure_method": structure_method,
        "structure_components": int(len(score_columns)),
        "structure_variance_explained": [float(value) for value in variance],
        "pcs": int(len(score_columns)),
        "pc_variance_explained": [float(value) for value in variance],
        "lambda_naive": float(genomic_inflation(naive["pvalue"].to_numpy())),
        "lambda_adjusted": float(genomic_inflation(adjusted["pvalue"].to_numpy())),
        "significant_features_naive": int((naive["qvalue"] < 0.05).sum()),
        "significant_features_adjusted": int((adjusted["qvalue"] < 0.05).sum()),
    }
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    scores.to_csv(outdir / "structure_covariates.tsv", sep="\t", index=False)
    if structure_method == "feature_pca":
        scores.to_csv(outdir / "pca_scores.tsv", sep="\t", index=False)
    naive.to_csv(outdir / "association_naive.tsv", sep="\t", index=False)
    adjusted.to_csv(outdir / "association_adjusted.tsv", sep="\t", index=False)
    (outdir / "gwas_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary


def run_prepared_paths(
    features_path: str | Path,
    metadata_path: str | Path,
    outdir: str | Path,
    pcs: int = 2,
    distance_matrix_path: str | Path | None = None,
) -> dict[str, object]:
    features = pd.read_csv(features_path, sep="\t")
    metadata = pd.read_csv(metadata_path, sep="\t")
    if distance_matrix_path is None:
        return run_binary_gwas(features, metadata, outdir=outdir, pcs=pcs)

    distances = read_distance_matrix(distance_matrix_path)
    scores, variance = classical_mds_scores(
        distances,
        sample_ids=metadata["sample_id"].astype(str).tolist(),
        n_components=pcs,
    )
    return run_binary_gwas(
        features,
        metadata,
        outdir=outdir,
        pcs=pcs,
        covariates=scores,
        structure_method="distance_mds",
        structure_variance=[float(value) for value in variance],
    )
