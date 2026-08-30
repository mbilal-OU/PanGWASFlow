from __future__ import annotations

from pathlib import Path
import json

import pandas as pd

from .association import adjusted_logistic_scan, fisher_scan
from .structure import genomic_inflation, pca_scores


def run_binary_gwas(
    features: pd.DataFrame,
    metadata: pd.DataFrame,
    outdir: str | Path,
    pcs: int = 2,
) -> dict[str, object]:
    """Run the current baseline and PCA-adjusted GWAS models on a prepared binary feature matrix."""
    if set(features["sample_id"].astype(str)) != set(metadata["sample_id"].astype(str)):
        raise ValueError("features and metadata must contain the same prepared sample set")
    if metadata["label"].nunique() != 2:
        raise ValueError("both phenotype classes must be present")
    scores, variance = pca_scores(features, n_components=pcs)
    naive = fisher_scan(features, metadata)
    adjusted = adjusted_logistic_scan(features, metadata, scores)
    summary: dict[str, object] = {
        "samples": int(metadata.shape[0]),
        "features": int(features.shape[1] - 1),
        "cases": int(metadata["label"].sum()),
        "controls": int((metadata["label"] == 0).sum()),
        "pcs": int(scores.shape[1] - 1),
        "pc_variance_explained": [float(value) for value in variance],
        "lambda_naive": float(genomic_inflation(naive["pvalue"].to_numpy())),
        "lambda_adjusted": float(genomic_inflation(adjusted["pvalue"].to_numpy())),
        "significant_features_naive": int((naive["qvalue"] < 0.05).sum()),
        "significant_features_adjusted": int((adjusted["qvalue"] < 0.05).sum()),
    }
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
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
) -> dict[str, object]:
    features = pd.read_csv(features_path, sep="\t")
    metadata = pd.read_csv(metadata_path, sep="\t")
    return run_binary_gwas(features, metadata, outdir=outdir, pcs=pcs)
