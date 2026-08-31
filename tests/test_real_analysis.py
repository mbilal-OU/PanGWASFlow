from pathlib import Path

import numpy as np
import pandas as pd

from pangwasflow.benchmark import make_demo
from pangwasflow.real_analysis import run_binary_gwas, run_prepared_paths


def test_prepared_binary_gwas_writes_expected_outputs(tmp_path: Path):
    features, metadata, _ = make_demo(seed=42, samples=160, features=50)
    summary = run_binary_gwas(features, metadata[["sample_id", "label"]], tmp_path, pcs=2)
    assert summary["samples"] == 160
    assert summary["features"] == 50
    assert summary["structure_method"] == "feature_pca"
    assert (tmp_path / "pca_scores.tsv").exists()
    assert (tmp_path / "structure_covariates.tsv").exists()
    assert (tmp_path / "association_naive.tsv").exists()
    assert (tmp_path / "association_adjusted.tsv").exists()
    assert (tmp_path / "gwas_summary.json").exists()


def test_prepared_paths_can_use_distance_matrix(tmp_path: Path):
    features, metadata, _ = make_demo(seed=7, samples=80, features=40)
    feature_path = tmp_path / "features.tsv"
    metadata_path = tmp_path / "metadata.tsv"
    distance_path = tmp_path / "distances.tsv"
    features.to_csv(feature_path, sep="\t", index=False)
    metadata[["sample_id", "label"]].to_csv(metadata_path, sep="\t", index=False)

    latent = metadata["latent_group"].to_numpy(dtype=float)
    positions = np.column_stack([latent * 2.0, np.arange(latent.size, dtype=float) * 1e-3])
    distance_values = np.linalg.norm(positions[:, None, :] - positions[None, :, :], axis=2)
    sample_ids = metadata["sample_id"].astype(str).tolist()
    pd.DataFrame(distance_values, index=sample_ids, columns=sample_ids).to_csv(distance_path, sep="\t")

    outdir = tmp_path / "gwas"
    summary = run_prepared_paths(feature_path, metadata_path, outdir, pcs=2, distance_matrix_path=distance_path)
    assert summary["structure_method"] == "distance_mds"
    assert summary["structure_components"] == 2
    covariates = pd.read_csv(outdir / "structure_covariates.tsv", sep="\t")
    assert covariates.columns.tolist() == ["sample_id", "MDS1", "MDS2"]
    assert not (outdir / "pca_scores.tsv").exists()
