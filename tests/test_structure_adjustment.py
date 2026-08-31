import json

import numpy as np
import pandas as pd
import pytest

from pangwasflow.benchmark import make_demo, run_demo
from pangwasflow.structure import classical_mds_scores, pca_scores, read_distance_matrix


def test_pca_recovers_latent_structure():
    features, metadata, _ = make_demo()
    scores, _ = pca_scores(features, n_components=2)
    correlation = np.corrcoef(scores["PC1"], metadata["latent_group"])[0, 1]
    assert abs(correlation) > 0.9


def test_classical_mds_reconstructs_euclidean_distances():
    samples = ["s1", "s2", "s3", "s4"]
    coordinates = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 2.0], [1.0, 2.0]])
    distances = np.linalg.norm(coordinates[:, None, :] - coordinates[None, :, :], axis=2)
    matrix = pd.DataFrame(distances, index=samples, columns=samples)
    requested = ["s4", "s2", "s1"]
    scores, variance = classical_mds_scores(matrix, sample_ids=requested, n_components=2)
    recovered = scores[["MDS1", "MDS2"]].to_numpy()
    recovered_distances = np.linalg.norm(recovered[:, None, :] - recovered[None, :, :], axis=2)
    expected = matrix.loc[requested, requested].to_numpy()
    assert scores["sample_id"].tolist() == requested
    assert np.allclose(recovered_distances, expected, atol=1e-8)
    assert variance.size == 2


def test_distance_matrix_reader_rejects_asymmetry(tmp_path):
    path = tmp_path / "bad.tsv"
    pd.DataFrame([[0.0, 0.1], [0.2, 0.0]], index=["a", "b"], columns=["a", "b"]).to_csv(path, sep="\t")
    with pytest.raises(ValueError, match="symmetric"):
        read_distance_matrix(path)


def test_structure_adjustment_removes_lineage_false_positives(tmp_path):
    run_demo(tmp_path)
    summary = json.loads((tmp_path / "summary.json").read_text())
    assert summary["lambda_naive"] > 1.5
    assert 0.7 < summary["lambda_adjusted"] < 1.3
    assert summary["significant_lineage_markers_naive"] >= 30
    assert summary["significant_lineage_markers_adjusted"] == 0
    assert summary["causal_qvalue_adjusted"] < 0.05
    assert (tmp_path / "figures" / "pca_structure.svg").exists()
    assert (tmp_path / "figures" / "association_comparison.svg").exists()
    assert (tmp_path / "figures" / "qq_comparison.svg").exists()
