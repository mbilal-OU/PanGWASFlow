import json

import numpy as np

from pangwasflow.benchmark import make_demo, run_demo
from pangwasflow.structure import pca_scores


def test_pca_recovers_latent_structure():
    features, metadata, _ = make_demo()
    scores, _ = pca_scores(features, n_components=2)
    correlation = np.corrcoef(scores["PC1"], metadata["latent_group"])[0, 1]
    assert abs(correlation) > 0.9


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
