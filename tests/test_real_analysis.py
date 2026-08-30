from pathlib import Path

from pangwasflow.benchmark import make_demo
from pangwasflow.real_analysis import run_binary_gwas


def test_prepared_binary_gwas_writes_expected_outputs(tmp_path: Path):
    features, metadata, _ = make_demo(seed=42, samples=160, features=50)
    summary = run_binary_gwas(features, metadata[["sample_id", "label"]], tmp_path, pcs=2)
    assert summary["samples"] == 160
    assert summary["features"] == 50
    assert (tmp_path / "pca_scores.tsv").exists()
    assert (tmp_path / "association_naive.tsv").exists()
    assert (tmp_path / "association_adjusted.tsv").exists()
    assert (tmp_path / "gwas_summary.json").exists()
