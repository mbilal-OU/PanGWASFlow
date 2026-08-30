from pathlib import Path

from pangwasflow.benchmark import make_demo, run_demo


def test_demo_is_reproducible():
    x1, y1, truth1 = make_demo(seed=42)
    x2, y2, truth2 = make_demo(seed=42)
    assert x1.equals(x2)
    assert y1.equals(y2)
    assert truth1.equals(truth2)


def test_demo_writes_expected_outputs(tmp_path: Path):
    run_demo(tmp_path, seed=42, samples=120, features=40, pcs=2)
    expected = [
        "summary.json",
        "features.tsv",
        "metadata.tsv",
        "feature_truth.tsv",
        "pca_scores.tsv",
        "association_naive.tsv",
        "association_adjusted.tsv",
        "figures/pca_structure.svg",
        "figures/association_comparison.svg",
        "figures/qq_comparison.svg",
    ]
    for relative in expected:
        assert (tmp_path / relative).exists()
