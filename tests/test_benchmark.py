from pathlib import Path

from pangwasflow.benchmark import make_demo, run_demo


def test_demo_is_reproducible():
    x1, y1 = make_demo(seed=42)
    x2, y2 = make_demo(seed=42)
    assert x1.equals(x2)
    assert y1.equals(y2)


def test_demo_writes_expected_outputs(tmp_path: Path):
    run_demo(tmp_path, seed=42, samples=120, features=40)
    assert (tmp_path / "summary.json").exists()
    assert (tmp_path / "association_naive.tsv").exists()
