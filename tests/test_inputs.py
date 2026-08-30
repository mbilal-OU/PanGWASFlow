from pathlib import Path

import pandas as pd
import pytest

from pangwasflow.inputs import (
    align_samples,
    filter_features_by_prevalence,
    prepare_inputs,
    read_binary_feature_matrix,
    read_gene_presence_absence,
)


def test_binary_matrix_rejects_non_binary_values(tmp_path: Path):
    path = tmp_path / "features.tsv"
    pd.DataFrame({"sample_id": ["a", "b"], "snp_1": [0, 2]}).to_csv(path, sep="\t", index=False)
    with pytest.raises(ValueError, match="outside 0/1"):
        read_binary_feature_matrix(path)


def test_gene_presence_absence_parser(tmp_path: Path):
    path = tmp_path / "gene_presence_absence.csv"
    pd.DataFrame(
        {
            "Gene": ["geneA", "geneB"],
            "No. isolates": [1, 2],
            "No. sequences": [1, 2],
            "Avg group size nuc": [800, 900],
            "sample_A": ["locus_1", ""],
            "sample_B": ["", "locus_2"],
        }
    ).to_csv(path, index=False)
    matrix = read_gene_presence_absence(path)
    assert matrix["sample_id"].tolist() == ["sample_A", "sample_B"]
    assert matrix[["geneA", "geneB"]].values.tolist() == [[1, 0], [0, 1]]


def test_alignment_and_prevalence_qc():
    features = pd.DataFrame({"sample_id": [f"s{i}" for i in range(6)], "common": [1, 1, 1, 1, 1, 1], "variable": [0, 0, 0, 1, 1, 1]})
    metadata = pd.DataFrame({"sample_id": [f"s{i}" for i in range(1, 7)], "label": [0, 0, 1, 1, 0, 1]})
    x, y, dropped = align_samples(features, metadata, min_samples=5)
    assert x.shape[0] == 5
    assert dropped["metadata_only"] == ["s6"]
    filtered, qc = filter_features_by_prevalence(x, min_prevalence=0.2, max_prevalence=0.8)
    assert filtered.columns.tolist() == ["sample_id", "variable"]
    assert qc.loc[qc["feature"] == "common", "retained"].item() == False


def test_prepare_inputs_writes_auditable_outputs(tmp_path: Path):
    samples = [f"s{i:02d}" for i in range(24)]
    feature_path = tmp_path / "features.tsv"
    metadata_path = tmp_path / "metadata.tsv"
    pd.DataFrame({"sample_id": samples, "snp_1": [0] * 12 + [1] * 12, "snp_2": [0, 1] * 12}).to_csv(feature_path, sep="\t", index=False)
    pd.DataFrame({"sample_id": samples, "label": [0, 1] * 12}).to_csv(metadata_path, sep="\t", index=False)
    outdir = tmp_path / "prepared"
    summary = prepare_inputs(feature_path, metadata_path, outdir, min_samples=20)
    assert summary.samples_retained == 24
    assert (outdir / "features_qc.tsv").exists()
    assert (outdir / "feature_qc.tsv").exists()
    assert (outdir / "sample_alignment.json").exists()
