from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from pangwasflow.inputs import prepare_inputs, read_rtab_presence_absence


def test_rtab_transposes_gene_by_sample_matrix(tmp_path: Path):
    path = tmp_path / "presence_absence.Rtab"
    pd.DataFrame(
        {
            "Gene": ["geneA", "geneB", "geneC"],
            "sample_1": [1, 0, 1],
            "sample_2": [0, 1, 1],
            "sample_3": [1, 1, 0],
        }
    ).to_csv(path, sep="\t", index=False)
    matrix = read_rtab_presence_absence(path)
    assert matrix.columns.tolist() == ["sample_id", "geneA", "geneB", "geneC"]
    assert matrix["sample_id"].tolist() == ["sample_1", "sample_2", "sample_3"]
    assert matrix[["geneA", "geneB", "geneC"]].values.tolist() == [[1, 0, 1], [0, 1, 1], [1, 1, 0]]


def test_rtab_preserves_dot_and_blank_as_missing(tmp_path: Path):
    path = tmp_path / "missing.Rtab"
    path.write_text(
        "Gene\tsample_1\tsample_2\tsample_3\n"
        "geneA\t1\t.\t0\n"
        "geneB\t0\t\t1\n",
        encoding="utf-8",
    )
    matrix = read_rtab_presence_absence(path)
    assert matrix.loc[matrix["sample_id"] == "sample_2", "geneA"].isna().item()
    assert matrix.loc[matrix["sample_id"] == "sample_2", "geneB"].isna().item()
    assert matrix.loc[matrix["sample_id"] == "sample_1", "geneA"].item() == 1


def test_rtab_rejects_non_binary_values(tmp_path: Path):
    path = tmp_path / "bad.Rtab"
    pd.DataFrame({"Gene": ["geneA"], "sample_1": [2], "sample_2": [0]}).to_csv(path, sep="\t", index=False)
    with pytest.raises(ValueError, match="outside 0/1"):
        read_rtab_presence_absence(path)


def test_prepare_drops_missing_rtab_feature_without_treating_missing_as_absence(tmp_path: Path):
    samples = [f"sample_{index}" for index in range(1, 25)]
    rtab = tmp_path / "genes_missing.Rtab"
    phenotype = tmp_path / "phenotype.tsv"
    rows = [
        ["complete_variable", *[str(index % 2) for index in range(24)]],
        ["incomplete_variable", *(["."] + [str(index % 2) for index in range(1, 24)])],
    ]
    rtab.write_text(
        "Gene\t" + "\t".join(samples) + "\n" + "\n".join("\t".join(row) for row in rows) + "\n",
        encoding="utf-8",
    )
    pd.DataFrame({"samples": samples, "binary": [index % 2 for index in range(24)]}).to_csv(
        phenotype, sep="\t", index=False
    )
    outdir = tmp_path / "prepared_missing"
    summary = prepare_inputs(
        rtab,
        phenotype,
        outdir,
        feature_format="rtab",
        sample_id_column="samples",
        label_column="binary",
        min_prevalence=0.05,
        max_prevalence=0.95,
        min_samples=20,
    )
    assert summary.features_retained == 1
    retained = pd.read_csv(outdir / "features_qc.tsv", sep="\t")
    assert retained.columns.tolist() == ["sample_id", "complete_variable"]
    qc = pd.read_csv(outdir / "feature_qc.tsv", sep="\t")
    incomplete = qc.loc[qc["feature"] == "incomplete_variable"].iloc[0]
    assert np.isclose(incomplete["missingness"], 1 / 24)
    assert incomplete["reason"] == "missingness_filter"


def test_rtab_can_enter_standard_preparation_path(tmp_path: Path):
    samples = [f"sample_{index}" for index in range(1, 25)]
    rtab = tmp_path / "genes.Rtab"
    phenotype = tmp_path / "phenotype.tsv"
    pd.DataFrame(
        {
            "Gene": ["core_gene", "variable_gene", "rare_gene"],
            **{
                sample: [1, int(index % 2 == 0), int(index == 1)]
                for index, sample in enumerate(samples, start=1)
            },
        }
    ).to_csv(rtab, sep="\t", index=False)
    pd.DataFrame({"samples": samples, "binary": [index % 2 for index in range(24)]}).to_csv(
        phenotype, sep="\t", index=False
    )
    outdir = tmp_path / "prepared"
    summary = prepare_inputs(
        rtab,
        phenotype,
        outdir,
        feature_format="rtab",
        sample_id_column="samples",
        label_column="binary",
        min_prevalence=0.05,
        max_prevalence=0.95,
        min_samples=20,
    )
    assert summary.samples_retained == 24
    assert summary.features_input == 3
    assert summary.features_retained == 1
    retained = pd.read_csv(outdir / "features_qc.tsv", sep="\t")
    assert retained.columns.tolist() == ["sample_id", "variable_gene"]
