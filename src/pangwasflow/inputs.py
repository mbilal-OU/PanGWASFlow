from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import json

import numpy as np
import pandas as pd


ROARY_METADATA_COLUMNS = {
    "Gene",
    "Non-unique Gene name",
    "Annotation",
    "No. isolates",
    "No. sequences",
    "Avg sequences per isolate",
    "Genome Fragment",
    "Order within Fragment",
    "Accessory Fragment",
    "Accessory Order with Fragment",
    "QC",
    "Min group size nuc",
    "Max group size nuc",
    "Avg group size nuc",
}


@dataclass(frozen=True)
class InputSummary:
    samples_metadata: int
    samples_features: int
    samples_retained: int
    features_input: int
    features_retained: int
    cases: int
    controls: int
    min_prevalence: float
    max_prevalence: float

    def to_dict(self) -> dict[str, int | float]:
        return self.__dict__.copy()


def _read_table(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path)
    return pd.read_csv(path, sep="\t")


def _validate_sample_ids(table: pd.DataFrame, sample_id_column: str, table_name: str) -> pd.Series:
    if sample_id_column not in table.columns:
        raise ValueError(f"{table_name} must contain sample ID column '{sample_id_column}'")
    values = table[sample_id_column].astype(str).str.strip()
    if (values == "").any():
        raise ValueError(f"{table_name} contains empty sample IDs")
    duplicated = values[values.duplicated()].unique().tolist()
    if duplicated:
        preview = ", ".join(duplicated[:5])
        raise ValueError(f"{table_name} contains duplicate sample IDs: {preview}")
    return values


def _unique_feature_names(raw_names: pd.Series, prefix: str = "feature") -> list[str]:
    names: list[str] = []
    seen: dict[str, int] = {}
    for index, raw in enumerate(raw_names.astype(str)):
        base = raw.strip() or f"{prefix}_{index:05d}"
        count = seen.get(base, 0)
        seen[base] = count + 1
        names.append(base if count == 0 else f"{base}__{count + 1}")
    return names


def read_metadata(
    path: str | Path,
    sample_id_column: str = "sample_id",
    label_column: str = "label",
) -> pd.DataFrame:
    """Read strict binary phenotype metadata and normalize key column names."""
    table = _read_table(path).copy()
    table[sample_id_column] = _validate_sample_ids(table, sample_id_column, "metadata")
    if label_column not in table.columns:
        raise ValueError(f"metadata must contain phenotype column '{label_column}'")
    label = pd.to_numeric(table[label_column], errors="coerce")
    if label.isna().any() or not set(label.unique()).issubset({0, 1}):
        raise ValueError("phenotype must be encoded strictly as 0/1; recode labels before analysis")
    table = table.rename(columns={sample_id_column: "sample_id", label_column: "label"})
    table["sample_id"] = table["sample_id"].astype(str)
    table["label"] = label.astype(int)
    return table


def read_binary_feature_matrix(path: str | Path, sample_id_column: str = "sample_id") -> pd.DataFrame:
    """Read a sample-by-feature matrix encoded strictly as complete 0/1 data."""
    table = _read_table(path).copy()
    table[sample_id_column] = _validate_sample_ids(table, sample_id_column, "feature matrix")
    feature_columns = [column for column in table.columns if column != sample_id_column]
    if not feature_columns:
        raise ValueError("feature matrix contains no feature columns")
    numeric = table[feature_columns].apply(pd.to_numeric, errors="coerce")
    if numeric.isna().any().any():
        raise ValueError("feature matrix contains missing or non-numeric values; current matrix adapter requires complete 0/1 data")
    invalid = ~numeric.isin([0, 1])
    if invalid.any().any():
        bad_column = invalid.any(axis=0).idxmax()
        raise ValueError(f"feature matrix column '{bad_column}' contains values outside 0/1")
    result = numeric.astype(np.int8)
    result.insert(0, "sample_id", table[sample_id_column].astype(str))
    return result


def read_gene_presence_absence(path: str | Path) -> pd.DataFrame:
    """Convert a Roary/Panaroo-style gene_presence_absence.csv to a binary sample-by-gene matrix."""
    table = pd.read_csv(path, dtype=str, keep_default_na=False)
    if "Gene" not in table.columns:
        raise ValueError("gene presence/absence table must contain a 'Gene' column")
    sample_columns = [column for column in table.columns if column not in ROARY_METADATA_COLUMNS]
    if not sample_columns:
        raise ValueError("no sample columns detected in gene presence/absence table")

    gene_names = _unique_feature_names(table["Gene"], prefix="gene_cluster")
    presence = table[sample_columns].apply(lambda column: column.astype(str).str.strip().ne(""))
    matrix = presence.to_numpy(dtype=np.int8).T
    result = pd.DataFrame(matrix, columns=gene_names)
    result.insert(0, "sample_id", [str(column).strip() for column in sample_columns])
    return result


def read_rtab_presence_absence(path: str | Path) -> pd.DataFrame:
    """Read SEER/pyseer-style Rtab data, preserving '.' and blank cells as missing."""
    table = pd.read_csv(
        path,
        sep="\t",
        compression="infer",
        dtype=str,
        keep_default_na=False,
    )
    if table.empty or table.shape[1] < 2:
        raise ValueError("Rtab input must contain a feature column and at least one sample column")

    feature_column = str(table.columns[0])
    feature_names = _unique_feature_names(table[feature_column], prefix="feature")
    sample_columns = [str(column).strip() for column in table.columns[1:]]
    if any(not sample for sample in sample_columns):
        raise ValueError("Rtab input contains an empty sample column name")
    if len(sample_columns) != len(set(sample_columns)):
        raise ValueError("Rtab input contains duplicate sample column names")

    raw = table.iloc[:, 1:].apply(lambda column: column.astype(str).str.strip())
    allowed = raw.isin(["0", "1", ".", ""])
    if not allowed.all().all():
        bad_columns = allowed.all(axis=0)
        bad_column = str(bad_columns[~bad_columns].index[0])
        bad_values = sorted(set(raw.loc[~allowed[bad_column], bad_column].astype(str)))
        preview = ", ".join(repr(value) for value in bad_values[:5])
        raise ValueError(
            f"Rtab sample column '{bad_column}' contains values outside 0/1/./blank: {preview}"
        )

    numeric = raw.replace({"0": 0.0, "1": 1.0, ".": np.nan, "": np.nan}).astype(float)
    result = pd.DataFrame(numeric.to_numpy(dtype=float).T, columns=feature_names)
    result.insert(0, "sample_id", sample_columns)
    return result


def align_samples(
    features: pd.DataFrame,
    metadata: pd.DataFrame,
    min_samples: int = 20,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, list[str]]]:
    """Intersect feature and metadata sample IDs while preserving metadata order."""
    feature_ids = set(features["sample_id"].astype(str))
    metadata_ids = metadata["sample_id"].astype(str).tolist()
    metadata_id_set = set(metadata_ids)
    retained = [sample for sample in metadata_ids if sample in feature_ids]
    if len(retained) < int(min_samples):
        raise ValueError(f"only {len(retained)} samples overlap; at least {min_samples} are required")
    dropped = {
        "metadata_only": [sample for sample in metadata_ids if sample not in feature_ids],
        "features_only": [sample for sample in features["sample_id"].astype(str) if sample not in metadata_id_set],
    }
    metadata_out = metadata.set_index("sample_id").loc[retained].reset_index()
    features_out = features.set_index("sample_id").loc[retained].reset_index()
    return features_out, metadata_out, dropped


def filter_features_by_prevalence(
    features: pd.DataFrame,
    min_prevalence: float = 0.01,
    max_prevalence: float = 0.99,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Filter incomplete and prevalence-extreme features, returning an auditable QC table.

    The current association/PCA layer requires complete binary features. Missing Rtab
    observations are therefore preserved during parsing and explicitly removed here at
    the feature level rather than being interpreted as absence or silently imputed.
    """
    if not 0 <= min_prevalence < max_prevalence <= 1:
        raise ValueError("prevalence thresholds must satisfy 0 <= min < max <= 1")
    columns = [column for column in features.columns if column != "sample_id"]
    values = features[columns].apply(pd.to_numeric, errors="coerce")
    prevalence = values.mean(axis=0, skipna=True)
    missingness = values.isna().mean(axis=0)
    complete = missingness.eq(0.0)
    within_prevalence = (prevalence >= min_prevalence) & (prevalence <= max_prevalence)
    keep = complete & within_prevalence
    reason = np.select(
        [~complete.to_numpy(), ~within_prevalence.to_numpy()],
        ["missingness_filter", "prevalence_filter"],
        default="retained",
    )
    qc = pd.DataFrame(
        {
            "feature": columns,
            "prevalence": prevalence.to_numpy(dtype=float),
            "missingness": missingness.to_numpy(dtype=float),
            "retained": keep.to_numpy(dtype=bool),
            "reason": reason,
        }
    )
    retained_columns = [column for column in columns if bool(keep[column])]
    if not retained_columns:
        raise ValueError("no complete variable features remain after missingness and prevalence filtering")
    retained_values = values[retained_columns].astype(np.int8)
    result = retained_values.copy()
    result.insert(0, "sample_id", features["sample_id"].astype(str).to_numpy())
    return result, qc


def prepare_inputs(
    feature_path: str | Path,
    metadata_path: str | Path,
    outdir: str | Path,
    feature_format: str = "matrix",
    sample_id_column: str = "sample_id",
    label_column: str = "label",
    min_prevalence: float = 0.01,
    max_prevalence: float = 0.99,
    min_samples: int = 20,
) -> InputSummary:
    """Normalize, align, QC, and write binary genomic features plus phenotype metadata."""
    if feature_format == "matrix":
        features = read_binary_feature_matrix(feature_path, sample_id_column=sample_id_column)
    elif feature_format in {"roary", "panaroo", "gene-pa"}:
        features = read_gene_presence_absence(feature_path)
    elif feature_format == "rtab":
        features = read_rtab_presence_absence(feature_path)
    else:
        raise ValueError("feature_format must be one of: matrix, roary, panaroo, gene-pa, rtab")
    metadata = read_metadata(metadata_path, sample_id_column=sample_id_column, label_column=label_column)
    features_input = features.shape[1] - 1
    features_aligned, metadata_aligned, dropped = align_samples(features, metadata, min_samples=min_samples)
    features_qc, qc = filter_features_by_prevalence(
        features_aligned,
        min_prevalence=min_prevalence,
        max_prevalence=max_prevalence,
    )
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    features_qc.to_csv(outdir / "features_qc.tsv", sep="\t", index=False)
    metadata_aligned.to_csv(outdir / "metadata_qc.tsv", sep="\t", index=False)
    qc.to_csv(outdir / "feature_qc.tsv", sep="\t", index=False)
    (outdir / "sample_alignment.json").write_text(json.dumps(dropped, indent=2) + "\n", encoding="utf-8")
    summary = InputSummary(
        samples_metadata=int(metadata.shape[0]),
        samples_features=int(features.shape[0]),
        samples_retained=int(metadata_aligned.shape[0]),
        features_input=int(features_input),
        features_retained=int(features_qc.shape[1] - 1),
        cases=int(metadata_aligned["label"].sum()),
        controls=int((metadata_aligned["label"] == 0).sum()),
        min_prevalence=float(min_prevalence),
        max_prevalence=float(max_prevalence),
    )
    (outdir / "input_summary.json").write_text(json.dumps(summary.to_dict(), indent=2) + "\n", encoding="utf-8")
    return summary
