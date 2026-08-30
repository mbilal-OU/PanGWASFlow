from __future__ import annotations

from pathlib import Path
import json

import numpy as np
import pandas as pd
from scipy.special import expit
from scipy.stats import fisher_exact

from .multiple_testing import benjamini_hochberg


def make_demo(seed: int = 42, samples: int = 240, features: int = 120) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Create a deterministic binary-feature teaching benchmark."""
    rng = np.random.default_rng(seed)
    group = rng.binomial(1, 0.5, samples)
    matrix = rng.binomial(1, 0.20, size=(samples, features))

    matrix[:, 0] = rng.binomial(1, 0.30, samples)
    matrix[:, 1] = rng.binomial(1, np.where(group == 1, 0.82, 0.08))

    probability = expit(-2.2 + 1.8 * matrix[:, 0] + 1.4 * group)
    label = rng.binomial(1, probability)

    feature_names = [f"feature_{index:03d}" for index in range(features)]
    feature_table = pd.DataFrame(matrix, columns=feature_names)
    feature_table.insert(0, "sample_id", [f"sample_{index:03d}" for index in range(samples)])
    metadata = pd.DataFrame(
        {
            "sample_id": feature_table["sample_id"],
            "label": label,
            "group": group,
        }
    )
    return feature_table, metadata


def fisher_scan(features: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    """Run a simple unadjusted 2x2 association scan for the teaching benchmark."""
    merged = metadata.merge(features, on="sample_id", how="inner")
    rows = []
    for column in features.columns:
        if column == "sample_id":
            continue
        table = pd.crosstab(merged[column], merged["label"]).reindex(index=[0, 1], columns=[0, 1], fill_value=0)
        odds_ratio, pvalue = fisher_exact(table.to_numpy())
        rows.append({"feature": column, "odds_ratio": float(odds_ratio), "pvalue": float(pvalue)})
    result = pd.DataFrame(rows)
    result["qvalue"] = benjamini_hochberg(result["pvalue"].to_numpy())
    return result.sort_values("pvalue").reset_index(drop=True)


def run_demo(outdir: str | Path, seed: int = 42, samples: int = 240, features: int = 120) -> None:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    feature_table, metadata = make_demo(seed=seed, samples=samples, features=features)
    results = fisher_scan(feature_table, metadata)

    feature_table.to_csv(outdir / "features.tsv", sep="\t", index=False)
    metadata.to_csv(outdir / "metadata.tsv", sep="\t", index=False)
    results.to_csv(outdir / "association_naive.tsv", sep="\t", index=False)

    summary = {
        "seed": seed,
        "samples": samples,
        "features": features,
        "cases": int(metadata["label"].sum()),
        "controls": int((metadata["label"] == 0).sum()),
        "top_feature": str(results.iloc[0]["feature"]),
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
