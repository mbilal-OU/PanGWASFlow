from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def render(artifact: Path, assets: Path, preview: Path) -> None:
    assets.mkdir(parents=True, exist_ok=True)
    preview.mkdir(parents=True, exist_ok=True)

    metadata = pd.read_csv(artifact / "metadata_qc.tsv", sep="\t")
    structure = pd.read_csv(artifact / "structure_covariates.tsv", sep="\t")
    naive = pd.read_csv(artifact / "association_naive.tsv", sep="\t")
    adjusted = pd.read_csv(artifact / "association_adjusted.tsv", sep="\t")
    input_summary = json.loads((artifact / "input_summary.json").read_text())
    gwas_summary = json.loads((artifact / "gwas_summary.json").read_text())

    (preview / "input_summary.json").write_text(json.dumps(input_summary, indent=2) + "\n")
    (preview / "gwas_summary.json").write_text(json.dumps(gwas_summary, indent=2) + "\n")
    adjusted.head(30).to_csv(preview / "top30_adjusted.tsv", sep="\t", index=False)

    plot = structure.merge(metadata, on="sample_id", validate="one_to_one")
    fig, ax = plt.subplots(figsize=(6.2, 4.8))
    for label, name, marker in [(0, "Susceptible", "o"), (1, "Resistant", "^")]:
        subset = plot[plot["label"] == label]
        ax.scatter(
            subset["MDS1"],
            subset["MDS2"],
            s=22,
            alpha=0.72,
            label=f"{name} (n={len(subset)})",
            marker=marker,
            linewidths=0.25,
            edgecolors="black",
        )
    ax.axhline(0, linewidth=0.5, color="0.8")
    ax.axvline(0, linewidth=0.5, color="0.8")
    variance = gwas_summary["structure_variance_explained"]
    ax.set_xlabel(f"MDS1 ({variance[0] * 100:.1f}% of positive-eigenvalue variation)")
    ax.set_ylabel(f"MDS2 ({variance[1] * 100:.1f}%)")
    ax.set_title("Mash-derived population structure")
    ax.legend(frameon=False, fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(assets / "pneumococcus_mds_structure.svg", format="svg")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.1))
    panels = [
        (naive, "Unadjusted Fisher scan", gwas_summary["lambda_naive"]),
        (adjusted, "Mash-MDS adjusted LRT", gwas_summary["lambda_adjusted"]),
    ]
    for ax, (table, title, inflation) in zip(axes, panels):
        pvalues = np.clip(table["pvalue"].to_numpy(float), np.finfo(float).tiny, 1.0)
        observed = -np.log10(np.sort(pvalues))
        expected = -np.log10((np.arange(1, len(pvalues) + 1) - 0.5) / len(pvalues))
        xmax = float(expected.max()) * 1.05
        ymax = float(observed.max()) * 1.03
        ax.plot([0, float(expected.max())], [0, float(expected.max())], linewidth=1, color="0.55")
        ax.scatter(expected, observed, s=10, alpha=0.68, linewidths=0, rasterized=True)
        ax.set_xlim(0, xmax)
        ax.set_ylim(0, ymax)
        ax.set_xlabel("Expected −log10(p)")
        ax.set_ylabel("Observed −log10(p)")
        ax.set_title(title)
        ax.text(0.04, 0.94, f"λ = {inflation:.2f}", transform=ax.transAxes, va="top", fontsize=9)
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(assets / "pneumococcus_qq_comparison.svg", format="svg")
    plt.close(fig)

    ranked = adjusted.sort_values("pvalue").reset_index(drop=True).copy()
    ranked["rank"] = np.arange(1, len(ranked) + 1)
    ranked["mlog10p"] = -np.log10(np.clip(ranked["pvalue"].to_numpy(float), np.finfo(float).tiny, 1.0))
    stable = ranked["effect_stable"].astype(bool)
    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    ax.scatter(
        ranked.loc[stable, "rank"],
        ranked.loc[stable, "mlog10p"],
        s=13,
        alpha=0.68,
        label="Stable effect estimate",
        linewidths=0,
        rasterized=True,
    )
    ax.scatter(
        ranked.loc[~stable, "rank"],
        ranked.loc[~stable, "mlog10p"],
        s=13,
        alpha=0.55,
        marker="x",
        label="Effect estimate flagged unstable",
        rasterized=True,
    )
    significant = ranked[ranked["qvalue"] < 0.05]
    if not significant.empty:
        boundary = float(significant["pvalue"].max())
        ax.axhline(-np.log10(boundary), linewidth=0.9, linestyle="--", color="0.45", label="BH FDR < 0.05 boundary")
    for feature in ["group_4276", "group_4417", "cpsG"]:
        row = ranked[ranked["feature"] == feature]
        if not row.empty:
            item = row.iloc[0]
            ax.annotate(feature, (item["rank"], item["mlog10p"]), xytext=(5, 5), textcoords="offset points", fontsize=8)
    ax.set_xlabel("Association rank")
    ax.set_ylabel("−log10(LRT p-value)")
    ax.set_title("Accessory-gene associations with penicillin resistance")
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(assets / "pneumococcus_accessory_rank.svg", format="svg")
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifact", type=Path, required=True)
    parser.add_argument("--assets", type=Path, default=Path("docs/assets"))
    parser.add_argument(
        "--preview",
        type=Path,
        default=Path("case_studies/pneumococcus_penicillin/results_preview"),
    )
    args = parser.parse_args()
    render(args.artifact, args.assets, args.preview)
