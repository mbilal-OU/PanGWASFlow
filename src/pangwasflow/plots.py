from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def _prepare() -> None:
    plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9, "legend.fontsize": 8, "figure.dpi": 150, "savefig.bbox": "tight", "svg.fonttype": "none"})


def plot_structure(metadata: pd.DataFrame, scores: pd.DataFrame, variance: np.ndarray, output: Path) -> None:
    _prepare()
    merged = metadata.merge(scores, on="sample_id", how="inner")
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    scatter = ax.scatter(merged["PC1"], merged["PC2"], c=merged["latent_group"], s=22, alpha=0.75)
    ax.set_xlabel(f"PC1 ({variance[0] * 100:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({variance[1] * 100:.1f}% variance)")
    ax.set_title("Synthetic population structure")
    handles, _ = scatter.legend_elements()
    ax.legend(handles, ["Lineage 0", "Lineage 1"], title="Benchmark lineage", frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(output, format="svg")
    plt.close(fig)


def plot_association_comparison(naive: pd.DataFrame, adjusted: pd.DataFrame, truth: pd.DataFrame, output: Path) -> None:
    _prepare()
    base = truth.copy()
    base["feature_index"] = np.arange(len(base))
    left = base.merge(naive[["feature", "pvalue"]], on="feature", how="left")
    right = base.merge(adjusted[["feature", "pvalue"]], on="feature", how="left")
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 6.2), sharex=True)
    for ax, table, title in zip(axes, [left, right], ["Naive association", "PCA-adjusted association"]):
        y = -np.log10(np.clip(table["pvalue"].to_numpy(dtype=float), 1e-300, 1.0))
        ax.scatter(table["feature_index"], y, s=18, alpha=0.8)
        causal = table["role"] == "causal"
        lineage = table["role"] == "lineage_marker"
        ax.scatter(table.loc[lineage, "feature_index"], y[lineage], s=21, marker="s", label="Lineage marker")
        ax.scatter(table.loc[causal, "feature_index"], y[causal], s=38, marker="*", label="Causal feature")
        ax.axhline(-np.log10(0.05 / len(table)), linestyle="--", linewidth=0.9, label="Bonferroni 0.05")
        ax.set_ylabel("-log10(p)")
        ax.set_title(title, loc="left")
        ax.spines[["top", "right"]].set_visible(False)
    axes[-1].set_xlabel("Synthetic feature index")
    axes[0].legend(frameon=False, ncol=3)
    fig.suptitle("Population structure changes the association landscape", y=1.01, fontsize=11)
    fig.tight_layout()
    fig.savefig(output, format="svg")
    plt.close(fig)


def plot_qq(naive: pd.DataFrame, adjusted: pd.DataFrame, lambda_naive: float, lambda_adjusted: float, output: Path) -> None:
    _prepare()
    fig, ax = plt.subplots(figsize=(5.2, 5.0))
    for table, label in [(naive, f"Naive (lambda={lambda_naive:.2f})"), (adjusted, f"Adjusted (lambda={lambda_adjusted:.2f})")]:
        p = np.sort(np.clip(table["pvalue"].to_numpy(dtype=float), 1e-300, 1.0))
        observed = -np.log10(p)
        expected = -np.log10((np.arange(1, p.size + 1) - 0.5) / p.size)
        ax.plot(expected, observed, marker="o", markersize=3, linewidth=0.9, label=label)
    limit = max(ax.get_xlim()[1], ax.get_ylim()[1])
    ax.plot([0, limit], [0, limit], linestyle="--", linewidth=0.8)
    ax.set_xlim(0, limit)
    ax.set_ylim(0, limit)
    ax.set_xlabel("Expected -log10(p)")
    ax.set_ylabel("Observed -log10(p)")
    ax.set_title("QQ diagnostic")
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(output, format="svg")
    plt.close(fig)


def write_demo_figures(metadata: pd.DataFrame, scores: pd.DataFrame, variance: np.ndarray, naive: pd.DataFrame, adjusted: pd.DataFrame, truth: pd.DataFrame, lambda_naive: float, lambda_adjusted: float, outdir: str | Path) -> None:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    plot_structure(metadata, scores, variance, outdir / "pca_structure.svg")
    plot_association_comparison(naive, adjusted, truth, outdir / "association_comparison.svg")
    plot_qq(naive, adjusted, lambda_naive, lambda_adjusted, outdir / "qq_comparison.svg")
