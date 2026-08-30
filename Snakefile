configfile: "config/config.yaml"

include: "workflow/rules/demo.smk"

rule all:
    input:
        "results/demo/summary.json",
        "results/demo/features.tsv",
        "results/demo/metadata.tsv",
        "results/demo/feature_truth.tsv",
        "results/demo/pca_scores.tsv",
        "results/demo/association_naive.tsv",
        "results/demo/association_adjusted.tsv",
        "results/demo/figures/pca_structure.svg",
        "results/demo/figures/association_comparison.svg",
        "results/demo/figures/qq_comparison.svg"
