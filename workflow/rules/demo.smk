rule synthetic_demo:
    output:
        summary="results/demo/summary.json",
        features="results/demo/features.tsv",
        metadata="results/demo/metadata.tsv",
        truth="results/demo/feature_truth.tsv",
        pca="results/demo/pca_scores.tsv",
        naive="results/demo/association_naive.tsv",
        adjusted="results/demo/association_adjusted.tsv",
        pca_figure="results/demo/figures/pca_structure.svg",
        association_figure="results/demo/figures/association_comparison.svg",
        qq_figure="results/demo/figures/qq_comparison.svg"
    params:
        seed=lambda wc: config["demo"]["seed"],
        samples=lambda wc: config["demo"]["samples"],
        features=lambda wc: config["demo"]["features"],
        pcs=lambda wc: config["demo"]["pcs"]
    shell:
        "pangwasflow demo --outdir results/demo --seed {params.seed} "
        "--samples {params.samples} --features {params.features} --pcs {params.pcs}"
