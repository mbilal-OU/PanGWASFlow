rule synthetic_demo:
    output:
        summary="results/demo/summary.json",
        features="results/demo/features.tsv",
        metadata="results/demo/metadata.tsv",
        association="results/demo/association_naive.tsv"
    params:
        seed=lambda wc: config["demo"]["seed"],
        samples=lambda wc: config["demo"]["samples"],
        features=lambda wc: config["demo"]["features"]
    shell:
        "pangwasflow demo --outdir results/demo --seed {params.seed} "
        "--samples {params.samples} --features {params.features}"
