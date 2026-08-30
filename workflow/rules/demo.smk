rule synthetic_demo:
    output:
        summary="results/demo/summary.json",
        naive="results/demo/association_naive.tsv",
        adjusted="results/demo/association_adjusted.tsv",
        manhattan_naive="results/demo/manhattan_naive.svg",
        manhattan_adjusted="results/demo/manhattan_adjusted.svg",
        qq_naive="results/demo/qq_naive.svg",
        qq_adjusted="results/demo/qq_adjusted.svg"
    params:
        seed=lambda wc: config["demo"]["seed"],
        samples=lambda wc: config["demo"]["samples"],
        features=lambda wc: config["demo"]["features"],
        pcs=lambda wc: config["demo"]["pcs"]
    shell:
        "pangwasflow demo --outdir results/demo --seed {params.seed} "
        "--samples {params.samples} --features {params.features} --pcs {params.pcs}"
