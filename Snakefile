configfile: "config/config.yaml"

include: "workflow/rules/demo.smk"

rule all:
    input:
        "results/demo/summary.json",
        "results/demo/association_naive.tsv",
        "results/demo/association_adjusted.tsv",
        "results/demo/manhattan_naive.svg",
        "results/demo/manhattan_adjusted.svg",
        "results/demo/qq_naive.svg",
        "results/demo/qq_adjusted.svg"
