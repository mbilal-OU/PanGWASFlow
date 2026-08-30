configfile: "config/config.yaml"

include: "workflow/rules/demo.smk"

rule all:
    input:
        "results/demo/summary.json",
        "results/demo/features.tsv",
        "results/demo/metadata.tsv",
        "results/demo/association_naive.tsv"
