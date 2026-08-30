from __future__ import annotations

import argparse

from .benchmark import run_demo
from .inputs import prepare_inputs
from .real_analysis import run_prepared_paths


def main() -> None:
    parser = argparse.ArgumentParser(prog="pangwasflow")
    subparsers = parser.add_subparsers(dest="command", required=True)

    demo = subparsers.add_parser("demo", help="run the deterministic population-structure benchmark")
    demo.add_argument("--outdir", default="results/demo")
    demo.add_argument("--seed", type=int, default=42)
    demo.add_argument("--samples", type=int, default=400)
    demo.add_argument("--features", type=int, default=120)
    demo.add_argument("--pcs", type=int, default=2)

    prepare = subparsers.add_parser("prepare", help="normalize and QC real binary genomic feature inputs")
    prepare.add_argument("--features", required=True, help="sample-by-feature matrix, Roary/Panaroo gene table, or Rtab file")
    prepare.add_argument("--metadata", required=True, help="phenotype metadata TSV/CSV")
    prepare.add_argument("--outdir", default="results/prepared")
    prepare.add_argument("--format", choices=["matrix", "roary", "panaroo", "gene-pa", "rtab"], default="matrix")
    prepare.add_argument("--sample-id-column", default="sample_id")
    prepare.add_argument("--label-column", default="label")
    prepare.add_argument("--min-prevalence", type=float, default=0.01)
    prepare.add_argument("--max-prevalence", type=float, default=0.99)
    prepare.add_argument("--min-samples", type=int, default=20)

    analyze = subparsers.add_parser("analyze", help="run baseline and structure-adjusted GWAS on prepared binary inputs")
    analyze.add_argument("--features", required=True, help="prepared features_qc.tsv")
    analyze.add_argument("--metadata", required=True, help="prepared metadata_qc.tsv")
    analyze.add_argument("--outdir", default="results/gwas")
    analyze.add_argument("--pcs", type=int, default=2, help="number of PCA or MDS structure components")
    analyze.add_argument(
        "--structure-distance",
        default=None,
        help="optional square TSV genomic distance matrix; when supplied, classical MDS replaces feature-derived PCA",
    )

    args = parser.parse_args()
    if args.command == "demo":
        run_demo(args.outdir, seed=args.seed, samples=args.samples, features=args.features, pcs=args.pcs)
    elif args.command == "prepare":
        prepare_inputs(
            args.features,
            args.metadata,
            outdir=args.outdir,
            feature_format=args.format,
            sample_id_column=args.sample_id_column,
            label_column=args.label_column,
            min_prevalence=args.min_prevalence,
            max_prevalence=args.max_prevalence,
            min_samples=args.min_samples,
        )
    elif args.command == "analyze":
        run_prepared_paths(
            args.features,
            args.metadata,
            outdir=args.outdir,
            pcs=args.pcs,
            distance_matrix_path=args.structure_distance,
        )


if __name__ == "__main__":
    main()
