from __future__ import annotations

import argparse

from .benchmark import run_demo


def main() -> None:
    parser = argparse.ArgumentParser(prog="pangwasflow")
    subparsers = parser.add_subparsers(dest="command", required=True)
    demo = subparsers.add_parser("demo", help="run the deterministic teaching benchmark")
    demo.add_argument("--outdir", default="results/demo")
    demo.add_argument("--seed", type=int, default=42)
    demo.add_argument("--samples", type=int, default=240)
    demo.add_argument("--features", type=int, default=120)
    args = parser.parse_args()

    if args.command == "demo":
        run_demo(args.outdir, seed=args.seed, samples=args.samples, features=args.features)


if __name__ == "__main__":
    main()
