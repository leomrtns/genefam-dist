"""Command line entry point for the Python tree-signature workflow."""

import argparse
import csv
from pathlib import Path

from .analysis import LENGTH_METRICS, analyze_files, generate_references
from . import _treesignalc


def _mapping(path):
    mapping = {}
    with Path(path).open(encoding="utf-8", newline="") as handle:
        for index, row in enumerate(csv.reader(handle, delimiter="\t"), 1):
            if not row or row[0].startswith("#"):
                continue
            if len(row) != 2 or not all(row):
                raise ValueError(f"{path}:{index}: expected two nonempty tab-separated fields (gene, species)")
            if row[0] in mapping:
                raise ValueError(f"{path}:{index}: duplicate gene label {row[0]!r}")
            mapping[row[0]] = row[1]
    return mapping


def main(argv=None):
    parser = argparse.ArgumentParser(description="Generate reference trees and compare gene-family trees.")
    commands = parser.add_subparsers(dest="command", required=True)
    analyze = commands.add_parser("analyze", help="compare one Newick gene tree per file")
    analyze.add_argument("--references", required=True, type=Path)
    analyze.add_argument("--genes", required=True, nargs="+", help="gene-tree paths or quoted glob patterns")
    analyze.add_argument("--output", required=True, type=Path)
    analyze.add_argument("--metrics", nargs="+", choices=(*LENGTH_METRICS, "topology"), default=LENGTH_METRICS)
    analyze.add_argument("--normalize", action="store_true", help="remove length scale (or rescale topology counts)")
    analyze.add_argument("--mapping", type=Path, help="two-column gene/species TSV, without header")
    analyze.add_argument("--aggregation", choices=("min", "mean"), default="min")
    references = commands.add_parser("references", help="generate a reference SPR chain using the native API")
    references.add_argument("--species", required=True, type=Path, help="one species label per line")
    references.add_argument("--count", type=int, default=20)
    references.add_argument("--spr", type=int, default=1)
    references.add_argument("--seed", type=int, default=0)
    references.add_argument("--branch-length", type=float)
    references.add_argument("--output", required=True, type=Path)
    args = parser.parse_args(argv)
    try:
        if args.command == "references":
            if args.output.resolve() == args.species.resolve():
                raise ValueError("output must not overwrite the species-label input file")
            labels = [line.strip() for line in args.species.read_text(encoding="utf-8").splitlines()
                      if line.strip() and not line.lstrip().startswith("#")]
            trees = generate_references(labels, count=args.count, n_spr=args.spr,
                                        seed=args.seed, branch_length=args.branch_length)
            args.output.write_text(trees.newick + "\n", encoding="utf-8")
        else:
            protected = [args.references, *([args.mapping] if args.mapping else [])]
            if args.output.resolve() in {path.resolve() for path in protected}:
                raise ValueError("output must not overwrite the references or mapping input file")
            analyze_files(args.genes, args.references.read_text(encoding="utf-8"), args.output,
                          metrics=args.metrics, normalize=args.normalize,
                          gene_to_species=_mapping(args.mapping) if args.mapping else None,
                          aggregation=args.aggregation)
    except (OSError, ValueError, TypeError, OverflowError, _treesignalc.error) as error:
        parser.error(str(error))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
