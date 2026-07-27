"""Command-line interface: ``python -m phylocass``."""

from __future__ import annotations

import argparse
import sys

from .cass import CassOptions, cass, cass_from_trees
from .clusters import clusters_of_trees
from .newick import read_trees

__all__ = ["main"]


def _read_clusters(text: str):
    """Read one cluster per line, taxa separated by whitespace or commas."""
    clusters = set()
    taxa = set()
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        members = [t for t in line.replace(",", " ").split() if t]
        if members:
            clusters.add(frozenset(members))
            taxa.update(members)
    return clusters, frozenset(taxa)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="phylocass",
        description=(
            "Build a rooted phylogenetic network of minimum level representing "
            "all clusters of a set of rooted trees (the Cass algorithm)."
        ),
    )
    p.add_argument(
        "input",
        nargs="?",
        default="-",
        help="input file, or '-' for stdin (default)",
    )
    p.add_argument(
        "--format",
        choices=("newick", "clusters"),
        default="newick",
        help="'newick': one rooted tree per ';'. 'clusters': one cluster per line.",
    )
    p.add_argument(
        "--max-level",
        type=int,
        default=4,
        help="highest level to attempt (default: 4; the paper proves exactness up to 2)",
    )
    p.add_argument(
        "--time-limit",
        type=float,
        default=None,
        help="seconds allowed per biconnected component (default: unlimited)",
    )
    p.add_argument(
        "--max-networks",
        type=int,
        default=20000,
        help="cap on intermediate networks per subproblem (default: 20000)",
    )
    p.add_argument(
        "--show-clusters",
        action="store_true",
        help="also print the input clusters that were collected",
    )
    p.add_argument("--quiet", "-q", action="store_true", help="print only the network")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    text = sys.stdin.read() if args.input == "-" else open(args.input, encoding="utf-8").read()

    options = CassOptions(
        max_level=args.max_level,
        time_limit=args.time_limit,
        max_networks=args.max_networks,
    )

    if args.format == "newick":
        trees = read_trees(text)
        if not trees:
            print("error: no trees found in input", file=sys.stderr)
            return 2
        clusters, taxa = clusters_of_trees(trees)
        result = cass_from_trees(trees, options)
    else:
        clusters, taxa = _read_clusters(text)
        if not clusters:
            print("error: no clusters found in input", file=sys.stderr)
            return 2
        result = cass(clusters, taxa, options)

    if args.show_clusters and not args.quiet:
        print(f"# {len(clusters)} clusters on {len(taxa)} taxa", file=sys.stderr)
        for c in sorted(clusters, key=lambda s: (len(s), sorted(map(str, s)))):
            print("#   {" + ", ".join(sorted(map(str, c))) + "}", file=sys.stderr)

    if not args.quiet:
        print(
            f"# taxa={len(result.taxa)} clusters={len(result.clusters)} "
            f"level={result.level} reticulations={result.reticulation_number} "
            f"verified={result.represents_input()}",
            file=sys.stderr,
        )

    print(result.to_enewick())
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
