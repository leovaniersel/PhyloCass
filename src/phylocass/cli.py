"""Command-line interface: ``python -m phylocass`` or ``phylocass``."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from phylozoo.utils.exceptions.base import PhyloZooError

from . import __version__
from .cass import CassOptions, cass, cass_from_trees
from .io import read_cluster_file, read_trees

__all__ = ["main"]


def _phylozoo_version() -> str:
    try:
        import phylozoo

        return getattr(phylozoo, "__version__", "?")
    except Exception:  # pragma: no cover - phylozoo is a hard dependency
        return "not installed"


def _read_stdin() -> str:
    """Read stdin as UTF-8, whatever the platform's default encoding is.

    ``sys.stdin.read()`` decodes with the locale encoding, which on Windows is
    typically cp1252 -- so a UTF-8 taxon name, or the byte-order mark that
    PowerShell prepends when piping, would come through mangled.  Decoding the
    raw bytes ourselves makes piped input behave the same everywhere.
    """
    buffer = getattr(sys.stdin, "buffer", None)
    if buffer is None:  # stdin replaced, e.g. by a test harness
        return sys.stdin.read()
    data = buffer.read()
    try:
        return data.decode("utf-8-sig")
    except UnicodeDecodeError:
        return data.decode("latin-1")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="phylocass",
        description=(
            "Build a rooted phylogenetic network of minimum level representing "
            "all clusters of a set of rooted trees (the Cass algorithm)."
        ),
    )
    p.add_argument(
        "--version",
        action="version",
        version=f"phylocass {__version__} (phylozoo {_phylozoo_version()})",
    )
    p.add_argument("input", nargs="?", default="-", help="input file, or '-' for stdin")
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
        help="seconds per conflicting component, shared across levels (default: none)",
    )
    p.add_argument(
        "--max-networks",
        type=int,
        default=20000,
        help="cap on intermediate networks per subproblem (default: 20000)",
    )
    p.add_argument(
        "--z-closure",
        action="store_true",
        help=(
            "accept input trees with different but overlapping taxon sets, by "
            "completing their partial clusters with Z-closure first (heuristic)"
        ),
    )
    p.add_argument(
        "--display-trees",
        action="store_true",
        help=(
            "require the network to display the input trees, not just their "
            "clusters; makes the reticulation number a heuristic upper bound "
            "on the hybridization number (Newick input only)"
        ),
    )
    p.add_argument(
        "--out",
        default=None,
        help="write the network here via PhyloZoo (.enewick/.nwk, or .dot)",
    )
    p.add_argument(
        "--plot",
        metavar="IMAGE",
        default=None,
        help=(
            "draw the input trees and the resulting network to an image file "
            '(.png, .pdf, .svg); needs pip install "phylocass[viz]"'
        ),
    )
    p.add_argument(
        "--plot-dpi",
        type=int,
        default=150,
        help="resolution for --plot (default: 150)",
    )
    p.add_argument(
        "--show",
        action="store_true",
        help="open the drawing in a window as well as, or instead of, saving it",
    )
    p.add_argument(
        "--show-clusters", action="store_true", help="list the input clusters on stderr"
    )
    p.add_argument("--quiet", "-q", action="store_true", help="print only the network")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    options = CassOptions(
        max_level=args.max_level,
        time_limit=args.time_limit,
        max_networks=args.max_networks,
        display_trees=args.display_trees,
        z_closure=args.z_closure,
    )

    if args.format == "clusters":
        if args.display_trees:
            print(
                "error: --display-trees needs tree input; a bare cluster list does "
                "not say which cluster came from which tree",
                file=sys.stderr,
            )
            return 2
        if args.z_closure:
            print(
                "error: --z-closure needs tree input; a bare cluster list does "
                "not say which taxon set each cluster is known on",
                file=sys.stderr,
            )
            return 2

    try:
        if args.format == "newick":
            text = (
                _read_stdin()
                if args.input == "-"
                else Path(args.input).read_text(encoding="utf-8-sig")
            )
            trees = read_trees(text)
            if not trees:
                print("error: no trees found in input", file=sys.stderr)
                return 2
            result = cass_from_trees(trees, options)
        else:
            if args.input == "-":
                print("error: --format clusters needs a file, not stdin", file=sys.stderr)
                return 2
            trees = []
            clusters, taxa = read_cluster_file(args.input)
            if not clusters:
                print("error: no clusters found in input", file=sys.stderr)
                return 2
            result = cass(clusters, taxa, options)
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    except FileNotFoundError:
        print(f"error: no such file: {args.input}", file=sys.stderr)
        return 2
    except PhyloZooError as exc:
        # malformed (e)Newick and the like: report it, do not dump a traceback
        print(f"error: could not read {args.input}: {exc}", file=sys.stderr)
        return 2
    except RuntimeError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    if args.show_clusters and not args.quiet:
        print(
            f"# {len(result.clusters)} clusters on {len(result.taxa)} taxa",
            file=sys.stderr,
        )
        for c in sorted(result.clusters, key=lambda s: (len(s), sorted(map(str, s)))):
            print("#   {" + ", ".join(sorted(map(str, c))) + "}", file=sys.stderr)

    if not args.quiet and result.z_closure is not None:
        z = result.z_closure
        print(
            f"# z-closure: {z.partial_total} clusters known, {len(z.clusters)} "
            f"complete, {z.dropped} dropped as still partial "
            f"({z.rounds} rounds{', hit the cap' if z.hit_limit else ''})",
            file=sys.stderr,
        )

    if not args.quiet:
        displays = result.displays_input_trees()
        extra = "" if displays is None else f" displays-trees={displays}"
        print(
            f"# taxa={len(result.taxa)} clusters={len(result.clusters)} "
            f"level={result.level} reticulations={result.reticulation_number} "
            f"verified={result.represents_input()}{extra}",
            file=sys.stderr,
        )

    if args.out:
        result.network.save(args.out, overwrite=True)

    if args.plot or args.show:
        from .plot import PlotUnavailable, save_plot

        try:
            written = save_plot(
                result,
                args.plot or "phylocass.png",
                trees=trees,
                dpi=args.plot_dpi,
                show=args.show,
            )
        except PlotUnavailable as exc:
            print(f"error: {exc}", file=sys.stderr)
            return 2
        if not args.quiet:
            print(f"# wrote {written}", file=sys.stderr)

    print(result.to_enewick())
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
