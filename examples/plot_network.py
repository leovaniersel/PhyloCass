"""Build a network from input trees with Cass and draw it.

A thin wrapper around ``phylocass.plot.save_plot``; the CLI's ``--plot`` does
the same thing without writing any Python::

    phylocass examples/double_conflict.newick --plot network.png

Needs PhyloZoo's plotting extra::

    pip install "phylocass[viz]"

Usage::

    python examples/plot_network.py                              # built-in example
    python examples/plot_network.py examples/double_conflict.newick
    python examples/plot_network.py my_trees.newick -o out.png --show
"""

from __future__ import annotations

import argparse

from phylocass import CassOptions, cass_from_trees, read_tree_file, read_trees
from phylocass.plot import save_plot

DEFAULT_TREES = """
    (((a,b),c),(d,e));
    (((a,c),b),(d,e));
"""


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("input", nargs="?", help="Newick file; omit for a built-in example")
    parser.add_argument("-o", "--out", default="network.png", help="image to write")
    parser.add_argument("--show", action="store_true", help="also open a window")
    parser.add_argument(
        "--display-trees",
        action="store_true",
        help="require the network to display the input trees, not just their clusters",
    )
    args = parser.parse_args()

    trees = read_tree_file(args.input) if args.input else read_trees(DEFAULT_TREES)
    result = cass_from_trees(trees, CassOptions(display_trees=args.display_trees))

    print(f"taxa               {len(result.taxa)}")
    print(f"clusters           {len(result.clusters)}")
    print(f"level              {result.level}")
    print(f"reticulations      {result.reticulation_number}")
    print(f"displays all input {result.represents_input()}")
    print(f"displays the trees {result.displays_input_trees()}")
    print(result.to_enewick())

    written = save_plot(result, args.out, trees=trees, show=args.show)
    print(f"\nwrote {written}")


if __name__ == "__main__":
    main()
