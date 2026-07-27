"""Build a network from input trees with Cass and draw it.

Plots each input tree alongside the network Cass builds from their clusters,
using PhyloZoo's matplotlib backend.  Reticulation edges are drawn in red.

Needs PhyloZoo's plotting extra::

    pip install "phylozoo[viz]"

Usage::

    python examples/plot_network.py                              # built-in example
    python examples/plot_network.py examples/double_conflict.newick
    python examples/plot_network.py my_trees.newick -o out.png
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from phylozoo.viz import plot

from phylocass import cass_from_trees, read_tree_file, read_trees

DEFAULT_TREES = """
    (((a,b),c),(d,e));
    (((a,c),b),(d,e));
"""


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("input", nargs="?", help="Newick file; omit for a built-in example")
    parser.add_argument("-o", "--out", default="network.png", help="image to write")
    parser.add_argument("--show", action="store_true", help="also open a window")
    args = parser.parse_args()

    trees = read_tree_file(args.input) if args.input else read_trees(DEFAULT_TREES)

    result = cass_from_trees(trees)
    print(f"taxa               {len(result.taxa)}")
    print(f"clusters           {len(result.clusters)}")
    print(f"level              {result.level}")
    print(f"reticulations      {result.reticulation_number}")
    print(f"displays all input {result.represents_input()}")
    print(result.to_enewick())

    # one panel per input tree, then a wider one for the network
    n = len(trees)
    fig, axes = plt.subplots(
        1, n + 1, figsize=(3.2 * n + 5, 4.5), gridspec_kw={"width_ratios": [1] * n + [1.8]}
    )
    axes = list(axes) if n + 1 > 1 else [axes]

    for i, (tree, ax) in enumerate(zip(trees, axes), start=1):
        plot(tree, ax=ax)
        ax.set_title(f"input tree {i}", fontsize=11)

    plot(result.network, ax=axes[-1])
    axes[-1].set_title(
        f"Cass network — level {result.level}, "
        f"{result.reticulation_number} reticulation"
        f"{'s' if result.reticulation_number != 1 else ''}",
        fontsize=11,
    )

    for ax in axes:
        ax.set_axis_off()

    fig.suptitle("Conflicting clusters combined into one network", fontsize=13)
    fig.tight_layout()
    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    print(f"\nwrote {Path(args.out).resolve()}")

    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
