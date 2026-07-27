"""Building the unique tree that represents a compatible set of clusters."""

from __future__ import annotations

from itertools import combinations
from typing import Iterable

from .clusters import is_compatible, round_up
from .network import Network

__all__ = ["build_tree", "add_root_edge", "graft"]


def build_tree(
    clusters: Iterable[frozenset],
    blocks: Iterable[frozenset],
    return_map: bool = False,
):
    """The unique tree on ``blocks`` representing exactly ``clusters``, or ``None``.

    ``blocks`` is the current taxon set (each block acts as one taxon); each
    cluster must be a union of blocks.  Returns ``None`` when the clusters are
    not pairwise compatible, i.e. when no such tree exists.

    The tree is the Hasse diagram of the containment order: every non-root node
    is one input cluster, so the tree represents the input clusters and nothing
    else, and it contains no contraction-safe edges.
    """
    blocks = list(blocks)
    if not blocks:
        return (None, {}) if return_map else None
    ground: frozenset = frozenset().union(*blocks)

    node_sets: set[frozenset] = set(blocks)
    for c in clusters:
        c = round_up(c, blocks)
        if c and c != ground:
            node_sets.add(c)

    for a, b in combinations(node_sets, 2):
        if not is_compatible(a, b):
            return (None, {}) if return_map else None

    tree = Network()
    block_set = set(blocks)
    ordered = sorted(node_sets, key=lambda s: (len(s), sorted(map(str, s))))
    node_of: dict[frozenset, int] = {}
    for s in ordered:
        node_of[s] = tree.new_node(label=s if s in block_set else None)

    # The root is always a separate node.  It must not be looked up through
    # ``node_of``: with a single block that block *is* the whole taxon set, and
    # reusing its node would wire the root to itself.
    root = tree.new_node()

    # parent = smallest strictly larger node set (ordered is ascending by size)
    for s in ordered:
        parent = None
        for t in ordered:
            if len(t) > len(s) and s < t:
                parent = t
                break
        tree.add_edge(root if parent is None else node_of[parent], node_of[s])

    node_of.setdefault(ground, root)

    if return_map:
        return tree, node_of
    return tree


def add_root_edge(tree: Network) -> None:
    """Prepend a dummy root above the current root.

    Cass needs this so that a reticulation edge can also be hung from above the
    old root; the extra node is suppressed again on output.
    """
    old = tree.root
    new = tree.new_node()
    tree.add_edge(new, old)


def graft(host: Network, at: int, guest: Network, identify: dict[frozenset, int]) -> None:
    """Merge ``guest`` into ``host``.

    ``guest``'s root becomes ``at``; a guest leaf whose label appears in
    ``identify`` becomes the corresponding host node instead of a fresh one.
    Remaining guest nodes are copied.
    """
    mapping: dict[int, int] = {guest.root: at}
    for v in guest.nodes():
        if v in mapping:
            continue
        lbl = guest.label.get(v)
        if lbl is not None and lbl in identify:
            mapping[v] = identify[lbl]
        else:
            mapping[v] = host.new_node(label=lbl)
    for u, v in guest.edges():
        host.add_edge(mapping[u], mapping[v])
