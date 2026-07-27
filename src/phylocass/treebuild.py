"""Building the unique tree that represents a compatible set of clusters."""

from __future__ import annotations

from itertools import combinations
from typing import Iterable

from .clusters import is_compatible, round_up
from .workgraph import WorkGraph

__all__ = ["build_tree", "add_root_edge", "graft"]


def build_tree(
    clusters: Iterable[frozenset],
    blocks: Iterable[frozenset],
    return_map: bool = False,
):
    """The unique tree on ``blocks`` representing exactly ``clusters``, or ``None``.

    ``blocks`` is the current taxon set, each block acting as a single taxon;
    every cluster must be a union of blocks.  Returns ``None`` when the
    clusters are not pairwise compatible, i.e. when no such tree exists.

    The tree is the Hasse diagram of the containment order, so every non-root
    node is exactly one input cluster: the tree represents the input clusters
    and nothing else, and has no contraction-safe edges.
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

    tree = WorkGraph()
    block_set = set(blocks)
    ordered = sorted(node_sets, key=lambda s: (len(s), sorted(map(str, s))))
    node_of: dict[frozenset, int] = {
        s: tree.new_node(block=s if s in block_set else None) for s in ordered
    }

    # The root is always its own node.  It must not be looked up through
    # ``node_of``: with a single block that block *is* the whole taxon set, and
    # reusing its node would wire the root to itself.
    root = tree.new_node()

    for s in ordered:  # ordered ascending by size, so the first hit is the parent
        parent = next((t for t in ordered if len(t) > len(s) and s < t), None)
        tree.add_edge(root if parent is None else node_of[parent], node_of[s])

    node_of.setdefault(ground, root)

    if return_map:
        return tree, node_of
    return tree


def add_root_edge(tree: WorkGraph) -> None:
    """Prepend a dummy root above the current root.

    Cass needs this so a reticulation edge can also be hung from above the old
    root; the extra node is suppressed again on output.
    """
    old = tree.root
    tree.add_edge(tree.new_node(), old)


def graft(host: WorkGraph, at: int, guest: WorkGraph, identify: dict[frozenset, int]) -> None:
    """Merge ``guest`` into ``host``.

    ``guest``'s root becomes ``at``; a guest leaf whose block appears in
    ``identify`` becomes the corresponding host node instead of a fresh one.
    Every other guest node is copied.
    """
    mapping: dict[int, int] = {guest.root: at}
    for v in guest.nodes():
        if v in mapping:
            continue
        block = guest.block.get(v)
        if block is not None and block in identify:
            mapping[v] = identify[block]
        else:
            mapping[v] = host.new_node(block=block)
    for u, v in guest.edges():
        host.add_edge(mapping[u], mapping[v])
