"""Reading input and extracting clusters, on top of PhyloZoo's I/O.

PhyloZoo parses (e)Newick, so PhyloCass has no parser of its own.  What is
added here is reading *several* trees from one file -- PhyloZoo's
``from_string`` deliberately accepts a single tree -- and turning networks into
the cluster sets Cass consumes.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Iterator

from phylozoo import DirectedPhyNetwork
from phylozoo.core.network.dnetwork.derivations import displayed_trees, subnetwork

__all__ = [
    "read_trees",
    "read_tree_file",
    "read_cluster_file",
    "hardwired_clusters",
    "softwired_clusters",
    "softwired_clusters_via_displayed_trees",
    "displays_trees",
    "displays_trees_via_displayed_trees",
    "displays_partial_trees",
    "clusters_of_trees",
    "per_tree_clusters",
]


#: Byte-order mark.  Windows tooling sprinkles these liberally -- Notepad
#: writes one, and PowerShell prepends one when piping into a program -- and
#: left in place it parses as a taxon and stops ``#`` comment lines from being
#: recognised.  It is never meaningful in Newick, so it is always dropped.
_BOM = "﻿"


def _split_trees(text: str) -> Iterator[str]:
    """Yield each ``;``-terminated tree, dropping blank and ``#``-comment lines."""
    kept = []
    for line in text.replace(_BOM, "").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        kept.append(stripped)
    joined = " ".join(kept)
    for chunk in joined.split(";"):
        if chunk.strip():
            yield chunk.strip() + ";"


def read_trees(text: str) -> list[DirectedPhyNetwork]:
    """Parse every (e)Newick tree in ``text`` via PhyloZoo."""
    return [DirectedPhyNetwork.from_string(chunk) for chunk in _split_trees(text)]


def read_tree_file(path: str | Path) -> list[DirectedPhyNetwork]:
    # utf-8-sig transparently drops a leading byte-order mark if there is one
    return read_trees(Path(path).read_text(encoding="utf-8-sig"))


def read_cluster_file(path: str | Path) -> tuple[set[frozenset], frozenset]:
    """Read one cluster per line, taxa separated by whitespace or commas."""
    clusters: set[frozenset] = set()
    taxa: set[str] = set()
    text = Path(path).read_text(encoding="utf-8-sig").replace(_BOM, "")
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        members = [t for t in line.replace(",", " ").split() if t]
        if members:
            clusters.add(frozenset(members))
            taxa.update(members)
    return clusters, frozenset(taxa)


def hardwired_clusters(network: DirectedPhyNetwork) -> set[frozenset]:
    """Descendant taxon set of every node of ``network``.

    For a tree this is precisely the set of clusters it displays.
    """
    children = {v: [] for v in network.nodes}
    parents = {v: [] for v in network.nodes}
    for u, v in network.edges:
        children[u].append(v)
        parents[v].append(u)

    remaining = {v: len(parents[v]) for v in children}
    queue = [v for v, d in remaining.items() if d == 0]
    order: list = []
    while queue:
        v = queue.pop()
        order.append(v)
        for w in children[v]:
            remaining[w] -= 1
            if remaining[w] == 0:
                queue.append(w)

    below: dict = {}
    for v in reversed(order):
        label = network.get_label(v)
        acc = frozenset({label}) if (label and not children[v]) else frozenset()
        for w in children[v]:
            acc |= below[w]
        below[v] = acc
    return {s for s in below.values() if s}


def softwired_clusters(network: DirectedPhyNetwork) -> set[frozenset]:
    """Every cluster ``network`` represents in the softwired sense.

    Computed blob by blob (see
    :meth:`phylocass.workgraph.WorkGraph.softwired_clusters`), which costs a
    sum of ``2 ** (reticulations in one blob)`` instead of
    ``2 ** (total reticulations)``. On a network whose conflicts are local that
    is the difference between instant and hopeless.

    :func:`softwired_clusters_via_displayed_trees` is the same thing computed
    the obvious way; the tests check the two agree.
    """
    from .workgraph import WorkGraph

    return WorkGraph.from_phylozoo(network).softwired_clusters()


def softwired_clusters_via_displayed_trees(
    network: DirectedPhyNetwork,
) -> set[frozenset]:
    """Softwired clusters straight from PhyloZoo's ``displayed_trees``.

    Exponential in the *total* number of reticulations, so unusable on large
    networks -- but it goes through PhyloZoo rather than any of PhyloCass's own
    machinery, which makes it the reference the fast version is checked against.
    """
    found: set[frozenset] = set()
    for tree in displayed_trees(network):
        found |= hardwired_clusters(tree)
    return found


def per_tree_clusters(
    trees: Iterable[DirectedPhyNetwork], taxa: frozenset | None = None
) -> list[set[frozenset]]:
    """The non-trivial clusters of each tree, kept separate per tree.

    This is the provenance the display-mode search needs: which clusters have
    to end up in one and the same switching.
    """
    trees = list(trees)
    if taxa is None:
        taxa = frozenset().union(*(frozenset(t.taxa) for t in trees)) if trees else frozenset()
    return [
        {c for c in hardwired_clusters(t) if 1 < len(c) < len(taxa)} for t in trees
    ]


def displays_trees(
    network: DirectedPhyNetwork, tree_clusters: Iterable[Iterable[frozenset]]
) -> bool:
    """Does ``network`` display each tree, not merely its clusters?

    For every tree there must be a single switching whose clusters include all
    of it.  Decided blob by blob (see
    :meth:`phylocass.workgraph.WorkGraph.displays`) rather than by enumerating
    every displayed tree.

    :func:`displays_trees_via_displayed_trees` is the same question answered
    the obvious way; the tests check the two agree.
    """
    from .workgraph import WorkGraph

    return WorkGraph.from_phylozoo(network).displays(tree_clusters)


def displays_trees_via_displayed_trees(
    network: DirectedPhyNetwork, tree_clusters: Iterable[Iterable[frozenset]]
) -> bool:
    """As :func:`displays_trees`, straight from PhyloZoo's ``displayed_trees``.

    Exponential in the total number of reticulations; kept as the reference the
    fast version is checked against.
    """
    wanted = [{c for c in want if c} for want in tree_clusters]
    wanted = [w for w in wanted if w]
    if not wanted:
        return True
    available = [hardwired_clusters(t) for t in displayed_trees(network)]
    return all(any(want <= have for have in available) for want in wanted)


def displays_partial_trees(
    network: DirectedPhyNetwork, trees: Iterable[DirectedPhyNetwork]
) -> list[bool]:
    """Does ``network``, restricted to each tree's taxa, display that tree?

    This is the right question when the input trees have different taxon sets:
    a tree on a subset of the taxa can only be displayed by the corresponding
    sub-network.  Comparing raw cluster sets would ask something stronger and
    slightly wrong.

    Restriction is PhyloZoo's ``subnetwork``; the check is then the usual one.
    Returns one verdict per tree, in order.
    """
    verdicts = []
    for tree in trees:
        taxa = sorted(tree.taxa)
        wanted = {c for c in hardwired_clusters(tree) if 1 < len(c) < len(taxa)}
        if not wanted:
            verdicts.append(True)
            continue
        try:
            restricted = subnetwork(network, taxa)
        except Exception:
            verdicts.append(False)
            continue
        available = [hardwired_clusters(t) for t in displayed_trees(restricted)]
        verdicts.append(any(wanted <= have for have in available))
    return verdicts


def clusters_of_trees(
    trees: Iterable[DirectedPhyNetwork],
) -> tuple[set[frozenset], frozenset]:
    """Collect the non-trivial clusters displayed by a collection of trees.

    Returns ``(clusters, taxa)``.  Singletons and the full taxon set are
    dropped: they are representable by anything, and Cass treats them as
    implicit.
    """
    taxa: frozenset = frozenset()
    raw: set[frozenset] = set()
    for tree in trees:
        taxa |= frozenset(tree.taxa)
        raw |= hardwired_clusters(tree)
    return {c for c in raw if 1 < len(c) < len(taxa)}, taxa
