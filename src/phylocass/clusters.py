"""Clusters, compatibility, incompatibility graphs, ST-sets.

Throughout PhyloCass a *taxon* is an arbitrary hashable label (normally a
string) and a *cluster* is a ``frozenset`` of taxa.  Following van Iersel et
al. (2010), a cluster is a *proper* subset of the full taxon set X.

During the Cass recursion taxa get collapsed into groups.  Rather than
inventing composite taxon names, PhyloCass keeps every cluster expressed in
terms of the *original* taxa and tracks the current level of collapsing with a
separate partition of the taxon set into *blocks*.  A block is itself a
``frozenset`` of original taxa and plays the role of a single taxon.  This
keeps ``C|S``, ``C \\ S`` and the "does this network represent C?" test free of
any translation layer.
"""

from __future__ import annotations

from itertools import combinations
from typing import Iterable, Sequence

Cluster = frozenset
Block = frozenset

__all__ = [
    "is_compatible",
    "is_separated",
    "restrict",
    "remove_taxa",
    "incompatibility_graph",
    "nontrivial_components",
    "is_st_set",
    "maximal_st_sets",
    "maximal_unseparated_sets",
    "unseparated_closure",
    "collapse",
    "round_up",
    "clusters_of_trees",
]


def is_compatible(a: frozenset, b: frozenset) -> bool:
    """Two clusters are compatible if they are disjoint or one contains the other."""
    return a.isdisjoint(b) or a <= b or b <= a


def is_separated(x: frozenset, clusters: Iterable[frozenset]) -> bool:
    """A taxon set ``x`` is separated by ``clusters`` if some cluster is incompatible with it."""
    return any(not is_compatible(c, x) for c in clusters)


def restrict(clusters: Iterable[frozenset], s: frozenset) -> set[frozenset]:
    """``C|S``: restrict every cluster to ``s``, dropping empty and improper results."""
    out = set()
    for c in clusters:
        r = c & s
        if r and r != s:
            out.add(r)
    return out


def remove_taxa(clusters: Iterable[frozenset], s: frozenset, ground: frozenset) -> set[frozenset]:
    """``C \\ S``: remove the taxa in ``s`` from every cluster.

    ``ground`` is the taxon set *after* the removal; clusters that become empty
    or that swallow the whole remaining taxon set are dropped, since clusters
    are by definition proper non-empty subsets.
    """
    out = set()
    for c in clusters:
        r = c - s
        if r and r != ground:
            out.add(r)
    return out


def incompatibility_graph(clusters: Sequence[frozenset]) -> list[set[int]]:
    """Adjacency list of ``IG(C)``, indexed by position in ``clusters``."""
    adj: list[set[int]] = [set() for _ in clusters]
    for i, j in combinations(range(len(clusters)), 2):
        if not is_compatible(clusters[i], clusters[j]):
            adj[i].add(j)
            adj[j].add(i)
    return adj


def nontrivial_components(clusters: Sequence[frozenset]) -> list[list[frozenset]]:
    """Connected components of ``IG(C)`` that contain more than one cluster."""
    adj = incompatibility_graph(clusters)
    seen = [False] * len(clusters)
    components = []
    for start in range(len(clusters)):
        if seen[start]:
            continue
        stack = [start]
        seen[start] = True
        comp = []
        while stack:
            v = stack.pop()
            comp.append(v)
            for w in adj[v]:
                if not seen[w]:
                    seen[w] = True
                    stack.append(w)
        if len(comp) > 1:
            components.append([clusters[i] for i in sorted(comp)])
    return components


def is_st_set(s: frozenset, clusters: Iterable[frozenset], ground: frozenset) -> bool:
    """Is ``s`` an ST-set (strict tree set) w.r.t. ``clusters``?

    ``s != ground``, ``s`` is not separated, and ``C|s`` is pairwise compatible.
    """
    if s == ground or not s:
        return False
    clusters = list(clusters)
    if is_separated(s, clusters):
        return False
    inside = list(restrict(clusters, s))
    return all(is_compatible(a, b) for a, b in combinations(inside, 2))


def maximal_st_sets(clusters: Iterable[frozenset], blocks: Iterable[frozenset]) -> list[frozenset]:
    """Partition ``blocks`` into the maximal ST-sets w.r.t. ``clusters``.

    Merging two blocks at a time suffices.  Inside a maximal ST-set ``S`` the
    restriction ``C|S`` is compatible, so it is a tree, and any two sibling
    blocks of that tree can be merged without leaving the ST-set family --
    so greedy merging never stalls before reaching the maximal sets.
    """
    clusters = list(clusters)
    parts = list(blocks)
    ground: frozenset = frozenset().union(*parts) if parts else frozenset()

    merged = True
    while merged:
        merged = False
        for i, j in combinations(range(len(parts)), 2):
            union = parts[i] | parts[j]
            if is_st_set(union, clusters, ground):
                parts = [p for idx, p in enumerate(parts) if idx not in (i, j)]
                parts.append(union)
                merged = True
                break
    return parts


def unseparated_closure(
    seed: frozenset, clusters: Iterable[frozenset], ground: frozenset
) -> frozenset | None:
    """The smallest unseparated superset of ``seed``, or ``None`` if that is ``ground``.

    Whenever a cluster ``C`` is incompatible with the current set ``S``, every
    unseparated superset of ``S`` must also contain ``C``, so absorbing ``C``
    keeps the result minimal.
    """
    clusters = list(clusters)
    s = frozenset(seed)
    changed = True
    while changed:
        changed = False
        for c in clusters:
            if not is_compatible(c, s):
                s = s | c
                if s == ground:
                    return None
                changed = True
    return None if s == ground else s


def maximal_unseparated_sets(
    clusters: Iterable[frozenset], blocks: Iterable[frozenset]
) -> list[frozenset]:
    """Partition ``blocks`` into the maximal proper subsets not separated by ``clusters``.

    This is the weaker notion used in Step 1 of the decomposition: only the
    "unseparated" half of the ST-set condition is required.

    Unlike ST-sets, a subset of an unseparated set need not itself be
    unseparated, so plain pairwise merging can stall.  Each candidate merge is
    therefore grown to the *smallest* unseparated superset of the pair first,
    which is safe precisely because unseparated sets are closed under union of
    intersecting members.
    """
    clusters = list(clusters)
    parts = list(blocks)
    ground: frozenset = frozenset().union(*parts) if parts else frozenset()

    merged = True
    while merged:
        merged = False
        for i, j in combinations(range(len(parts)), 2):
            grown = unseparated_closure(parts[i] | parts[j], clusters, ground)
            if grown is None:
                continue
            # the merge target must be a union of whole blocks
            changed = True
            while changed and grown is not None:
                changed = False
                for p in parts:
                    if p & grown and not p <= grown:
                        grown = grown | p
                        if grown == ground:
                            grown = None
                            break
                        changed = True
            if grown is None or is_separated(grown, clusters):
                continue
            parts = [p for p in parts if not p <= grown]
            parts.append(grown)
            merged = True
            break
    return parts


def round_up(cluster: frozenset, blocks: Iterable[frozenset]) -> frozenset:
    """Express ``cluster`` in terms of ``blocks``: the union of every block it meets."""
    out: frozenset = frozenset()
    for b in blocks:
        if b & cluster:
            out |= b
    return out


def collapse(
    clusters: Iterable[frozenset], blocks: Iterable[frozenset]
) -> tuple[list[frozenset], dict[frozenset, list[frozenset]], set[frozenset]]:
    """``Collapse(C)`` of the paper, in block form.

    Returns ``(new_blocks, decollapse_map, new_clusters)`` where
    ``decollapse_map`` sends each new block to the old blocks it swallowed and
    ``new_clusters`` is ``C`` rewritten over the new blocks (clusters that end
    up inside a single block become trivial and are dropped).
    """
    clusters = list(clusters)
    blocks = list(blocks)
    new_blocks = maximal_st_sets(clusters, blocks)
    ground = frozenset().union(*new_blocks) if new_blocks else frozenset()

    mapping = {nb: [b for b in blocks if b <= nb] for nb in new_blocks}

    new_clusters = set()
    for c in clusters:
        r = round_up(c, new_blocks)
        if r and r != ground and not any(r == nb for nb in new_blocks):
            new_clusters.add(r)
    return new_blocks, mapping, new_clusters


def clusters_of_trees(trees: Iterable["object"]) -> tuple[set[frozenset], frozenset]:
    """Collect every (non-trivial) cluster displayed by a collection of trees.

    Each tree must expose ``.clusters()`` and ``.taxa()`` (see
    :mod:`phylocass.network`).  Returns ``(clusters, taxa)`` with the taxon set
    being the union over all trees; clusters equal to the whole taxon set or to
    a singleton are excluded, matching the convention that Cass adds singletons
    back implicitly.
    """
    taxa: frozenset = frozenset()
    raw: set[frozenset] = set()
    for t in trees:
        taxa |= t.taxa()
        raw |= t.clusters()
    clusters = {c for c in raw if 1 < len(c) < len(taxa)}
    return clusters, taxa
