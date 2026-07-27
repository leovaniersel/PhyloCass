"""Clusters, compatibility, incompatibility graphs, ST-sets.

This is the combinatorial core of Cass, and the part PhyloZoo has no
equivalent for -- it supplies networks, splits and partitions, but not the
cluster machinery the algorithm is built on.

A *cluster* is a ``frozenset`` of taxa and, following van Iersel et al.
(2010), a proper subset of the full taxon set X.

The Cass recursion collapses taxa into groups over and over.  Rather than
minting composite taxon names, every cluster here stays written in terms of
the *original* taxa, and the current level of collapsing is tracked separately
as a partition into **blocks** -- a block being a ``frozenset`` of original
taxa that currently acts as one taxon.  Restriction, deletion and the
"does this network represent C?" test then need no translation layer.
"""

from __future__ import annotations

from itertools import combinations
from typing import Iterable, Sequence

from phylozoo.core.primitives.partition import Partition

__all__ = [
    "canonical",
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
    "block_partition",
]


def _key(s: frozenset) -> tuple:
    return (len(s), sorted(map(str, s)))


def canonical(sets: Iterable[frozenset]) -> list[frozenset]:
    """Order taxon sets deterministically.

    Cass explores taxa and blocks in whatever order it is handed them, and
    several equally valid networks usually exist, so iterating a ``set``
    directly would make the output depend on Python's hash randomisation.
    Sorting at every point that drives the search makes runs reproducible.
    """
    return sorted(sets, key=_key)


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


def remove_taxa(
    clusters: Iterable[frozenset], s: frozenset, ground: frozenset
) -> set[frozenset]:
    """``C \\ S``: remove the taxa in ``s`` from every cluster.

    ``ground`` is the taxon set *after* the removal; clusters that become empty
    or swallow the whole remaining taxon set are dropped, since clusters are by
    definition proper non-empty subsets.
    """
    out = set()
    for c in clusters:
        r = c - s
        if r and r != ground:
            out.add(r)
    return out


def round_up(cluster: frozenset, blocks: Iterable[frozenset]) -> frozenset:
    """Express ``cluster`` over ``blocks``: the union of every block it meets."""
    out: frozenset = frozenset()
    for b in blocks:
        if b & cluster:
            out |= b
    return out


def block_partition(blocks: Iterable[frozenset]) -> Partition:
    """Wrap blocks as a PhyloZoo :class:`~phylozoo.core.primitives.partition.Partition`.

    Handy for interoperating with the rest of PhyloZoo; the algorithm itself
    works with plain lists of blocks.
    """
    return Partition([set(b) for b in blocks])


# ----------------------------------------------------------------------
# incompatibility graph
# ----------------------------------------------------------------------


def incompatibility_graph(clusters: Sequence[frozenset]) -> list[set[int]]:
    """Adjacency list of ``IG(C)``, indexed by position in ``clusters``."""
    adjacency: list[set[int]] = [set() for _ in clusters]
    for i, j in combinations(range(len(clusters)), 2):
        if not is_compatible(clusters[i], clusters[j]):
            adjacency[i].add(j)
            adjacency[j].add(i)
    return adjacency


def nontrivial_components(clusters: Sequence[frozenset]) -> list[list[frozenset]]:
    """Connected components of ``IG(C)`` containing more than one cluster."""
    adjacency = incompatibility_graph(clusters)
    seen = [False] * len(clusters)
    components: list[list[frozenset]] = []
    for start in range(len(clusters)):
        if seen[start]:
            continue
        seen[start] = True
        stack, group = [start], []
        while stack:
            v = stack.pop()
            group.append(v)
            for w in adjacency[v]:
                if not seen[w]:
                    seen[w] = True
                    stack.append(w)
        if len(group) > 1:
            components.append([clusters[i] for i in sorted(group)])
    return components


# ----------------------------------------------------------------------
# ST-sets
# ----------------------------------------------------------------------


def is_st_set(s: frozenset, clusters: Iterable[frozenset], ground: frozenset) -> bool:
    """Is ``s`` an ST-set (strict tree set) w.r.t. ``clusters``?

    ``s != ground``, ``s`` is not separated, and ``C|s`` is pairwise compatible.
    """
    if not s or s == ground:
        return False
    clusters = list(clusters)
    if is_separated(s, clusters):
        return False
    inside = list(restrict(clusters, s))
    return all(is_compatible(a, b) for a, b in combinations(inside, 2))


def maximal_st_sets(
    clusters: Iterable[frozenset], blocks: Iterable[frozenset]
) -> list[frozenset]:
    """Partition ``blocks`` into the maximal ST-sets w.r.t. ``clusters``.

    Merging two blocks at a time suffices.  Inside a maximal ST-set ``S`` the
    restriction ``C|S`` is compatible, hence a tree, and any two sibling blocks
    of that tree can be merged without leaving the ST-set family -- so greedy
    merging never stalls before reaching the maximal sets.
    """
    clusters = list(clusters)
    parts = canonical(blocks)
    ground: frozenset = frozenset().union(*parts) if parts else frozenset()

    merged = True
    while merged:
        merged = False
        for i, j in combinations(range(len(parts)), 2):
            union = parts[i] | parts[j]
            if is_st_set(union, clusters, ground):
                parts = canonical(
                    [p for idx, p in enumerate(parts) if idx not in (i, j)] + [union]
                )
                merged = True
                break
    return parts


def collapse(
    clusters: Iterable[frozenset], blocks: Iterable[frozenset]
) -> tuple[list[frozenset], dict[frozenset, list[frozenset]], set[frozenset]]:
    """``Collapse(C)`` of the paper, in block form.

    Returns ``(new_blocks, decollapse_map, new_clusters)``, where
    ``decollapse_map`` sends each new block to the old blocks it swallowed and
    ``new_clusters`` is ``C`` rewritten over the new blocks.  A cluster that
    ends up inside a single block has become trivial and is dropped.
    """
    clusters = list(clusters)
    blocks = list(blocks)
    new_blocks = maximal_st_sets(clusters, blocks)
    ground: frozenset = frozenset().union(*new_blocks) if new_blocks else frozenset()

    decollapse_map = {nb: canonical(b for b in blocks if b <= nb) for nb in new_blocks}
    new_block_set = set(new_blocks)

    new_clusters = {
        r
        for c in clusters
        if (r := round_up(c, new_blocks)) and r != ground and r not in new_block_set
    }
    return new_blocks, decollapse_map, new_clusters


# ----------------------------------------------------------------------
# unseparated sets
# ----------------------------------------------------------------------


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

    Step 1 of the decomposition uses this weaker notion: only the
    "unseparated" half of the ST-set condition is required.

    Unlike ST-sets, a subset of an unseparated set need not itself be
    unseparated, so plain pairwise merging can stall -- with
    ``C = {{a,b},{a,c}}`` on ``{a,b,c,d}`` every pair inside ``{a,b,c}`` is
    separated while ``{a,b,c}`` is not.  Each candidate merge is therefore
    grown to the *smallest* unseparated superset of the pair first, which is
    sound because unseparated sets are closed under union of intersecting
    members.
    """
    clusters = list(clusters)
    parts = canonical(blocks)
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
            parts = canonical([p for p in parts if not p <= grown] + [grown])
            merged = True
            break
    return parts
