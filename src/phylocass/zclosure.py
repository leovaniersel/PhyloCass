"""Z-closure, for input trees whose taxon sets differ but overlap.

Cass needs clusters of one taxon set.  A tree on a subset of the taxa only
gives *partial* clusters: ``{a,b}`` seen on ``{a,b,c}`` says a and b group
together against c, and says nothing at all about the taxa that tree is
missing.  Feeding such a cluster to Cass as though it were a cluster of the
full taxon set asserts more than the tree does.

The Z-closure of Huson, Dezulian, Kloepper and Steel (2004) fills the gaps in:
partial splits are combined by a rule until no new ones appear, and whichever
have grown to cover the whole taxon set are then genuine full clusters.

The rule, as restated by Gruenewald, Huber and Moulton, is for two partial
splits ``S1 = A1|A1~`` and ``S2 = A2|A2~``::

    if  A1 n A2 != {},  A2 n A1~ != {},  A1~ n A2~ != {},  and  A1 n A2~ == {}
    then produce   A1 | (A1~ u A2~)   and   (A1 u A2) | A2~

Exactly one of the four intersections is empty, which is where the name comes
from.

Rooted clusters need one adjustment.  A split does not say which side is which,
but a cluster does, and all the input trees share a root -- so conceptually
every tree gets an extra root taxon that always sits opposite the cluster.
That pins the orientation: the cluster is the side without the root taxon.

The rule lets ``A_i`` be *either* part of split ``i``, so with ``C`` the cluster
and ``R = X \\ C`` its root side (which the root taxon joins), there are four
combinations.  One of them, taking the root side of both splits and asking
``A1 n A2~ == {}``, can never fire, because the root taxon lies in both root
sides.  Another is the first with the two clusters swapped.  What is left, for
a cluster ``C`` known on taxon set ``X``, are two rules:

**Overlapping.**  If ``C1 n C2 != {}``, ``C2 n R1 != {}`` and ``C1 n R2 == {}``
-- so everything C1 knows about X2 lies inside C2 -- then::

    C1 is known on X1 u R2      and      C1 u C2 is a cluster known on C1 u X2

**Disjoint.**  If ``C1 n C2 == {}``, ``C1 n R2 != {}`` and ``C2 n R1 != {}`` --
two clades that are apart, each with a witness outside the other -- then::

    C1 is known on X1 u C2

(the matching ``C2 is known on C1 u X2`` is the same rule on the swapped pair).

Widening the taxon set is what does the real work: it is how a partial cluster
eventually becomes full.  The disjoint rule matters more than it looks -- most
clusters of an input tree are disjoint from most clusters of another, and
without it very little reaches the full taxon set.

This is a heuristic and an incomplete one -- Z-closure does not derive every
cluster implied by the input, which is what motivated the later M- and Y-rules.
Clusters that never grow to the full taxon set are dropped.

Reference: D. Huson, T. Dezulian, T. Kloepper, M. Steel, "Phylogenetic
super-networks from partial trees", IEEE/ACM TCBB 1(4):151-158, 2004.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Sequence

from .clusters import canonical

__all__ = ["PartialCluster", "ZClosureResult", "partial_clusters", "z_closure"]


@dataclass(frozen=True)
class PartialCluster:
    """A cluster together with the taxon set it is actually known on."""

    cluster: frozenset
    known_on: frozenset

    def __post_init__(self) -> None:
        if not self.cluster or not (self.cluster < self.known_on):
            raise ValueError("a partial cluster must be a proper non-empty subset")


@dataclass
class ZClosureResult:
    """What the closure produced."""

    clusters: set[frozenset]
    """Clusters that grew to cover the whole taxon set."""

    taxa: frozenset
    tree_clusters: list[set[frozenset]] = field(default_factory=list)
    """Per input tree, its own clusters that became full (provenance for display mode)."""

    partial_total: int = 0
    """How many distinct clusters the closure knew about when it stopped."""

    dropped: int = 0
    """How many never reached the full taxon set and were discarded."""

    rounds: int = 0
    hit_limit: bool = False


def partial_clusters(
    trees: Sequence["object"],
) -> tuple[list[tuple[frozenset, frozenset, int]], frozenset]:
    """Every non-trivial cluster of every tree, tagged with its taxon set and tree.

    Trees are expected to expose ``.taxa``; clusters come from
    :func:`phylocass.io.hardwired_clusters`.  A cluster equal to a tree's whole
    taxon set is dropped: the root of a partial tree says nothing about whether
    its taxa form a clade of the full tree.
    """
    from .io import hardwired_clusters

    out: list[tuple[frozenset, frozenset, int]] = []
    taxa: frozenset = frozenset()
    for index, tree in enumerate(trees):
        tree_taxa = frozenset(tree.taxa)
        taxa |= tree_taxa
        for cluster in hardwired_clusters(tree):
            if 1 < len(cluster) < len(tree_taxa):
                out.append((cluster, tree_taxa, index))
    return out, taxa


def z_closure(
    partials: Iterable[tuple[frozenset, frozenset, int]],
    taxa: frozenset,
    max_clusters: int = 20000,
    max_rounds: int = 64,
) -> ZClosureResult:
    """Close a set of partial clusters under the (rooted) Z-rule.

    ``partials`` is a sequence of ``(cluster, taxon set it is known on, tree
    index)``.  Returns the clusters that ended up known on all of ``taxa``.

    Only the widest taxon set per cluster is kept: if a cluster is known on
    ``X`` and on ``X'`` it is known on ``X u X'``, so the state is a map from
    cluster to taxon set that only ever grows.  That makes the closure a
    genuine fixed point -- independent of the order rules are applied in, and
    therefore deterministic -- unlike the "replace the two inputs" reading of
    the rule, which is order-dependent.
    """
    known: dict[frozenset, frozenset] = {}
    origin: dict[frozenset, set[int]] = {}
    n_trees = 0

    for cluster, known_on, index in partials:
        known[cluster] = known.get(cluster, frozenset()) | known_on
        origin.setdefault(cluster, set()).add(index)
        n_trees = max(n_trees, index + 1)

    hit_limit = False
    rounds = 0
    for rounds in range(1, max_rounds + 1):
        changed = False
        items = canonical(known)
        for c1 in items:
            x1 = known[c1]
            r1 = x1 - c1
            for c2 in items:
                if c1 == c2:
                    continue
                x2 = known[c2]
                r2 = x2 - c2

                if c1 & c2:
                    # overlapping: C2 n R1 != {} and C1 n R2 == {}
                    if not (c2 & r1) or (c1 & r2):
                        continue
                    widened = x1 | r2
                    merged, merged_on = c1 | c2, c1 | x2
                else:
                    # disjoint: C1 n R2 != {} and C2 n R1 != {}
                    if not (c1 & r2) or not (c2 & r1):
                        continue
                    widened = x1 | c2
                    merged = merged_on = None

                if widened != known[c1]:
                    known[c1] = widened
                    x1, r1 = widened, widened - c1
                    changed = True

                if merged is not None and 1 < len(merged) < len(taxa):
                    before = known.get(merged)
                    if before is None or not (merged_on <= before):
                        if before is None and len(known) >= max_clusters:
                            hit_limit = True
                        else:
                            known[merged] = (before or frozenset()) | merged_on
                            changed = True
        if not changed or hit_limit:
            break

    full = {c for c, x in known.items() if x == taxa and 1 < len(c) < len(taxa)}

    # provenance: a full cluster belongs to tree i if tree i contributed it.
    # Clusters the rule merged into existence belong to no single tree, so
    # display mode does not pin them to one switching.
    tree_clusters: list[set[frozenset]] = [set() for _ in range(n_trees)]
    for cluster in full:
        for index in origin.get(cluster, ()):
            tree_clusters[index].add(cluster)

    return ZClosureResult(
        clusters=full,
        taxa=taxa,
        tree_clusters=tree_clusters,
        partial_total=len(known),
        dropped=len(known) - len(full),
        rounds=rounds,
        hit_limit=hit_limit,
    )
