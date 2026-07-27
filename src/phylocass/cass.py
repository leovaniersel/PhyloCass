"""The Cass algorithm.

Implements van Iersel, Kelk, Rupp and Huson, *Phylogenetic Networks Do not
Need to Be Complex: Using Fewer Reticulations to Represent Conflicting
Clusters*, Bioinformatics 26(12):i124-i131, 2010 (arXiv:0910.3082).

The two halves of the paper map onto the two public entry points here:

``cass_simple``
    Algorithm 1 of the paper -- build a *simple* level-<=k network for a
    cluster set whose incompatibility graph is connected, by ``k`` rounds of
    "remove a leaf, collapse the maximal ST-sets", then decollapsing and
    hanging the removed leaves back below new reticulations.

``cass``
    The four-step decomposition of Section 3 -- split the clusters by the
    connected components of the incompatibility graph, solve each component
    with ``cass_simple``, and stitch the results into the tree on the
    remaining clusters.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Iterable, Sequence

from .clusters import (
    collapse,
    maximal_unseparated_sets,
    nontrivial_components,
    remove_taxa,
    restrict,
    round_up,
)
from .network import DUMMY, Network
from .treebuild import add_root_edge, build_tree, graft

__all__ = ["CassOptions", "CassResult", "cass", "cass_simple"]


class CassTimeout(RuntimeError):
    """Raised internally when the search budget is exhausted."""


@dataclass
class CassOptions:
    """Tuning knobs for the search."""

    max_level: int = 4
    """Highest level to try before giving up.

    The paper proves Cass is exact for level <= 2; higher levels are run as a
    heuristic, which is how the published implementation is used in practice.
    """

    time_limit: float | None = None
    """Seconds allowed per connected component, or ``None`` for no limit."""

    max_networks: int | None = 20000
    """Cap on the intermediate networks kept per recursive subproblem."""


@dataclass
class CassResult:
    """What Cass produced, plus the numbers worth reporting."""

    network: Network
    level: int
    reticulation_number: int
    clusters: set[frozenset] = field(repr=False, default_factory=set)
    taxa: frozenset = field(repr=False, default_factory=frozenset)

    def represents_input(self) -> bool:
        """Independent re-check that the output really displays every input cluster."""
        return self.network.represents(self.clusters)

    def to_enewick(self, namer=None) -> str:
        return self.network.to_enewick(namer)


# ----------------------------------------------------------------------
# Algorithm 1: simple level-k networks
# ----------------------------------------------------------------------


class _SimpleSearch:
    def __init__(self, clusters, blocks, k, options: CassOptions):
        self.k = k
        self.options = options
        self.deadline = (
            None if options.time_limit is None else time.monotonic() + options.time_limit
        )
        self.memo: dict[tuple, list[Network]] = {}
        self.top_clusters = set(clusters)
        self.top_blocks = list(blocks)

    def check_budget(self) -> None:
        if self.deadline is not None and time.monotonic() > self.deadline:
            raise CassTimeout

    def run(self) -> Network | None:
        """Return a cleaned simple level-<=k network, or ``None``."""
        for net in self._recurse(self.top_clusters, self.top_blocks, self.k):
            candidate = net.copy()
            candidate.clean()
            if not candidate.represents(self.top_clusters):
                continue
            if candidate.level() > self.k:
                continue
            if not candidate.is_simple():
                continue
            return candidate
        return None

    def _recurse(self, clusters, blocks, kp):
        """Yield networks on ``blocks`` representing ``clusters`` with ``kp`` reticulations.

        Yields lazily so that :meth:`run` can stop at the first usable answer.
        The results of complete sub-searches are memoised, since removing
        ``x`` then ``y`` reaches the same subproblem as removing ``y`` then
        ``x``.
        """
        key = (frozenset(clusters), frozenset(blocks), kp)
        cached = self.memo.get(key)
        if cached is not None:
            yield from cached
            return

        produced: list[Network] = []
        for net in self._generate(clusters, blocks, kp):
            produced.append(net)
            yield net
            if self.options.max_networks is not None and len(produced) >= self.options.max_networks:
                break
        self.memo[key] = produced

    def _generate(self, clusters, blocks, kp):
        self.check_budget()
        blocks = list(blocks)

        if kp == 0:
            tree = build_tree(clusters, blocks)
            if tree is not None:
                add_root_edge(tree)
                yield tree
            return

        ground: frozenset = frozenset().union(*blocks) if blocks else frozenset()

        # x ranges over the current taxa plus a dummy taxon d.  Removing d
        # removes nothing and -- per the paper -- deliberately skips the
        # collapse, which is what lets Cass reach reticulations of indegree 3.
        for x in list(blocks) + [DUMMY]:
            self.check_budget()
            if x is DUMMY or not x:
                clusters_removed = set(clusters)
                blocks_removed = list(blocks)
                blocks_sub = blocks_removed
                decollapse_map = {b: [b] for b in blocks_removed}
                clusters_sub = clusters_removed
            else:
                blocks_removed = [b for b in blocks if b != x]
                if not blocks_removed:
                    continue
                clusters_removed = remove_taxa(clusters, x, ground - x)
                blocks_sub, decollapse_map, clusters_sub = collapse(
                    clusters_removed, blocks_removed
                )

            for sub in self._recurse(clusters_sub, blocks_sub, kp - 1):
                self.check_budget()
                base = sub.copy()
                if not _decollapse(base, decollapse_map, clusters_removed):
                    continue
                edges = base.edges()
                for i in range(len(edges)):
                    for j in range(i, len(edges)):
                        self.check_budget()
                        net = _hang_below_reticulation(base, edges[i], edges[j], x)
                        if net is not None and net.represents(clusters):
                            yield net


def _decollapse(net: Network, decollapse_map, clusters_before) -> bool:
    """Expand every collapsed leaf back into its strict subtree.

    A leaf labelled by a maximal ST-set ``S`` is replaced by the tree on the
    blocks inside ``S`` representing exactly ``C'|S``.  Returns ``False`` if
    some such tree does not exist (which cannot happen for a genuine ST-set,
    but keeps the search total).
    """
    for v in list(net.leaves()):
        label = net.label.get(v, DUMMY)
        parts = decollapse_map.get(label)
        if not parts or len(parts) == 1:
            continue
        subtree = build_tree(restrict(clusters_before, label), parts)
        if subtree is None:
            return False
        del net.label[v]
        graft(net, v, subtree, {})
    return True


def _hang_below_reticulation(base: Network, e1, e2, label: frozenset) -> Network | None:
    """Add a leaf labelled ``label`` below a new reticulation fed by ``e1`` and ``e2``."""
    net = base.copy()
    u1, v1 = e1
    if e1 == e2:
        # both reticulation edges come off the same edge: subdivide it twice
        w1 = net.subdivide(u1, v1)
        w2 = net.subdivide(w1, v1)
    else:
        u2, v2 = e2
        w1 = net.subdivide(u1, v1)
        w2 = net.subdivide(u2, v2)
    ret = net.new_node()
    leaf = net.new_node(label=label)
    net.add_edge(ret, leaf)
    net.add_edge(w1, ret)
    net.add_edge(w2, ret)
    return net


def cass_simple(
    clusters: Iterable[frozenset],
    blocks: Iterable[frozenset],
    k: int,
    options: CassOptions | None = None,
) -> Network | None:
    """Construct a simple level-<=``k`` network representing ``clusters``.

    ``blocks`` is the taxon set, each taxon given as a ``frozenset`` of the
    original taxa it stands for.  Returns ``None`` if the search found none.
    """
    options = options or CassOptions()
    search = _SimpleSearch(set(clusters), list(blocks), k, options)
    try:
        return search.run()
    except CassTimeout:
        return None


# ----------------------------------------------------------------------
# Steps 1-4: the full algorithm
# ----------------------------------------------------------------------


def cass(
    clusters: Iterable[frozenset],
    taxa: Iterable[str] | None = None,
    options: CassOptions | None = None,
) -> CassResult:
    """Build a rooted phylogenetic network representing all of ``clusters``.

    Follows Steps 1-4 of the paper: decompose along the incompatibility graph,
    solve each non-trivial component with :func:`cass_simple` at the lowest
    level that works, then merge the pieces into the tree on the remaining
    clusters.
    """
    options = options or CassOptions()
    clusters = {frozenset(c) for c in clusters if c}
    if taxa is None:
        taxon_set: frozenset = frozenset().union(*clusters) if clusters else frozenset()
    else:
        taxon_set = frozenset(taxa)
    proper = {c for c in clusters if 1 < len(c) < len(taxon_set)}

    if len(taxon_set) < 2:
        net = Network()
        root = net.new_node()
        for t in sorted(taxon_set, key=str):
            net.add_edge(root, net.new_node(label=frozenset({t})))
        return CassResult(net, 0, 0, proper, taxon_set)

    # ---- Step 1: decompose along the incompatibility graph ----------------
    components = nontrivial_components(sorted(proper, key=lambda c: (len(c), sorted(map(str, c)))))

    component_data = []
    for comp in components:
        support: frozenset = frozenset().union(*comp)
        blocks = maximal_unseparated_sets(comp, [frozenset({t}) for t in support])
        collapsed = set()
        for c in comp:
            r = round_up(c, blocks)
            if r and r != support and r not in blocks:
                collapsed.add(r)
        component_data.append((comp, support, blocks, collapsed))

    # ---- Step 2: a simple level-<=k network per component ----------------
    simple_networks = []
    for comp, support, blocks, collapsed in component_data:
        net = None
        for k in range(1, options.max_level + 1):
            net = cass_simple(collapsed, blocks, k, options)
            if net is not None:
                break
        if net is None:
            raise RuntimeError(
                "Cass found no network up to level "
                f"{options.max_level} for a component on {len(support)} taxa; "
                "raise CassOptions.max_level or the time limit"
            )
        simple_networks.append(net)

    # ---- Step 3: the tree on everything that is left ---------------------
    used = set()
    for comp, _, _, _ in component_data:
        used.update(comp)
    star = {c for c in proper if c not in used}
    for _, support, blocks, _ in component_data:
        if len(support) < len(taxon_set):
            star.add(support)
        star.update(b for b in blocks if len(b) > 1)

    tree, node_of = build_tree(star, [frozenset({t}) for t in taxon_set], return_map=True)
    if tree is None:
        raise RuntimeError(
            "internal error: the clusters left after decomposition are incompatible"
        )

    # ---- Step 4: splice each component network into the tree -------------
    for (comp, support, blocks, _), net in zip(component_data, simple_networks):
        anchor = node_of[support] if support != taxon_set else tree.root
        for child in list(tree.children[anchor]):
            tree.remove_edge(anchor, child)
        identify = {b: node_of[b] for b in blocks}
        if anchor in tree.label:
            del tree.label[anchor]
        graft(tree, anchor, net, identify)

    tree.clean()
    return CassResult(tree, tree.level(), tree.reticulation_number(), proper, taxon_set)


def cass_from_trees(trees: Sequence[Network], options: CassOptions | None = None) -> CassResult:
    """Run Cass on all clusters displayed by a collection of input trees."""
    from .clusters import clusters_of_trees

    clusters, taxa = clusters_of_trees(trees)
    return cass(clusters, taxa, options)
