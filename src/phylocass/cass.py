"""The Cass algorithm.

Implements van Iersel, Kelk, Rupp and Huson, *Phylogenetic Networks Do not
Need to Be Complex: Using Fewer Reticulations to Represent Conflicting
Clusters*, Bioinformatics 26(12):i124-i131, 2010 (arXiv:0910.3082).

The two halves of the paper map onto the two entry points here:

``cass_simple``
    Algorithm 1 -- build a *simple* level-<=k network for a cluster set whose
    incompatibility graph is connected, by ``k`` rounds of "remove a leaf,
    collapse the maximal ST-sets", then decollapsing and hanging the removed
    leaves back below new reticulations.

``cass``
    The four-step decomposition of Section 3 -- split the clusters by the
    connected components of the incompatibility graph, solve each component
    with ``cass_simple``, and stitch the results into the tree on the
    remaining clusters.

The search itself runs on :class:`~phylocass.workgraph.WorkGraph`; the answer
comes back as a PhyloZoo ``DirectedPhyNetwork``.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Callable, Iterable, Sequence

from phylozoo import DirectedPhyNetwork
from phylozoo.core.network.dnetwork.classifications import (
    level as pz_level,
    reticulation_number as pz_reticulation_number,
)

from .clusters import (
    canonical,
    collapse,
    maximal_unseparated_sets,
    nontrivial_components,
    project,
    remove_taxa,
    restrict,
    round_up,
)
from .io import (
    clusters_of_trees,
    displays_trees,
    per_tree_clusters,
    softwired_clusters,
)
from .treebuild import add_root_edge, build_tree, graft
from .workgraph import DUMMY, WorkGraph, default_namer

__all__ = ["CassOptions", "CassResult", "cass", "cass_from_trees", "cass_simple"]


class CassTimeout(RuntimeError):
    """Raised internally when the search budget is exhausted."""


def _freeze_trees(tree_clusters) -> tuple[frozenset, ...] | None:
    """Normalise per-tree cluster sets into a hashable, memo-key-safe form.

    Empty tree entries are dropped: a tree whose clusters have all become
    trivial at this point in the recursion constrains nothing.
    """
    if tree_clusters is None:
        return None
    frozen = tuple(frozenset(want) for want in tree_clusters)
    return tuple(want for want in frozen if want)


@dataclass
class CassOptions:
    """Tuning knobs for the search."""

    max_level: int = 4
    """Highest level to try before giving up.

    The paper proves Cass is exact for level <= 2; higher levels are a
    heuristic, which is how the published implementation is used in practice.
    """

    time_limit: float | None = None
    """Seconds allowed per conflicting component, or ``None`` for no limit.

    One budget for the whole component, shared across every level the search
    climbs through -- not a fresh allowance per level.
    """

    max_networks: int | None = 20000
    """Cap on the intermediate networks kept per recursive subproblem."""

    display_trees: bool = False
    """Require the network to display the input *trees*, not just their clusters.

    Off (the default), Cass does what the paper describes: every input cluster
    must be represented by some switching, but two clusters of the same tree
    may need different switchings, so the network need not display that tree.

    On, each input tree's clusters must all appear in one and the same
    switching, which means that switching's tree displays the input tree.  The
    network then displays the trees themselves, and its reticulation number is
    an upper bound on their hybridization number -- so this turns Cass into a
    heuristic for Hybridization Number on multiple trees.

    Requires knowing which cluster came from which tree, so it needs
    :func:`cass_from_trees`, or :func:`cass` with ``tree_clusters=``.

    The guarantees in the paper are about the cluster-mode search and do not
    carry over: this mode is a heuristic, and it may need to climb higher than
    ``max_level``.
    """


@dataclass
class CassResult:
    """What Cass produced, plus the numbers worth reporting."""

    network: DirectedPhyNetwork
    level: int
    reticulation_number: int
    clusters: set[frozenset] = field(repr=False, default_factory=set)
    taxa: frozenset = field(repr=False, default_factory=frozenset)
    tree_clusters: list[set[frozenset]] | None = field(repr=False, default=None)
    display_trees: bool = field(repr=False, default=False)

    def represents_input(self) -> bool:
        """Re-check the finished network against the input clusters.

        Goes through PhyloZoo's ``displayed_trees`` rather than the search's
        own cluster routine, so it is an independent verification.
        """
        wanted = {c for c in self.clusters if c}
        if not wanted:
            return True
        return wanted <= softwired_clusters(self.network)

    def displays_input_trees(self) -> bool | None:
        """Does the network display each input tree, not merely its clusters?

        Returns ``None`` when the per-tree provenance is unknown, which is the
        case whenever Cass was handed a bare cluster set.  Like
        :meth:`represents_input`, this verifies through PhyloZoo rather than
        through the search.

        Note that this is worth calling even in cluster mode: the answer is
        often ``True`` there too, just not guaranteed.
        """
        if self.tree_clusters is None:
            return None
        return displays_trees(self.network, self.tree_clusters)

    def to_enewick(self) -> str:
        return self.network.to_string()


# ----------------------------------------------------------------------
# Algorithm 1: simple level-k networks
# ----------------------------------------------------------------------


class _SimpleSearch:
    def __init__(
        self, clusters, blocks, k, options: CassOptions, deadline=None, tree_clusters=None
    ):
        self.k = k
        self.options = options
        if deadline is None and options.time_limit is not None:
            deadline = time.monotonic() + options.time_limit
        self.deadline = deadline
        self.memo: dict[tuple, list[WorkGraph]] = {}
        self.top_clusters = set(clusters)
        self.top_blocks = list(blocks)
        #: per-input-tree cluster sets, or None to only require the clusters
        self.top_trees = _freeze_trees(tree_clusters)

    def check_budget(self) -> None:
        if self.deadline is not None and time.monotonic() > self.deadline:
            raise CassTimeout

    def _accepts(self, net: WorkGraph, clusters, trees) -> bool:
        """The acceptance test, in whichever mode the search is running.

        Cluster mode asks only that every cluster appear somewhere; display
        mode additionally pins each input tree to a single switching.  Display
        mode implies cluster mode, so it is checked on its own.
        """
        if trees is None:
            return net.represents(clusters)
        return net.displays(trees)

    def run(self) -> WorkGraph | None:
        """Return a cleaned simple level-<=k working graph, or ``None``."""
        for candidate in self._recurse(
            self.top_clusters, self.top_trees, self.top_blocks, self.k
        ):
            net = candidate.copy()
            net.clean()
            if not self._accepts(net, self.top_clusters, self.top_trees):
                continue
            if net.level() > self.k:
                continue
            if not net.is_simple():
                continue
            return net
        return None

    def _recurse(self, clusters, trees, blocks, kp):
        """Yield networks on ``blocks`` representing ``clusters`` with ``kp`` reticulations.

        Yields lazily so :meth:`run` can stop at the first usable answer.
        Completed sub-searches are memoised: removing ``x`` then ``y`` reaches
        the same subproblem as removing ``y`` then ``x``.
        """
        key = (frozenset(clusters), trees, frozenset(blocks), kp)
        cached = self.memo.get(key)
        if cached is not None:
            yield from cached
            return

        produced: list[WorkGraph] = []
        for net in self._generate(clusters, trees, blocks, kp):
            produced.append(net)
            yield net
            cap = self.options.max_networks
            if cap is not None and len(produced) >= cap:
                break
        self.memo[key] = produced

    def _generate(self, clusters, trees, blocks, kp):
        self.check_budget()
        # canonical order keeps runs reproducible: several valid networks
        # usually exist, and set iteration order would otherwise pick between
        # them according to Python's hash randomisation
        blocks = canonical(blocks)

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
        for x in blocks + [DUMMY]:
            self.check_budget()
            if not x:
                clusters_removed = set(clusters)
                blocks_removed = list(blocks)
                blocks_sub = blocks_removed
                clusters_sub = clusters_removed
                trees_sub = trees
                decollapse_map = {b: [b] for b in blocks_removed}
            else:
                blocks_removed = [b for b in blocks if b != x]
                if not blocks_removed:
                    continue
                clusters_removed = remove_taxa(clusters, x, ground - x)
                blocks_sub, decollapse_map, clusters_sub = collapse(
                    clusters_removed, blocks_removed
                )
                # each tree's clusters take exactly the same two steps, so a
                # cluster never loses track of the tree it came from
                trees_sub = (
                    None
                    if trees is None
                    else _freeze_trees(
                        project(remove_taxa(want, x, ground - x), blocks_sub)
                        for want in trees
                    )
                )

            for sub in self._recurse(clusters_sub, trees_sub, blocks_sub, kp - 1):
                self.check_budget()
                base = sub.copy()
                if not _decollapse(base, decollapse_map, clusters_removed):
                    continue
                edges = base.edges()
                for i in range(len(edges)):
                    for j in range(i, len(edges)):
                        self.check_budget()
                        net = _hang_below_reticulation(base, edges[i], edges[j], x)
                        if self._accepts(net, clusters, trees):
                            yield net


def _decollapse(net: WorkGraph, decollapse_map, clusters_before) -> bool:
    """Expand every collapsed leaf back into its strict subtree.

    A leaf labelled by a maximal ST-set ``S`` is replaced by the tree on the
    blocks inside ``S`` representing exactly ``C'|S``.  Returns ``False`` if
    such a tree does not exist -- which cannot happen for a genuine ST-set,
    but keeps the search total.
    """
    for v in list(net.leaves()):
        block = net.block.get(v, DUMMY)
        parts = decollapse_map.get(block)
        if not parts or len(parts) == 1:
            continue
        subtree = build_tree(restrict(clusters_before, block), parts)
        if subtree is None:
            return False
        del net.block[v]
        graft(net, v, subtree, {})
    return True


def _hang_below_reticulation(base: WorkGraph, e1, e2, block: frozenset) -> WorkGraph:
    """Add a leaf for ``block`` below a new reticulation fed by ``e1`` and ``e2``."""
    net = base.copy()
    u1, v1 = e1
    if e1 == e2:
        # both reticulation edges come off one edge: subdivide it twice
        w1 = net.subdivide(u1, v1)
        w2 = net.subdivide(w1, v1)
    else:
        u2, v2 = e2
        w1 = net.subdivide(u1, v1)
        w2 = net.subdivide(u2, v2)
    reticulation = net.new_node()
    leaf = net.new_node(block=block)
    net.add_edge(reticulation, leaf)
    net.add_edge(w1, reticulation)
    net.add_edge(w2, reticulation)
    return net


def cass_simple(
    clusters: Iterable[frozenset],
    blocks: Iterable[frozenset],
    k: int,
    options: CassOptions | None = None,
) -> WorkGraph | None:
    """Construct a simple level-<=``k`` network representing ``clusters``.

    ``blocks`` is the taxon set, each taxon a ``frozenset`` of the original
    taxa it stands for.  Returns ``None`` if the search found none.
    """
    net, _ = _cass_simple(clusters, blocks, k, options or CassOptions())
    return net


def _cass_simple(clusters, blocks, k, options, deadline=None, tree_clusters=None):
    """As :func:`cass_simple`, but also reporting whether the budget ran out."""
    search = _SimpleSearch(
        set(clusters), list(blocks), k, options, deadline, tree_clusters
    )
    try:
        return search.run(), False
    except CassTimeout:
        return None, True


# ----------------------------------------------------------------------
# Steps 1-4: the full algorithm
# ----------------------------------------------------------------------


def cass(
    clusters: Iterable[frozenset],
    taxa: Iterable[str] | None = None,
    options: CassOptions | None = None,
    namer: Callable[[frozenset], str] | None = None,
    tree_clusters: Iterable[Iterable[frozenset]] | None = None,
) -> CassResult:
    """Build a rooted phylogenetic network representing all of ``clusters``.

    Follows Steps 1-4 of the paper: decompose along the incompatibility graph,
    solve each non-trivial component with :func:`cass_simple` at the lowest
    level that works, then merge the pieces into the tree on the remaining
    clusters.

    ``tree_clusters`` records which clusters came from which input tree.  It is
    only consulted when ``options.display_trees`` is set, and then the network
    is additionally required to display each of those trees.
    """
    options = options or CassOptions()
    namer = namer or default_namer
    if tree_clusters is not None:
        tree_clusters = [{frozenset(c) for c in want if c} for want in tree_clusters]
    if options.display_trees and tree_clusters is None:
        raise ValueError(
            "display_trees needs to know which cluster came from which tree; "
            "use cass_from_trees(), or pass tree_clusters= explicitly"
        )
    clusters = {frozenset(c) for c in clusters if c}
    if taxa is None:
        taxon_set: frozenset = frozenset().union(*clusters) if clusters else frozenset()
    else:
        taxon_set = frozenset(taxa)
    proper = {c for c in clusters if 1 < len(c) < len(taxon_set)}

    if len(taxon_set) < 2:
        trivial = WorkGraph()
        root = trivial.new_node()
        for t in sorted(taxon_set, key=str):
            trivial.add_edge(root, trivial.new_node(block=frozenset({t})))
        return _finish(
            trivial, 0, 0, proper, taxon_set, namer, tree_clusters, options.display_trees
        )

    # ---- Step 1: decompose along the incompatibility graph ----------------
    component_data = []
    for comp in nontrivial_components(canonical(proper)):
        support: frozenset = frozenset().union(*comp)
        blocks = maximal_unseparated_sets(
            comp, [frozenset({t}) for t in sorted(support, key=str)]
        )
        collapsed = project(comp, blocks)
        # Each component is solved on its own, so the display requirement has
        # to be restricted to it: for tree i, the clusters of this component
        # that tree i contributed.  Switchings in different biconnected
        # components are independent, so satisfying every component with its
        # own switching yields one global switching per tree.  Clusters of
        # tree i that are not in any component sit on the backbone tree, where
        # they hold under every switching.
        component_trees = (
            None
            if not options.display_trees
            else _freeze_trees(
                project(set(comp) & want, blocks) for want in tree_clusters
            )
        )
        component_data.append((comp, support, blocks, collapsed, component_trees))

    # ---- Step 2: a simple level-<=k network per component ----------------
    simple_networks = []
    for comp, support, blocks, collapsed, component_trees in component_data:
        # one shared budget for the whole component, not one per level attempt
        deadline = (
            None if options.time_limit is None else time.monotonic() + options.time_limit
        )
        net = None
        timed_out = False
        reached = 0
        for k in range(1, options.max_level + 1):
            reached = k
            net, timed_out = _cass_simple(
                collapsed, blocks, k, options, deadline, component_trees
            )
            if net is not None or timed_out:
                break
        if net is None:
            reason = (
                f"the {options.time_limit}s budget ran out at level {reached}"
                if timed_out
                else f"no network exists up to level {options.max_level}"
            )
            extra = (
                " Displaying the trees themselves is a stronger requirement than "
                "representing their clusters, so it often needs a higher level."
                if options.display_trees
                else ""
            )
            raise RuntimeError(
                f"Cass gave up on a conflicting component of {len(support)} taxa "
                f"and {len(comp)} clusters: {reason}. Raise CassOptions.max_level "
                f"or CassOptions.time_limit.{extra}"
            )
        simple_networks.append(net)

    # ---- Step 3: the tree on everything that is left ---------------------
    used = {c for comp, *_ in component_data for c in comp}
    star = {c for c in proper if c not in used}
    for _, support, blocks, _, _ in component_data:
        if len(support) < len(taxon_set):
            star.add(support)
        star.update(b for b in blocks if len(b) > 1)

    tree, node_of = build_tree(
        star, [frozenset({t}) for t in sorted(taxon_set, key=str)], return_map=True
    )
    if tree is None:
        raise RuntimeError(
            "internal error: the clusters left after decomposition are incompatible"
        )

    # ---- Step 4: splice each component network into the tree -------------
    for (comp, support, blocks, _, _), net in zip(component_data, simple_networks):
        anchor = node_of[support] if support != taxon_set else tree.root
        for _, child in [(u, v) for u, v in tree.edges() if u == anchor]:
            tree.remove_edge(anchor, child)
        tree.block.pop(anchor, None)
        graft(tree, anchor, net, {b: node_of[b] for b in blocks})

    tree.clean()
    return _finish(
        tree, None, None, proper, taxon_set, namer, tree_clusters, options.display_trees
    )


def _finish(
    work: WorkGraph, level, retics, clusters, taxa, namer, tree_clusters, display_trees
) -> CassResult:
    """Hand the finished graph to PhyloZoo and read the statistics back off it."""
    network = work.to_phylozoo(namer)
    if level is None:
        level = pz_level(network)
    if retics is None:
        retics = pz_reticulation_number(network)
    return CassResult(
        network, level, retics, clusters, taxa, tree_clusters, display_trees
    )


def cass_from_trees(
    trees: Sequence[DirectedPhyNetwork], options: CassOptions | None = None
) -> CassResult:
    """Run Cass on all clusters displayed by a collection of input trees.

    With ``options.display_trees`` set, the network is required to display the
    input trees themselves, making the reticulation number an upper bound on
    their hybridization number.
    """
    trees = list(trees)
    clusters, taxa = clusters_of_trees(trees)
    return cass(
        clusters, taxa, options, tree_clusters=per_tree_clusters(trees, taxa)
    )
