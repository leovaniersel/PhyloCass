"""The working representation used while Cass searches.

PhyloZoo's :class:`~phylozoo.DirectedPhyNetwork` is immutable and validates
node degrees: internal nodes must be tree nodes (in-degree 1, out-degree >= 2)
or hybrid nodes (in-degree >= 2, out-degree 1).  Every intermediate state of
the Cass construction violates that -- subdividing an edge creates a degree-2
node, and the algorithm deliberately parks dummy leaves and a dummy root that
only make sense mid-search.

So the search runs on PhyloZoo's mutable ``DirectedMultiGraph`` primitive
instead, wrapped here as :class:`WorkGraph`, and a validated
``DirectedPhyNetwork`` is materialised only once a candidate has been cleaned
up.  PhyloZoo still supplies the graph engine and the structural analysis
(biconnected components, cut-edges, degree-2 suppression, parallel-edge
identification); what :class:`WorkGraph` adds is the phylogenetic bookkeeping
Cass needs on top of it.

Every leaf carries a *block*: the ``frozenset`` of original taxa it currently
stands for.  The empty block marks a dummy leaf, which therefore contributes
nothing to any descendant set for free.
"""

from __future__ import annotations

from itertools import product
from typing import Callable, Iterable, Iterator

from phylozoo import DirectedPhyNetwork
from phylozoo.core.primitives.d_multigraph.base import DirectedMultiGraph
from phylozoo.core.primitives.d_multigraph.features import (
    biconnected_components,
    cut_edges,
)
from phylozoo.core.primitives.d_multigraph.transformations import (
    identify_parallel_edge,
    suppress_degree2_node,
)

DUMMY: frozenset = frozenset()

__all__ = ["DUMMY", "WorkGraph", "default_namer"]


def default_namer(block: frozenset) -> str:
    """Name a leaf after the taxa in its block."""
    return "+".join(sorted(str(t) for t in block))


class WorkGraph:
    """A Cass-intermediate network: a PhyloZoo multigraph plus leaf blocks."""

    __slots__ = ("g", "block", "_next")

    def __init__(self) -> None:
        self.g = DirectedMultiGraph()
        self.block: dict[int, frozenset] = {}
        self._next = 0

    # ------------------------------------------------------------------
    # construction
    # ------------------------------------------------------------------
    def new_node(self, block: frozenset | None = None) -> int:
        v = self._next
        self._next += 1
        self.g.add_node(v)
        if block is not None:
            self.block[v] = block
        return v

    def add_edge(self, u: int, v: int) -> None:
        self.g.add_edge(u, v)

    def remove_edge(self, u: int, v: int, key: int | None = None) -> None:
        self.g.remove_edge(u, v, key)

    def remove_node(self, v: int) -> None:
        self.g.remove_node(v)
        self.block.pop(v, None)

    def subdivide(self, u: int, v: int) -> int:
        """Insert a fresh node in the middle of edge ``(u, v)`` and return it."""
        self.remove_edge(u, v)
        w = self.new_node()
        self.add_edge(u, w)
        self.add_edge(w, v)
        return w

    @classmethod
    def from_phylozoo(cls, network: DirectedPhyNetwork) -> "WorkGraph":
        """Wrap a PhyloZoo network so the blob-wise machinery applies to it.

        Each taxon becomes a one-element block, which is what the rest of this
        class expects.
        """
        work = cls()
        mapping = {}
        for v in network.nodes:
            label = network.get_label(v)
            is_leaf = not any(True for _ in network.children(v))
            mapping[v] = work.new_node(
                block=frozenset({label}) if (is_leaf and label is not None) else None
            )
        for u, v in network.edges:
            work.add_edge(mapping[u], mapping[v])
        return work

    def copy(self) -> "WorkGraph":
        other = WorkGraph.__new__(WorkGraph)
        other.g = self.g.copy()
        other.block = dict(self.block)
        other._next = self._next
        return other

    # ------------------------------------------------------------------
    # basic queries
    # ------------------------------------------------------------------
    def nodes(self) -> list[int]:
        return list(self.g.nodes)

    def edges(self) -> list[tuple[int, int]]:
        return list(self.g.edges)

    def adjacency(self) -> tuple[dict[int, list[int]], dict[int, list[int]]]:
        """Children and parents as plain dicts.

        Snapshotting once is much faster than repeatedly querying the graph,
        and it counts parallel edges correctly.
        """
        children: dict[int, list[int]] = {v: [] for v in self.g.nodes}
        parents: dict[int, list[int]] = {v: [] for v in self.g.nodes}
        for u, v in self.g.edges:
            children[u].append(v)
            parents[v].append(u)
        return children, parents

    def roots(self) -> list[int]:
        return [v for v in self.g.nodes if self.g.indegree(v) == 0]

    @property
    def root(self) -> int:
        rs = self.roots()
        if len(rs) != 1:
            raise ValueError(f"working graph has {len(rs)} roots, expected exactly 1")
        return rs[0]

    def leaves(self) -> list[int]:
        return [v for v in self.g.nodes if self.g.outdegree(v) == 0]

    def taxa(self) -> frozenset:
        out: frozenset = frozenset()
        for v in self.leaves():
            out |= self.block.get(v, DUMMY)
        return out

    def blocks(self) -> list[frozenset]:
        """Leaf blocks, dummies excluded."""
        return [b for v in self.leaves() if (b := self.block.get(v, DUMMY))]

    def reticulation_number(self) -> int:
        return sum(max(0, self.g.indegree(v) - 1) for v in self.g.nodes)

    def topological_order(self, children=None, parents=None) -> list[int]:
        if children is None or parents is None:
            children, parents = self.adjacency()
        remaining = {v: len(parents[v]) for v in children}
        queue = [v for v, d in remaining.items() if d == 0]
        order: list[int] = []
        while queue:
            v = queue.pop()
            order.append(v)
            for w in children[v]:
                remaining[w] -= 1
                if remaining[w] == 0:
                    queue.append(w)
        if len(order) != len(children):
            raise ValueError("working graph contains a directed cycle")
        return order

    # ------------------------------------------------------------------
    # softwired clusters
    # ------------------------------------------------------------------
    def switching_cluster_sets(self) -> list[set[frozenset]]:
        """The clusters of each switching separately, one set per switching.

        A *switching* keeps one incoming edge per reticulation; what is left is
        a tree, and its clusters are the descendant taxon sets of its nodes.

        Keeping the switchings apart is what distinguishes the two acceptance
        tests: :meth:`represents` asks whether every cluster appears in *some*
        switching, :meth:`displays` whether each tree's clusters all appear in
        *one and the same* switching.

        This is computed here rather than through PhyloZoo's ``displayed_trees``
        because mid-search graphs are not valid ``DirectedPhyNetwork``s, and
        because this is the hot loop of the whole algorithm.  The test suite
        cross-checks it against ``displayed_trees`` on finished networks.
        """
        children, parents = self.adjacency()
        try:
            order = self.topological_order(children, parents)
        except ValueError:
            return []
        roots = [v for v in children if not parents[v]]
        if len(roots) != 1:
            return []
        root = roots[0]

        reticulations = [v for v in order if len(parents[v]) >= 2]
        per_switching: list[set[frozenset]] = []

        for choice in product(*(parents[r] for r in reticulations)):
            on = dict(zip(reticulations, choice))

            reachable = {root}
            for v in order:
                if v not in reachable:
                    continue
                for w in children[v]:
                    if len(parents[w]) >= 2 and on.get(w) != v:
                        continue
                    reachable.add(w)

            found: set[frozenset] = set()
            below: dict[int, frozenset] = {}
            for v in reversed(order):
                if v not in reachable:
                    continue
                acc = self.block.get(v, DUMMY)
                for w in children[v]:
                    if w in reachable and not (
                        len(parents[w]) >= 2 and on.get(w) != v
                    ):
                        acc = acc | below[w]
                below[v] = acc
                if acc:
                    found.add(acc)
            per_switching.append(found)
        return per_switching

    def blobs(self) -> list[set[int]]:
        """Node sets of the biconnected components that contain a reticulation."""
        edges = self.edges()
        out = []
        for comp in biconnected_components(self.g):
            inside = sum(1 for u, v in edges if u in comp and v in comp)
            if inside - len(comp) + 1 > 0:
                out.append(set(comp))
        return out

    def softwired_clusters(self) -> set[frozenset]:
        """Every cluster this network represents in the softwired sense.

        A cluster is represented if it appears in some switching; different
        clusters may use different switchings.

        Computed blob by blob rather than by enumerating all switchings at
        once. Switchings in different biconnected components are independent,
        and -- the fact that makes this work -- the taxa below the head of a
        cut-edge are the *same* under every switching: both parents of any
        reticulation down there lie in that same blob, so no switching can cut
        it off. So a node's descendant set varies only with the switching of
        the one blob containing it, and the union over all switchings is the
        union over each blob's switchings taken separately.

        That turns ``2 ** (total reticulations)`` into a sum of
        ``2 ** (reticulations in one blob)``, which is what makes a network
        with many small conflicts tractable.
        """
        fixed, per_blob = self.blob_cluster_sets()
        found = set(fixed)
        for sets in per_blob:
            for one in sets:
                found |= one
        return found

    def blob_cluster_sets(self) -> tuple[set[frozenset], list[list[set[frozenset]]]]:
        """Clusters split into the switching-invariant ones and per-blob ones.

        Returns ``(fixed, per_blob)`` where ``fixed`` holds the clusters present
        under every switching, and ``per_blob[i]`` lists one cluster set per
        switching of blob ``i``.

        This is the decomposition both :meth:`softwired_clusters` and
        :meth:`displays` are built on.
        """
        children, parents = self.adjacency()
        try:
            order = self.topological_order(children, parents)
        except ValueError:
            return set(), []
        if len([v for v in children if not parents[v]]) != 1:
            return set(), []

        # descendant taxa ignoring switchings; correct as-is for every node
        # outside a blob, and the value to use for anything hanging below one
        reach: dict[int, frozenset] = {}
        for v in reversed(order):
            acc = self.block.get(v, DUMMY)
            for w in children[v]:
                acc = acc | reach[w]
            reach[v] = acc

        blob_list = self.blobs()
        in_blob: dict[int, int] = {}
        for index, nodes in enumerate(blob_list):
            for v in nodes:
                in_blob[v] = index

        # A node's descendant set is switching-invariant exactly when it heads a
        # cut-edge (or is the root): no switching can detach anything below it.
        # That includes the *top of each blob*, which the enclosing blob can
        # also produce -- by switching all its other children off -- so leaving
        # those out would let one cluster belong to two blobs and break the
        # independence :meth:`displays` relies on.
        invariant = {v for _, v in cut_edges(self.g)}
        invariant.add(order[0])
        fixed = {reach[v] for v in invariant if reach[v]}
        per_blob = [
            self._blob_clusters(nodes, index, in_blob, children, order, reach)
            for index, nodes in enumerate(blob_list)
        ]
        return fixed, per_blob

    def _blob_clusters(
        self, nodes, index, in_blob, children, order, reach
    ) -> list[set[frozenset]]:
        """One cluster set per switching of this blob."""
        inside_children = {
            v: [w for w in children[v] if in_blob.get(w) == index] for v in nodes
        }
        # everything hanging off the blob is switching-invariant
        pendant = {
            v: frozenset().union(
                self.block.get(v, DUMMY),
                *(reach[w] for w in children[v] if in_blob.get(w) != index),
            )
            for v in nodes
        }

        inside_parents: dict[int, list[int]] = {v: [] for v in nodes}
        for v in nodes:
            for w in inside_children[v]:
                inside_parents[w].append(v)
        reticulations = [v for v in nodes if len(inside_parents[v]) >= 2]
        local_order = [v for v in order if v in nodes]

        per_switching: list[set[frozenset]] = []
        for choice in product(*(inside_parents[r] for r in reticulations)):
            on = dict(zip(reticulations, choice))
            found: set[frozenset] = set()
            below: dict[int, frozenset] = {}
            for v in reversed(local_order):
                acc = pendant[v]
                for w in inside_children[v]:
                    if len(inside_parents[w]) >= 2 and on.get(w) != v:
                        continue
                    acc = acc | below[w]
                below[v] = acc
                if acc:
                    found.add(acc)
            per_switching.append(found)
        return per_switching

    def represents(self, clusters: Iterable[frozenset]) -> bool:
        """Is every cluster represented, each free to use its own switching?"""
        wanted = {c for c in clusters if c}
        if not wanted:
            return True
        return wanted <= self.softwired_clusters()

    def displays(self, tree_clusters: Iterable[Iterable[frozenset]]) -> bool:
        """Does one switching account for the whole of each input tree?

        ``tree_clusters`` gives the clusters of each input tree separately.
        For every tree there must be a *single* switching containing all of its
        clusters -- then that switching's tree displays the input tree, so the
        network displays the trees themselves and not merely their clusters.

        Different trees may of course use different switchings; that is exactly
        what the reticulations are for.

        Decided blob by blob, like :meth:`softwired_clusters`, which needs one
        extra fact: a cluster that is *not* switching-invariant is produced by
        exactly one blob. A node above a blob reaches either all of that blob's
        taxa or none of them -- the taxa below a cut-edge are fixed -- so it can
        never produce a proper subset of them, and blobs that are side by side
        span disjoint taxa. So each such cluster can be attributed to one blob,
        and the blobs can then be satisfied independently.
        """
        wanted_per_tree = [{c for c in want if c} for want in tree_clusters]
        wanted_per_tree = [w for w in wanted_per_tree if w]
        if not wanted_per_tree:
            return True

        fixed, per_blob = self.blob_cluster_sets()
        reachable = [set().union(*sets) if sets else set() for sets in per_blob]
        anywhere = set().union(*reachable) if reachable else set()

        for want in wanted_per_tree:
            outstanding = want - fixed
            if not outstanding:
                continue
            if not outstanding <= anywhere:
                return False
            for sets, whole in zip(per_blob, reachable):
                needed = outstanding & whole
                if needed and not any(needed <= one for one in sets):
                    return False
        return True

    # ------------------------------------------------------------------
    # structure, via PhyloZoo
    # ------------------------------------------------------------------
    def level(self) -> int:
        """Largest reticulation number over the biconnected components."""
        edges = self.edges()
        best = 0
        for comp in biconnected_components(self.g):
            # an edge belongs to the component containing both its endpoints;
            # two components share at most one node, so this is unambiguous
            inside = sum(1 for u, v in edges if u in comp and v in comp)
            best = max(best, inside - len(comp) + 1)
        return best

    def is_simple(self) -> bool:
        """Is the head of every cut-edge a leaf?"""
        children, _ = self.adjacency()
        return all(not children[v] for _, v in cut_edges(self.g))

    # ------------------------------------------------------------------
    # clean-up
    # ------------------------------------------------------------------
    def _drop_parallel_edges(self) -> bool:
        seen: set[tuple[int, int]] = set()
        duplicated: set[tuple[int, int]] = set()
        for u, v in self.g.edges:
            if (u, v) in seen:
                duplicated.add((u, v))
            seen.add((u, v))
        for u, v in duplicated:
            identify_parallel_edge(self.g, u, v)
        return bool(duplicated)

    def tidy(self) -> None:
        """The clean-up rules from the paper, applied to a fixpoint.

        Delete unlabelled dead ends, identify parallel edges, suppress
        degree-2 nodes, and drop an out-degree-1 root.
        """
        changed = True
        while changed:
            changed = False

            for v in self.nodes():
                if self.g.outdegree(v) == 0 and v not in self.block:
                    self.remove_node(v)
                    changed = True

            if self._drop_parallel_edges():
                changed = True

            for v in self.nodes():
                if v in self.block:
                    continue
                if self.g.indegree(v) == 1 and self.g.outdegree(v) == 1:
                    suppress_degree2_node(self.g, v)
                    changed = True

            for v in self.roots():
                if v not in self.block and self.g.outdegree(v) == 1:
                    self.remove_node(v)
                    changed = True
                    break

    def contract_reticulation_edges(self) -> None:
        """Contract every edge whose tail and head are both reticulations.

        This is what turns the dummy-taxon detour into a genuine indegree-3
        reticulation: hang ``d``, hang the real taxon, delete ``d``, and the
        edge left between the two reticulations collapses.
        """
        changed = True
        while changed:
            changed = False
            for u, v in self.edges():
                if self.g.indegree(u) >= 2 and self.g.indegree(v) >= 2:
                    ups = [p for p, _ in self._parent_edges(u)]
                    downs = [w for w in self._child_targets(u) if w != v]
                    self.remove_node(u)
                    for p in ups:
                        self.add_edge(p, v)
                    for w in downs:
                        self.add_edge(v, w)
                    changed = True
                    break

    def _parent_edges(self, v: int) -> list[tuple[int, int]]:
        return [(u, w) for u, w in self.g.edges if w == v]

    def _child_targets(self, v: int) -> list[int]:
        return [w for u, w in self.g.edges if u == v]

    def clean(self) -> None:
        """Full output clean-up: strip dummies, tidy, contract, tidy again."""
        for v in [v for v in self.leaves() if self.block.get(v, DUMMY) == DUMMY]:
            self.remove_node(v)
        self.tidy()
        self.contract_reticulation_edges()
        self.tidy()

    # ------------------------------------------------------------------
    # hand-off to PhyloZoo
    # ------------------------------------------------------------------
    def to_phylozoo(
        self, namer: Callable[[frozenset], str] | None = None
    ) -> DirectedPhyNetwork:
        """Materialise a validated :class:`~phylozoo.DirectedPhyNetwork`.

        Raises whatever PhyloZoo raises if the graph is not a legal
        phylogenetic network -- which is the point: it is the independent
        check that Cass produced something well formed.
        """
        namer = namer or default_namer
        nodes: list = []
        for v in self.nodes():
            block = self.block.get(v)
            if block:
                nodes.append((v, {"label": namer(block)}))
            else:
                nodes.append(v)
        return DirectedPhyNetwork(edges=self.edges(), nodes=nodes)

    def __repr__(self) -> str:  # pragma: no cover - debugging aid
        return (
            f"<WorkGraph nodes={self.g.number_of_nodes()} "
            f"edges={self.g.number_of_edges()} retic={self.reticulation_number()}>"
        )
