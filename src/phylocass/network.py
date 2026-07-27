"""Rooted phylogenetic networks: structure, softwired clusters, level.

A network is a rooted DAG whose leaves carry labels.  A label is a ``frozenset``
of original taxa -- a *block* in the sense of :mod:`phylocass.clusters`.  The
empty frozenset marks a *dummy* leaf: the placeholder taxon ``d`` that Cass
temporarily hangs below a reticulation and strips again before output.  Because
a dummy contributes no taxa, it drops out of every descendant set for free.
"""

from __future__ import annotations

from itertools import product
from typing import Iterable, Iterator

DUMMY: frozenset = frozenset()

__all__ = ["Network", "DUMMY"]


class Network:
    """A rooted DAG with labelled leaves.

    Nodes are integers.  Parallel edges are not stored twice; the tidy-up step
    collapses them, matching the convention of the paper.
    """

    __slots__ = ("_next_id", "children", "parents", "label")

    def __init__(self) -> None:
        self._next_id = 0
        self.children: dict[int, list[int]] = {}
        self.parents: dict[int, list[int]] = {}
        self.label: dict[int, frozenset] = {}

    # ------------------------------------------------------------------
    # construction
    # ------------------------------------------------------------------
    def new_node(self, label: frozenset | None = None) -> int:
        v = self._next_id
        self._next_id += 1
        self.children[v] = []
        self.parents[v] = []
        if label is not None:
            self.label[v] = label
        return v

    def add_edge(self, u: int, v: int) -> None:
        self.children[u].append(v)
        self.parents[v].append(u)

    def remove_edge(self, u: int, v: int) -> None:
        self.children[u].remove(v)
        self.parents[v].remove(u)

    def remove_node(self, v: int) -> None:
        for u in list(self.parents[v]):
            self.remove_edge(u, v)
        for w in list(self.children[v]):
            self.remove_edge(v, w)
        del self.children[v]
        del self.parents[v]
        self.label.pop(v, None)

    def subdivide(self, u: int, v: int) -> int:
        """Insert a new node in the middle of edge ``(u, v)`` and return it."""
        self.remove_edge(u, v)
        w = self.new_node()
        self.add_edge(u, w)
        self.add_edge(w, v)
        return w

    def copy(self) -> "Network":
        n = Network()
        n._next_id = self._next_id
        n.children = {v: list(cs) for v, cs in self.children.items()}
        n.parents = {v: list(ps) for v, ps in self.parents.items()}
        n.label = dict(self.label)
        return n

    # ------------------------------------------------------------------
    # basic queries
    # ------------------------------------------------------------------
    def nodes(self) -> list[int]:
        return list(self.children)

    def edges(self) -> list[tuple[int, int]]:
        return [(u, v) for u, cs in self.children.items() for v in cs]

    def roots(self) -> list[int]:
        return [v for v in self.children if not self.parents[v]]

    @property
    def root(self) -> int:
        rs = self.roots()
        if len(rs) != 1:
            raise ValueError(f"network has {len(rs)} roots, expected exactly 1")
        return rs[0]

    def leaves(self) -> list[int]:
        return [v for v in self.children if not self.children[v]]

    def is_reticulation(self, v: int) -> bool:
        return len(self.parents[v]) >= 2

    def reticulations(self) -> list[int]:
        return [v for v in self.children if len(self.parents[v]) >= 2]

    def reticulation_number(self) -> int:
        return sum(max(0, len(self.parents[v]) - 1) for v in self.children)

    def taxa(self) -> frozenset:
        out: frozenset = frozenset()
        for v in self.leaves():
            out |= self.label.get(v, DUMMY)
        return out

    def taxon_labels(self) -> list[frozenset]:
        """Leaf labels, dummies excluded."""
        return [self.label[v] for v in self.leaves() if self.label.get(v, DUMMY)]

    # ------------------------------------------------------------------
    # topological order
    # ------------------------------------------------------------------
    def topological_order(self) -> list[int]:
        indeg = {v: len(self.parents[v]) for v in self.children}
        queue = [v for v, d in indeg.items() if d == 0]
        order: list[int] = []
        while queue:
            v = queue.pop()
            order.append(v)
            for w in self.children[v]:
                indeg[w] -= 1
                if indeg[w] == 0:
                    queue.append(w)
        if len(order) != len(self.children):
            raise ValueError("network contains a directed cycle")
        return order

    # ------------------------------------------------------------------
    # softwired clusters
    # ------------------------------------------------------------------
    def hardwired_clusters(self) -> set[frozenset]:
        """Descendant taxon set of every node (the hardwired reading)."""
        below: dict[int, frozenset] = {}
        for v in reversed(self.topological_order()):
            s = self.label.get(v, DUMMY)
            for w in self.children[v]:
                s = s | below[w]
            below[v] = s
        return {s for s in below.values() if s}

    def clusters(self) -> set[frozenset]:
        """Clusters displayed by a *tree* (no reticulations): its edge clusters."""
        return self.hardwired_clusters()

    def switchings(self) -> Iterator[dict[int, int]]:
        """All ways of switching on exactly one incoming edge per reticulation."""
        rets = self.reticulations()
        if not rets:
            yield {}
            return
        for choice in product(*(self.parents[r] for r in rets)):
            yield dict(zip(rets, choice))

    def softwired_clusters(self) -> set[frozenset]:
        """Every cluster represented by ``self`` in the softwired sense.

        A cluster is represented if *some* displayed tree has a node whose
        descendant taxon set is exactly that cluster.  Different clusters may
        use different switchings.
        """
        found: set[frozenset] = set()
        order = self.topological_order()
        try:
            root = self.root
        except ValueError:
            return found
        for on in self.switchings():
            # nodes reachable from the root using switched-on edges only
            reachable = {root}
            for v in order:
                if v not in reachable:
                    continue
                for w in self.children[v]:
                    if len(self.parents[w]) >= 2 and on.get(w) != v:
                        continue
                    reachable.add(w)
            below: dict[int, frozenset] = {}
            for v in reversed(order):
                if v not in reachable:
                    continue
                s = self.label.get(v, DUMMY)
                for w in self.children[v]:
                    if w in reachable and not (
                        len(self.parents[w]) >= 2 and on.get(w) != v
                    ):
                        s = s | below[w]
                below[v] = s
                if s:
                    found.add(s)
        return found

    def represents(self, clusters: Iterable[frozenset]) -> bool:
        """Does the network represent every cluster in ``clusters``?"""
        wanted = {c for c in clusters if c}
        if not wanted:
            return True
        return wanted <= self.softwired_clusters()

    # ------------------------------------------------------------------
    # biconnected components and level
    # ------------------------------------------------------------------
    def _undirected_edges(self) -> list[tuple[int, int]]:
        return self.edges()

    def biconnected_components(self) -> list[list[tuple[int, int]]]:
        """Biconnected components of the underlying undirected multigraph.

        Returns a list of edge lists.  Iterative Hopcroft--Tarjan, with edges
        identified by index so parallel edges are handled correctly.
        """
        edges = self._undirected_edges()
        adj: dict[int, list[tuple[int, int]]] = {v: [] for v in self.children}
        for idx, (u, v) in enumerate(edges):
            adj[u].append((v, idx))
            adj[v].append((u, idx))

        disc: dict[int, int] = {}
        low: dict[int, int] = {}
        timer = 0
        components: list[list[tuple[int, int]]] = []
        edge_stack: list[int] = []

        for start in self.children:
            if start in disc:
                continue
            # stack frames: (node, parent_edge_index, iterator position)
            stack: list[tuple[int, int, int]] = [(start, -1, 0)]
            disc[start] = low[start] = timer
            timer += 1
            while stack:
                v, pe, i = stack.pop()
                if i < len(adj[v]):
                    stack.append((v, pe, i + 1))
                    w, eidx = adj[v][i]
                    if eidx == pe:
                        continue
                    if w not in disc:
                        edge_stack.append(eidx)
                        disc[w] = low[w] = timer
                        timer += 1
                        stack.append((w, eidx, 0))
                    elif disc[w] < disc[v]:
                        edge_stack.append(eidx)
                        low[v] = min(low[v], disc[w])
                else:
                    if stack:
                        parent = stack[-1][0]
                        low[parent] = min(low[parent], low[v])
                        if low[v] >= disc[parent]:
                            comp = []
                            while edge_stack and edge_stack[-1] != pe:
                                comp.append(edges[edge_stack.pop()])
                            if edge_stack:
                                comp.append(edges[edge_stack.pop()])
                            if comp:
                                components.append(comp)
        return components

    def level(self) -> int:
        """Maximum reticulation number over all biconnected components."""
        best = 0
        for comp in self.biconnected_components():
            verts = {x for e in comp for x in e}
            best = max(best, len(comp) - len(verts) + 1)
        return best

    def bridges(self) -> list[tuple[int, int]]:
        """Cut-edges: the biconnected components consisting of a single edge."""
        return [comp[0] for comp in self.biconnected_components() if len(comp) == 1]

    def is_simple(self) -> bool:
        """Is every cut-edge's head a leaf? (definition of a *simple* network)"""
        return all(not self.children[v] for _, v in self.bridges())

    # ------------------------------------------------------------------
    # tidying
    # ------------------------------------------------------------------
    def remove_dummy_leaves(self) -> None:
        for v in [v for v in self.leaves() if self.label.get(v, DUMMY) == DUMMY]:
            if v in self.children:
                self.remove_node(v)

    def tidy_up(self, keep_root_edge: bool = False) -> None:
        """Repeatedly apply the five clean-up rules from the paper.

        (1) delete unlabelled outdegree-0 nodes, (2) suppress indegree-1
        outdegree-1 nodes, (3) replace parallel edges by single edges,
        (4) drop an outdegree-1 root, (5) contract biconnected components with
        a single outgoing edge.
        """
        changed = True
        while changed:
            changed = False

            # (1) unlabelled dead ends
            for v in list(self.children):
                if v in self.children and not self.children[v] and v not in self.label:
                    self.remove_node(v)
                    changed = True

            # (3) parallel edges
            for u in list(self.children):
                if u not in self.children:
                    continue
                seen = set()
                for v in list(self.children[u]):
                    if v in seen:
                        self.remove_edge(u, v)
                        changed = True
                    else:
                        seen.add(v)

            # (2) suppress degree-2 nodes
            for v in list(self.children):
                if v not in self.children or v in self.label:
                    continue
                if len(self.parents[v]) == 1 and len(self.children[v]) == 1:
                    u = self.parents[v][0]
                    w = self.children[v][0]
                    self.remove_node(v)
                    if w not in self.children[u]:
                        self.add_edge(u, w)
                    changed = True

            # (4) outdegree-1 root
            if not keep_root_edge:
                for v in self.roots():
                    if v in self.children and len(self.children[v]) == 1 and v not in self.label:
                        self.remove_node(v)
                        changed = True
                        break

            # (5) biconnected component with a single outgoing edge
            if not changed:
                for comp in self.biconnected_components():
                    if len(comp) < 2:
                        continue
                    verts = {x for e in comp for x in e}
                    outgoing = [
                        (u, w)
                        for u in verts
                        for w in self.children[u]
                        if (u, w) not in comp and (w not in verts)
                    ]
                    if len(outgoing) == 1:
                        self._contract_component(verts, outgoing[0])
                        changed = True
                        break

    def _contract_component(self, verts: set[int], out_edge: tuple[int, int]) -> None:
        """Collapse a biconnected component into its single outgoing edge."""
        top = next((v for v in verts if any(p not in verts for p in self.parents[v])), None)
        _, below = out_edge
        incoming = []
        if top is not None:
            incoming = [p for p in self.parents[top] if p not in verts]
        for v in list(verts):
            if v in self.children:
                self.remove_node(v)
        for p in incoming:
            if below in self.children.get(p, []):
                continue
            self.add_edge(p, below)

    def contract_reticulation_edges(self) -> None:
        """Contract every edge whose tail and head are both reticulations."""
        changed = True
        while changed:
            changed = False
            for u, v in self.edges():
                if u in self.children and v in self.children:
                    if len(self.parents[u]) >= 2 and len(self.parents[v]) >= 2:
                        # merge u into v: v inherits u's parents
                        ups = list(self.parents[u])
                        downs = [w for w in self.children[u] if w != v]
                        self.remove_node(u)
                        for p in ups:
                            if v not in self.children[p]:
                                self.add_edge(p, v)
                        for w in downs:
                            if w not in self.children[v]:
                                self.add_edge(v, w)
                        changed = True
                        break

    def clean(self) -> None:
        """Full output clean-up: strip dummies, contract, tidy."""
        self.remove_dummy_leaves()
        self.tidy_up()
        self.contract_reticulation_edges()
        self.tidy_up()

    # ------------------------------------------------------------------
    # serialisation
    # ------------------------------------------------------------------
    def to_enewick(self, namer=None) -> str:
        """Extended Newick, with reticulations written as ``#H1`` and friends."""
        if namer is None:
            def namer(block: frozenset) -> str:
                return "+".join(sorted(str(t) for t in block))

        hybrid_id: dict[int, int] = {}
        for i, r in enumerate(sorted(self.reticulations()), start=1):
            hybrid_id[r] = i
        written: set[int] = set()

        def render(v: int) -> str:
            tag = ""
            if v in hybrid_id:
                tag = f"#H{hybrid_id[v]}"
                if v in written:
                    return tag
                written.add(v)
            lbl = self.label.get(v)
            name = namer(lbl) if lbl else ""
            if not self.children[v]:
                return f"{name}{tag}"
            inner = ",".join(render(w) for w in self.children[v])
            return f"({inner}){name}{tag}"

        return render(self.root) + ";"

    def __repr__(self) -> str:  # pragma: no cover - debugging aid
        return (
            f"<Network nodes={len(self.children)} edges={len(self.edges())} "
            f"retic={self.reticulation_number()} level={self.level()}>"
        )
