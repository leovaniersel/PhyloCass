"""The working graph, and its agreement with PhyloZoo.

The value of these tests is the cross-checking: PhyloCass computes softwired
clusters and level itself, because mid-search graphs are not valid
``DirectedPhyNetwork``s and because it is the hot loop. Those routines are
checked here against PhyloZoo's own ``displayed_trees`` and ``level`` on
finished networks, which is a genuinely independent path.
"""

import random

import pytest
from phylozoo.core.network.dnetwork.classifications import (
    level as pz_level,
    reticulation_number as pz_reticulation_number,
)

from phylocass.io import hardwired_clusters, softwired_clusters
from phylocass.treebuild import build_tree
from phylocass.workgraph import DUMMY, WorkGraph

from conftest import fs, random_network, random_tree, shown, singletons


class TestConstruction:
    def test_subdivide_inserts_a_node(self):
        g = WorkGraph()
        a, b = g.new_node(), g.new_node(block=fs("x"))
        g.add_edge(a, b)
        w = g.subdivide(a, b)
        assert sorted(g.edges()) == sorted([(a, w), (w, b)])

    def test_copy_is_independent(self):
        g = WorkGraph()
        a, b = g.new_node(), g.new_node(block=fs("x"))
        g.add_edge(a, b)
        h = g.copy()
        h.add_edge(a, h.new_node(block=fs("y")))
        assert len(g.edges()) == 1
        assert len(h.edges()) == 2

    def test_root_and_leaves(self):
        tree = build_tree([fs("ab")], singletons("abc"))
        assert tree.g.indegree(tree.root) == 0
        assert shown(tree.blocks()) == [("a",), ("b",), ("c",)]

    def test_multiple_roots_rejected(self):
        g = WorkGraph()
        g.new_node()
        g.new_node()
        with pytest.raises(ValueError, match="2 roots"):
            _ = g.root


class TestSoftwiredClusters:
    def _one_reticulation(self):
        """root -> (u, v); u -> a, r; v -> b, r; r -> c"""
        g = WorkGraph()
        root, u, v, r = g.new_node(), g.new_node(), g.new_node(), g.new_node()
        g.add_edge(root, u)
        g.add_edge(root, v)
        g.add_edge(u, g.new_node(block=fs("a")))
        g.add_edge(v, g.new_node(block=fs("b")))
        g.add_edge(u, r)
        g.add_edge(v, r)
        g.add_edge(r, g.new_node(block=fs("c")))
        return g

    def test_both_switchings_are_seen(self):
        g = self._one_reticulation()
        soft = g.softwired_clusters()
        assert fs("ac") in soft
        assert fs("bc") in soft
        assert g.represents([fs("ac"), fs("bc")])

    def test_tree_softwired_equals_its_clusters(self):
        tree = build_tree([fs("ab"), fs("abc")], singletons("abcd"))
        assert fs("ab") in tree.softwired_clusters()
        assert fs("abc") in tree.softwired_clusters()

    def test_represents_rejects_absent_cluster(self):
        tree = build_tree([fs("ab"), fs("cd")], singletons("abcd"))
        assert not tree.represents([fs("ac")])

    def test_dummy_leaves_contribute_nothing(self):
        g = WorkGraph()
        root = g.new_node()
        g.add_edge(root, g.new_node(block=fs("a")))
        g.add_edge(root, g.new_node(block=DUMMY))
        assert g.taxa() == fs("a")
        assert shown(g.softwired_clusters()) == [("a",)]


class TestAgreementWithPhyloZoo:
    @pytest.mark.parametrize("seed", range(30))
    def test_softwired_clusters_match_displayed_trees(self, seed):
        rng = random.Random(seed)
        work = random_network(rng, rng.randint(3, 5), rng.randint(0, 2))
        work.clean()
        network = work.to_phylozoo()
        assert work.softwired_clusters() == softwired_clusters(network)

    @pytest.mark.parametrize("seed", range(30))
    def test_level_and_reticulation_number_match(self, seed):
        rng = random.Random(100 + seed)
        work = random_network(rng, rng.randint(3, 5), rng.randint(0, 2))
        work.clean()
        network = work.to_phylozoo()
        assert work.level() == pz_level(network)
        assert work.reticulation_number() == pz_reticulation_number(network)

    @pytest.mark.parametrize("seed", range(20))
    def test_cleaned_graphs_are_valid_phylozoo_networks(self, seed):
        rng = random.Random(200 + seed)
        work = random_network(rng, rng.randint(3, 5), rng.randint(0, 2))
        work.clean()
        work.to_phylozoo().validate()

    def test_tree_hardwired_clusters_match(self):
        rng = random.Random(7)
        work = random_tree(rng, list("abcde"))
        network = work.to_phylozoo()
        assert work.softwired_clusters() == hardwired_clusters(network)


class TestTidying:
    def test_dummy_leaf_removal_and_suppression(self):
        g = WorkGraph()
        root, mid = g.new_node(), g.new_node()
        g.add_edge(root, mid)
        g.add_edge(mid, g.new_node(block=fs("a")))
        g.add_edge(mid, g.new_node(block=fs("b")))
        g.add_edge(root, g.new_node(block=DUMMY))
        g.clean()
        assert shown(g.blocks()) == [("a",), ("b",)]
        assert g.reticulation_number() == 0

    def test_parallel_edges_are_identified(self):
        g = WorkGraph()
        root, mid = g.new_node(), g.new_node()
        g.add_edge(root, mid)
        g.add_edge(root, mid)
        g.add_edge(mid, g.new_node(block=fs("a")))
        g.add_edge(mid, g.new_node(block=fs("b")))
        g.clean()
        assert len(g.edges()) == len(set(g.edges()))
        assert g.reticulation_number() == 0

    def test_reticulation_chain_contracts_to_indegree_three(self):
        """Two stacked reticulations become one of indegree 3.

        This is the shape the dummy taxon leaves behind once ``d`` is removed.
        """
        g = WorkGraph()
        root = g.new_node()
        p1, p2, p3 = g.new_node(), g.new_node(), g.new_node()
        r1, r2 = g.new_node(), g.new_node()
        for p in (p1, p2, p3):
            g.add_edge(root, p)
            g.add_edge(p, g.new_node(block=fs(chr(ord("w") + (p % 3)) * 1)))
        g.add_edge(p1, r1)
        g.add_edge(p2, r1)
        g.add_edge(r1, r2)
        g.add_edge(p3, r2)
        g.add_edge(r2, g.new_node(block=fs("z")))
        g.contract_reticulation_edges()
        survivors = [v for v in g.nodes() if g.g.indegree(v) >= 2]
        assert len(survivors) == 1
        assert g.g.indegree(survivors[0]) == 3


class TestSimplicity:
    def test_pendant_edges_are_the_only_cut_edges_in_a_blob(self):
        g = WorkGraph()
        root, u, v, r = g.new_node(), g.new_node(), g.new_node(), g.new_node()
        g.add_edge(root, u)
        g.add_edge(root, v)
        g.add_edge(u, g.new_node(block=fs("a")))
        g.add_edge(v, g.new_node(block=fs("b")))
        g.add_edge(u, r)
        g.add_edge(v, r)
        g.add_edge(r, g.new_node(block=fs("c")))
        assert g.level() == 1
        assert g.is_simple()

    def test_a_tree_is_not_simple_unless_trivial(self):
        tree = build_tree([fs("ab")], singletons("abc"))
        # the edge into the {a,b} node is a cut-edge whose head is not a leaf
        assert not tree.is_simple()
        assert tree.level() == 0
