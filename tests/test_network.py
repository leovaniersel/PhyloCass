import pytest

from phylocass.network import DUMMY, Network
from phylocass.treebuild import build_tree


def fs(s):
    return frozenset(s)


def singletons(s):
    return [frozenset({t}) for t in s]


def shown(sets):
    return sorted(tuple(sorted(s)) for s in sets)


def leafy(net, name):
    """Add a leaf labelled ``name`` and return its node id."""
    return net.new_node(label=frozenset({name}))


class TestTreeStructure:
    def test_balanced_tree_clusters(self):
        tree = build_tree([fs("ab"), fs("cd")], singletons("abcd"))
        assert tree.reticulation_number() == 0
        assert tree.level() == 0
        assert shown(tree.clusters()) == [
            ("a",), ("a", "b"), ("a", "b", "c", "d"), ("b",), ("c",), ("c", "d"), ("d",)
        ]

    def test_incompatible_clusters_have_no_tree(self):
        assert build_tree([fs("ab"), fs("bc")], singletons("abc")) is None

    def test_nested_clusters_make_a_caterpillar(self):
        tree = build_tree([fs("ab"), fs("abc")], singletons("abcd"))
        assert tree.level() == 0
        assert fs("ab") in tree.clusters()
        assert fs("abc") in tree.clusters()

    def test_tree_has_no_extra_clusters(self):
        wanted = {fs("ab"), fs("abc")}
        tree = build_tree(wanted, singletons("abcd"))
        nontrivial = {c for c in tree.clusters() if 1 < len(c) < 4}
        assert nontrivial == wanted

    def test_blocks_become_single_leaves(self):
        tree = build_tree([fs("abc")], [fs("ab"), fs("c"), fs("d")])
        assert shown(tree.taxon_labels()) == [("a", "b"), ("c",), ("d",)]
        assert len(tree.leaves()) == 3


class TestLevelAndBiconnectedComponents:
    def _one_reticulation(self):
        """root -> (u, v); u -> a, r ; v -> b, r ; r -> c"""
        n = Network()
        root, u, v, r = n.new_node(), n.new_node(), n.new_node(), n.new_node()
        n.add_edge(root, u)
        n.add_edge(root, v)
        n.add_edge(u, leafy(n, "a"))
        n.add_edge(v, leafy(n, "b"))
        n.add_edge(u, r)
        n.add_edge(v, r)
        n.add_edge(r, leafy(n, "c"))
        return n

    def test_level_one(self):
        n = self._one_reticulation()
        assert n.reticulation_number() == 1
        assert n.level() == 1
        assert n.is_simple()

    def test_bridges_are_the_pendant_edges(self):
        n = self._one_reticulation()
        # only the three leaf edges are cut-edges; the rest form the blob
        assert len(n.bridges()) == 3
        for _, head in n.bridges():
            assert not n.children[head]

    def test_tree_is_level_zero(self):
        tree = build_tree([fs("ab"), fs("cd")], singletons("abcd"))
        assert tree.level() == 0
        assert tree.biconnected_components()  # every edge is its own component
        assert all(len(c) == 1 for c in tree.biconnected_components())

    def test_two_separate_blobs_stay_level_one(self):
        left = self._one_reticulation()
        right = self._one_reticulation()
        # graft right under a fresh root next to left
        n = Network()
        root = n.new_node()
        from phylocass.treebuild import graft

        a = n.new_node()
        b = n.new_node()
        n.add_edge(root, a)
        n.add_edge(root, b)
        graft(n, a, left, {})
        graft(n, b, right, {})
        assert n.reticulation_number() == 2
        assert n.level() == 1


class TestSoftwiredClusters:
    def test_reticulation_gives_two_displayed_trees(self):
        n = Network()
        root, u, v, r = n.new_node(), n.new_node(), n.new_node(), n.new_node()
        n.add_edge(root, u)
        n.add_edge(root, v)
        n.add_edge(u, leafy(n, "a"))
        n.add_edge(v, leafy(n, "b"))
        n.add_edge(u, r)
        n.add_edge(v, r)
        n.add_edge(r, leafy(n, "c"))

        soft = n.softwired_clusters()
        assert fs("ac") in soft  # switch on the u edge
        assert fs("bc") in soft  # switch on the v edge
        assert n.represents([fs("ac"), fs("bc")])

    def test_tree_softwired_equals_hardwired(self):
        tree = build_tree([fs("ab"), fs("abc")], singletons("abcd"))
        assert tree.softwired_clusters() == tree.hardwired_clusters()

    def test_represents_rejects_absent_cluster(self):
        tree = build_tree([fs("ab"), fs("cd")], singletons("abcd"))
        assert not tree.represents([fs("ac")])

    def test_dummy_leaves_contribute_nothing(self):
        n = Network()
        root = n.new_node()
        n.add_edge(root, leafy(n, "a"))
        n.add_edge(root, n.new_node(label=DUMMY))
        assert n.taxa() == fs("a")
        assert shown(n.hardwired_clusters()) == [("a",)]


class TestTidying:
    def test_dummy_leaf_removal_and_suppression(self):
        n = Network()
        root = n.new_node()
        mid = n.new_node()
        n.add_edge(root, mid)
        n.add_edge(mid, leafy(n, "a"))
        n.add_edge(mid, leafy(n, "b"))
        n.add_edge(root, n.new_node(label=DUMMY))
        n.clean()
        assert shown(n.taxon_labels()) == [("a",), ("b",)]
        assert n.reticulation_number() == 0

    def test_parallel_edges_collapse(self):
        n = Network()
        root, mid = n.new_node(), n.new_node()
        n.add_edge(root, mid)
        n.add_edge(root, mid)
        n.add_edge(mid, leafy(n, "a"))
        n.add_edge(mid, leafy(n, "b"))
        n.clean()
        assert len(n.edges()) == len(set(n.edges()))
        assert n.reticulation_number() == 0


class TestENewick:
    def test_tree_roundtrip_shape(self):
        from phylocass.newick import parse_newick

        tree = parse_newick("((a,b),(c,d));")
        text = tree.to_enewick()
        again = parse_newick(text)
        assert again.clusters() == tree.clusters()

    def test_reticulation_is_written_once(self):
        n = Network()
        root, u, v, r = n.new_node(), n.new_node(), n.new_node(), n.new_node()
        n.add_edge(root, u)
        n.add_edge(root, v)
        n.add_edge(u, leafy(n, "a"))
        n.add_edge(v, leafy(n, "b"))
        n.add_edge(u, r)
        n.add_edge(v, r)
        n.add_edge(r, leafy(n, "c"))
        text = n.to_enewick()
        assert text.count("#H1") == 2
        assert text.count("c") == 1
        assert text.endswith(";")
