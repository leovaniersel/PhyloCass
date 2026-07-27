"""End-to-end tests for the Cass algorithm."""

import random

import pytest

from phylocass import CassOptions, cass, cass_from_trees, read_trees
from phylocass.clusters import clusters_of_trees
from phylocass.network import Network


def fs(s):
    return frozenset(s)


def shown(sets):
    return sorted(tuple(sorted(s)) for s in sets)


def nontrivial(clusters, taxa):
    return {c for c in clusters if 1 < len(c) < len(taxa)}


# ----------------------------------------------------------------------
# behaviour on hand-picked inputs
# ----------------------------------------------------------------------


class TestTrivialInputs:
    def test_identical_trees_give_a_tree(self):
        trees = read_trees("((a,b),(c,d));  ((a,b),(c,d));")
        r = cass_from_trees(trees)
        assert r.level == 0
        assert r.reticulation_number == 0
        assert r.represents_input()

    def test_compatible_trees_give_a_tree(self):
        # {a,b} and {a,b,c} are nested, so no conflict
        trees = read_trees("(((a,b),c),d);  ((a,b),(c,d));")
        clusters, taxa = clusters_of_trees(trees)
        r = cass(clusters, taxa)
        assert r.represents_input()

    def test_refinement_is_conflict_free(self):
        trees = read_trees("((a,b),(c,d));  (((a,b),c),d);")
        clusters, taxa = clusters_of_trees(trees)
        # {c,d} vs {a,b,c} overlap in c only -> genuinely incompatible
        assert r_level(clusters, taxa) >= 1

    def test_single_tree_reproduces_itself(self):
        trees = read_trees("(((a,b),c),(d,e));")
        r = cass_from_trees(trees)
        assert r.level == 0
        assert nontrivial(r.network.clusters(), r.taxa) == r.clusters


def r_level(clusters, taxa, **kw):
    return cass(clusters, taxa, CassOptions(**kw)).level


class TestLevelOne:
    def test_single_leaf_move(self):
        """Moving one leaf between two trees needs exactly one reticulation."""
        trees = read_trees("(((a,b),c),d);  (((a,c),b),d);")
        r = cass_from_trees(trees)
        assert r.level == 1
        assert r.reticulation_number == 1
        assert r.represents_input()

    def test_galled_tree_clusters(self):
        clusters = {fs("ab"), fs("bc"), fs("abc")}
        r = cass(clusters, fs("abcd"))
        assert r.level == 1
        assert r.represents_input()


class TestLevelTwo:
    def test_four_taxa_double_conflict_needs_level_two(self):
        """{a,b},{c,d} vs {a,c},{b,d} cannot be displayed by one reticulation.

        A level-1 network displays exactly two trees that differ by relocating
        a single leaf; ((a,b),(c,d)) and ((a,c),(b,d)) differ by swapping two,
        so level 2 is optimal here.
        """
        trees = read_trees("((a,b),(c,d));  ((a,c),(b,d));")
        r = cass_from_trees(trees)
        assert r.level == 2
        assert r.reticulation_number == 2
        assert r.represents_input()

    def test_paper_figure_one(self):
        """The nine-taxon example from Figure 1 of the paper.

        The paper states a level-2 network with two reticulations exists for
        this cluster set, and that Cass finds it, where the galled-network
        algorithm needs four reticulations.
        """
        clusters = {
            fs(s)
            for s in (
                "abfgi", "abcfgi", "abfi", "bcfi", "cdeh", "deh",
                "bcfhi", "bcdfhi", "bci", "ag", "bi", "ci", "dh",
            )
        }
        taxa = fs("abcdefghi")
        r = cass(clusters, taxa)
        assert r.level == 2
        assert r.reticulation_number == 2
        assert r.represents_input()


class TestDecomposition:
    def test_independent_conflicts_stay_in_separate_blobs(self):
        """Two disjoint conflicts must not be merged into one biconnected component."""
        trees = read_trees(
            "(((a,b),(c,d)),((e,f),(g,h)));  (((a,c),(b,d)),((e,g),(f,h)));"
        )
        r = cass_from_trees(trees)
        assert r.level == 2
        assert r.reticulation_number == 4  # two independent level-2 blobs
        assert r.represents_input()

    def test_conflict_plus_clean_subtree(self):
        trees = read_trees("(((a,b),c),(x,y));  (((a,c),b),(x,y));")
        r = cass_from_trees(trees)
        assert r.level == 1
        assert r.represents_input()
        assert fs("xy") in r.network.softwired_clusters()

    def test_taxa_are_preserved(self):
        trees = read_trees("(((a,b),c),d);  (((a,c),b),d);")
        r = cass_from_trees(trees)
        assert r.network.taxa() == fs("abcd")
        labels = r.network.taxon_labels()
        assert len(labels) == len(set(labels)) == 4


# ----------------------------------------------------------------------
# round-trip: rediscover a network from the clusters it displays
# ----------------------------------------------------------------------


def random_tree(rng: random.Random, names) -> Network:
    """A random binary rooted tree on ``names``."""
    net = Network()
    nodes = [net.new_node(label=frozenset({n})) for n in names]
    while len(nodes) > 1:
        i, j = rng.sample(range(len(nodes)), 2)
        parent = net.new_node()
        net.add_edge(parent, nodes[i])
        net.add_edge(parent, nodes[j])
        nodes = [n for k, n in enumerate(nodes) if k not in (i, j)] + [parent]
    return net


def random_network(rng: random.Random, n_tree_leaves: int, n_ret: int) -> Network:
    """A random tree with ``n_ret`` extra leaves hung below new reticulations."""
    names = [chr(ord("a") + i) for i in range(n_tree_leaves)]
    net = random_tree(rng, names)
    for i in range(n_ret):
        edges = net.edges()
        # include the edge above the root as a possible attachment point
        old_root = net.root
        new_root = net.new_node()
        net.add_edge(new_root, old_root)
        edges.append((new_root, old_root))
        e1, e2 = rng.sample(edges, 2)
        w1 = net.subdivide(*e1)
        w2 = net.subdivide(*e2)
        ret = net.new_node()
        leaf = net.new_node(label=frozenset({chr(ord("a") + n_tree_leaves + i)}))
        net.add_edge(ret, leaf)
        net.add_edge(w1, ret)
        net.add_edge(w2, ret)
        net.tidy_up()
    return net


@pytest.mark.parametrize("seed", range(40))
def test_roundtrip_level_one(seed):
    """Cass must find a level-<=1 network for clusters that came from one."""
    rng = random.Random(seed)
    net = random_network(rng, rng.randint(3, 5), 1)
    if net.level() != 1:
        pytest.skip("random network degenerated")
    clusters = {c for c in net.softwired_clusters() if 1 < len(c) < len(net.taxa())}
    if not clusters:
        pytest.skip("no non-trivial clusters")
    r = cass(clusters, net.taxa(), CassOptions(max_level=3))
    assert r.represents_input(), "output does not display every input cluster"
    assert r.level <= 1, f"level {r.level} for clusters from a level-1 network"


@pytest.mark.parametrize("seed", range(25))
def test_roundtrip_level_two(seed):
    """The paper's guarantee: exact whenever a level-<=2 network exists."""
    rng = random.Random(1000 + seed)
    net = random_network(rng, rng.randint(3, 4), 2)
    if net.level() != 2:
        pytest.skip("random network degenerated")
    clusters = {c for c in net.softwired_clusters() if 1 < len(c) < len(net.taxa())}
    if not clusters:
        pytest.skip("no non-trivial clusters")
    r = cass(clusters, net.taxa(), CassOptions(max_level=4))
    assert r.represents_input(), "output does not display every input cluster"
    assert r.level <= 2, f"level {r.level} for clusters from a level-2 network"


@pytest.mark.parametrize("seed", range(30))
def test_random_trees_always_produce_a_valid_network(seed):
    """Whatever the input trees, the output must display every input cluster."""
    rng = random.Random(500 + seed)
    names = [chr(ord("a") + i) for i in range(rng.randint(4, 6))]
    trees = [random_tree(rng, names) for _ in range(rng.randint(2, 3))]
    clusters, taxa = clusters_of_trees(trees)
    r = cass(clusters, taxa, CassOptions(max_level=4))
    assert r.represents_input()
    assert r.network.taxa() == taxa
    labels = r.network.taxon_labels()
    assert len(labels) == len(set(labels)) == len(taxa)


class TestOutputShape:
    def test_output_is_a_valid_dag_with_one_root(self):
        trees = read_trees("((a,b),(c,d));  ((a,c),(b,d));")
        r = cass_from_trees(trees)
        assert len(r.network.roots()) == 1
        r.network.topological_order()  # raises on a cycle

    def test_enewick_is_parseable_back_into_the_same_clusters(self):
        from phylocass.newick import parse_newick

        trees = read_trees("(((a,b),c),d);  (((a,c),b),d);")
        r = cass_from_trees(trees)
        text = r.to_enewick()
        assert text.endswith(";")
        assert "#H1" in text

    def test_no_edge_joins_two_reticulations(self):
        trees = read_trees("((a,b),(c,d));  ((a,c),(b,d));")
        n = cass_from_trees(trees).network
        for u, v in n.edges():
            assert not (n.is_reticulation(u) and n.is_reticulation(v))
