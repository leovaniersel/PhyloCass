"""Display mode: require the network to display the input trees themselves.

Cluster mode asks that every input cluster appear in *some* switching, which
is what the paper describes.  Display mode asks that each input tree's clusters
all appear in *one* switching, so that switching's tree displays the input
tree.  That makes the reticulation number an upper bound on the hybridization
number of the input trees.
"""

import random

import pytest

from phylocass import CassOptions, cass, cass_from_trees, read_trees
from phylocass.clusters import project
from phylocass.io import clusters_of_trees, displays_trees, per_tree_clusters
from phylocass.workgraph import WorkGraph

from conftest import fs, random_tree, shown, singletons

# three trees whose clusters fit a level-2 network that does not display them
THREE_TREES = "(((a,b),c),(d,e)); (((a,c),b),(d,e)); (((a,d),b),(c,e));"


def both_modes(text, **kw):
    trees = read_trees(text)
    opts = dict(max_level=6, time_limit=60)
    opts.update(kw)
    return (
        cass_from_trees(trees, CassOptions(display_trees=False, **opts)),
        cass_from_trees(trees, CassOptions(display_trees=True, **opts)),
    )


class TestTheDistinction:
    def test_cluster_mode_can_miss_the_trees_and_display_mode_cannot(self):
        """The case the option exists for."""
        cluster_mode, display_mode = both_modes(THREE_TREES)

        # both represent every input cluster
        assert cluster_mode.represents_input()
        assert display_mode.represents_input()

        # but only display mode guarantees the trees themselves
        assert cluster_mode.displays_input_trees() is False
        assert display_mode.displays_input_trees() is True

        # and that costs reticulations
        assert display_mode.reticulation_number > cluster_mode.reticulation_number

    def test_reported_numbers_for_that_case(self):
        cluster_mode, display_mode = both_modes(THREE_TREES)
        assert (cluster_mode.level, cluster_mode.reticulation_number) == (2, 2)
        assert (display_mode.level, display_mode.reticulation_number) == (3, 3)

    def test_when_cluster_mode_already_displays_nothing_changes(self):
        cluster_mode, display_mode = both_modes("(((a,b),c),d); (((a,c),b),d);")
        assert cluster_mode.displays_input_trees() is True
        assert display_mode.reticulation_number == cluster_mode.reticulation_number
        assert display_mode.level == cluster_mode.level

    def test_identical_trees_still_give_a_tree(self):
        _, display_mode = both_modes("((a,b),(c,d)); ((a,b),(c,d));")
        assert display_mode.level == 0
        assert display_mode.reticulation_number == 0
        assert display_mode.network.is_tree()
        assert display_mode.displays_input_trees() is True


class TestTwoBinaryTrees:
    """van Iersel & Kelk 2011: for two binary trees on the same taxon set the
    minimum reticulation number and level are the same in the tree model and
    the softwired-cluster model.  That is about the minimum, not about any
    particular network -- a minimum cluster network need not display the trees.
    """

    TWO = "((a,e),(c,(b,d))); (d,(b,(c,(a,e))));"

    def test_the_minimum_coincides_so_display_mode_is_free(self):
        cluster_mode, display_mode = both_modes(self.TWO)
        assert display_mode.reticulation_number == cluster_mode.reticulation_number
        assert display_mode.level == cluster_mode.level

    def test_but_cluster_mode_still_misses_the_trees_here(self):
        cluster_mode, display_mode = both_modes(self.TWO)
        assert cluster_mode.displays_input_trees() is False
        assert display_mode.displays_input_trees() is True

    @pytest.mark.parametrize("seed", range(25))
    def test_display_mode_is_free_on_random_pairs(self, seed):
        rng = random.Random(seed)
        names = [chr(ord("a") + i) for i in range(rng.randint(4, 7))]
        trees = [random_tree(rng, names).to_phylozoo() for _ in range(2)]
        opts = dict(max_level=5, time_limit=10)
        try:
            cluster_mode = cass_from_trees(trees, CassOptions(**opts))
            display_mode = cass_from_trees(
                trees, CassOptions(display_trees=True, **opts)
            )
        except RuntimeError:
            pytest.skip("search budget exhausted")
        assert display_mode.reticulation_number == cluster_mode.reticulation_number
        assert display_mode.level == cluster_mode.level
        assert display_mode.displays_input_trees() is True


class TestApi:
    def test_display_trees_needs_provenance(self):
        with pytest.raises(ValueError, match="which cluster came from which tree"):
            cass({fs("ab"), fs("bc")}, fs("abc"), CassOptions(display_trees=True))

    def test_explicit_tree_clusters_are_accepted(self):
        result = cass(
            {fs("ab"), fs("ac")},
            fs("abcd"),
            CassOptions(display_trees=True, max_level=4),
            tree_clusters=[{fs("ab")}, {fs("ac")}],
        )
        assert result.displays_input_trees() is True

    def test_displays_input_trees_is_none_without_provenance(self):
        result = cass({fs("ab"), fs("bc")}, fs("abcd"))
        assert result.displays_input_trees() is None

    def test_cass_from_trees_records_provenance_even_in_cluster_mode(self):
        result = cass_from_trees(read_trees("(((a,b),c),d); (((a,c),b),d);"))
        assert result.tree_clusters is not None
        assert result.displays_input_trees() in (True, False)

    def test_per_tree_clusters_keeps_trees_apart(self):
        trees = read_trees("((a,b),(c,d)); ((a,c),(b,d));")
        clusters, taxa = clusters_of_trees(trees)
        per_tree = per_tree_clusters(trees, taxa)
        assert len(per_tree) == 2
        assert shown(per_tree[0]) == [("a", "b"), ("c", "d")]
        assert shown(per_tree[1]) == [("a", "c"), ("b", "d")]
        assert set().union(*per_tree) == clusters


class TestWorkGraphDisplays:
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

    def test_switchings_are_kept_apart(self):
        g = self._one_reticulation()
        per_switching = g.switching_cluster_sets()
        assert len(per_switching) == 2
        assert g.softwired_clusters() == set().union(*per_switching)

    def test_two_clusters_from_one_tree_needing_different_switchings(self):
        g = self._one_reticulation()
        # {a,c} lives in one switching and {b,c} in the other
        assert g.represents([fs("ac"), fs("bc")])
        assert not g.displays([[fs("ac"), fs("bc")]])
        # but as two separate trees they are fine
        assert g.displays([[fs("ac")], [fs("bc")]])

    def test_displays_implies_represents(self):
        g = self._one_reticulation()
        assert g.displays([[fs("ac")]])
        assert g.represents([fs("ac")])

    def test_empty_requirements_are_satisfied(self):
        g = self._one_reticulation()
        assert g.displays([])
        assert g.displays([[]])


class TestProject:
    def test_drops_clusters_that_became_a_single_block(self):
        blocks = [fs("ab"), fs("c"), fs("d")]
        assert project([fs("a")], blocks) == set()

    def test_drops_the_whole_taxon_set(self):
        blocks = [fs("a"), fs("b")]
        assert project([fs("ab")], blocks) == set()

    def test_rounds_up_to_whole_blocks(self):
        blocks = [fs("ab"), fs("c"), fs("d")]
        assert project([fs("ac")], blocks) == {fs("abc")}


class TestRandomised:
    @pytest.mark.parametrize("seed", range(20))
    def test_display_mode_output_really_displays_the_trees(self, seed):
        """Verified through PhyloZoo, independently of the search."""
        rng = random.Random(3000 + seed)
        names = [chr(ord("a") + i) for i in range(rng.randint(4, 5))]
        trees = [random_tree(rng, names).to_phylozoo() for _ in range(rng.randint(2, 3))]
        try:
            result = cass_from_trees(
                trees, CassOptions(max_level=5, time_limit=10, display_trees=True)
            )
        except RuntimeError:
            pytest.skip("search budget exhausted")
        assert result.displays_input_trees() is True
        assert result.represents_input()
        result.network.validate()

    @pytest.mark.parametrize("seed", range(15))
    def test_display_mode_never_reaches_a_lower_level(self, seed):
        """Display mode is strictly stronger, so it cannot do better."""
        rng = random.Random(6000 + seed)
        names = [chr(ord("a") + i) for i in range(rng.randint(4, 5))]
        trees = [random_tree(rng, names).to_phylozoo() for _ in range(rng.randint(2, 3))]
        opts = dict(max_level=5, time_limit=10)
        try:
            cluster_mode = cass_from_trees(trees, CassOptions(**opts))
            display_mode = cass_from_trees(
                trees, CassOptions(display_trees=True, **opts)
            )
        except RuntimeError:
            pytest.skip("search budget exhausted")
        assert display_mode.level >= cluster_mode.level

    @pytest.mark.parametrize("seed", range(10))
    def test_multiple_incompatibility_components(self, seed):
        """Each component gets its own switching; they compose into one global one."""
        rng = random.Random(8000 + seed)
        left, right = list("abcd"), list("efg")
        texts = []
        for _ in range(rng.randint(2, 3)):
            lt = random_tree(rng, left).to_phylozoo().to_string().rstrip(";")
            rt = random_tree(rng, right).to_phylozoo().to_string().rstrip(";")
            texts.append(f"({lt},{rt});")
        trees = [read_trees(t)[0] for t in texts]
        try:
            result = cass_from_trees(
                trees, CassOptions(max_level=5, time_limit=15, display_trees=True)
            )
        except RuntimeError:
            pytest.skip("search budget exhausted")
        assert result.displays_input_trees() is True
        result.network.validate()


class TestDeterminism:
    def test_display_mode_is_deterministic(self):
        answers = {
            cass_from_trees(
                read_trees(THREE_TREES), CassOptions(display_trees=True, max_level=5)
            ).to_enewick()
            for _ in range(3)
        }
        assert len(answers) == 1
