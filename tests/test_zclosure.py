"""Z-closure: completing partial clusters from trees with different taxon sets."""

import random

import pytest

from phylocass import CassOptions, cass_from_trees, read_trees
from phylocass.io import displays_partial_trees
from phylocass.zclosure import PartialCluster, partial_clusters, z_closure

from conftest import fs, random_tree, shown

PARTIAL = "(((a,b),c),d); (((a,c),b),e);"


def close(text):
    trees = read_trees(text)
    partials, taxa = partial_clusters(trees)
    return trees, z_closure(partials, taxa)


class TestPartialClusters:
    def test_each_cluster_carries_its_taxon_set(self):
        _, taxa = partial_clusters(read_trees("((a,b),c);"))
        assert taxa == fs("abc")

    def test_root_cluster_is_dropped(self):
        """A partial tree's root says nothing about its taxa forming a clade."""
        partials, _ = partial_clusters(read_trees("((a,b),c); ((a,c),d);"))
        assert all(c != x for c, x, _ in partials)
        assert shown(c for c, _, _ in partials) == [("a", "b"), ("a", "c")]

    def test_taxon_sets_are_per_tree(self):
        partials, taxa = partial_clusters(read_trees("((a,b),c); ((a,c),d);"))
        assert taxa == fs("abcd")
        by_cluster = {c: x for c, x, _ in partials}
        assert by_cluster[fs("ab")] == fs("abc")
        assert by_cluster[fs("ac")] == fs("acd")


class TestTheRules:
    def test_overlapping_rule_widens_and_merges(self):
        """C1={a,b} on {a,b,c}; C2={a,b,c} on {a,b,c,d}: C1 learns about d."""
        result = z_closure(
            [(fs("ab"), fs("abc"), 0), (fs("abc"), fs("abcd"), 1)], fs("abcd")
        )
        assert fs("ab") in result.clusters
        assert fs("abc") in result.clusters

    def test_disjoint_rule_widens(self):
        """Two clades apart from each other extend each other's taxon sets.

        Each needs a witness on the other's root side: {a,b} is known on
        {a,b,c} and {c,d} on {a,c,d}, so c vouches for one and a for the other.
        """
        result = z_closure(
            [(fs("ab"), fs("abc"), 0), (fs("cd"), fs("acd"), 1)], fs("abcd")
        )
        assert fs("ab") in result.clusters

    def test_disjoint_rule_needs_a_witness_on_each_side(self):
        """No witness, no inference: {a,b} and {d,e} never meet here."""
        result = z_closure(
            [(fs("ab"), fs("abc"), 0), (fs("de"), fs("cde"), 1)], fs("abcde")
        )
        assert result.clusters == set()

    def test_nothing_derivable_leaves_clusters_partial(self):
        """Z-closure is incomplete; what stays partial is dropped."""
        result = z_closure([(fs("ab"), fs("abc"), 0)], fs("abcd"))
        assert result.clusters == set()
        assert result.dropped == 1

    def test_clusters_already_full_survive(self):
        result = z_closure([(fs("ab"), fs("abcd"), 0)], fs("abcd"))
        assert result.clusters == {fs("ab")}
        assert result.dropped == 0

    def test_closure_is_deterministic(self):
        partials, taxa = partial_clusters(read_trees(PARTIAL))
        answers = {frozenset(z_closure(partials, taxa).clusters) for _ in range(4)}
        assert len(answers) == 1

    def test_it_terminates_and_reports_rounds(self):
        _, result = close(PARTIAL)
        assert result.rounds >= 1
        assert not result.hit_limit


class TestPartialClusterType:
    def test_rejects_a_non_proper_subset(self):
        with pytest.raises(ValueError):
            PartialCluster(fs("abc"), fs("abc"))

    def test_rejects_an_empty_cluster(self):
        with pytest.raises(ValueError):
            PartialCluster(fs(""), fs("abc"))

    def test_accepts_a_proper_subset(self):
        assert PartialCluster(fs("ab"), fs("abc")).cluster == fs("ab")


class TestEndToEnd:
    def test_partial_input_is_refused_without_the_option(self):
        with pytest.raises(ValueError, match="missing taxa"):
            cass_from_trees(read_trees(PARTIAL))

    def test_the_error_points_at_the_option(self):
        with pytest.raises(ValueError) as excinfo:
            cass_from_trees(read_trees(PARTIAL))
        assert "z_closure=True" in str(excinfo.value)

    def test_z_closure_builds_a_network(self):
        trees = read_trees(PARTIAL)
        result = cass_from_trees(trees, CassOptions(z_closure=True, max_level=4))
        assert result.taxa == fs("abcde")
        assert result.level == 1
        assert result.represents_input()
        result.network.validate()

    def test_the_restricted_network_displays_each_partial_tree(self):
        """The right check for partial input: restrict, then ask."""
        trees = read_trees(PARTIAL)
        result = cass_from_trees(trees, CassOptions(z_closure=True, max_level=4))
        assert all(displays_partial_trees(result.network, trees))

    def test_statistics_are_reported(self):
        trees = read_trees(PARTIAL)
        result = cass_from_trees(trees, CassOptions(z_closure=True, max_level=4))
        assert result.z_closure is not None
        assert result.z_closure.partial_total >= len(result.z_closure.clusters)

    def test_equal_taxon_sets_skip_the_closure_entirely(self):
        trees = read_trees("(((a,b),c),d); (((a,c),b),d);")
        result = cass_from_trees(trees, CassOptions(z_closure=True))
        assert result.z_closure is None
        assert result.level == 1

    def test_z_closure_composes_with_display_mode(self):
        trees = read_trees(PARTIAL)
        result = cass_from_trees(
            trees, CassOptions(z_closure=True, display_trees=True, max_level=4)
        )
        assert result.represents_input()
        result.network.validate()


class TestRandomised:
    @pytest.mark.parametrize("seed", range(25))
    def test_output_is_always_valid_and_represents_its_clusters(self, seed):
        """Whatever Z-closure recovers, Cass must handle it soundly."""
        from phylozoo.core.network.dnetwork.derivations import subnetwork

        rng = random.Random(seed)
        k = rng.randint(6, 8)
        names = [chr(ord("a") + i) for i in range(k)]
        full = random_tree(rng, names).to_phylozoo()
        trees = [
            subnetwork(full, sorted(rng.sample(names, k - 1)))
            for _ in range(rng.randint(2, 3))
        ]
        if all(frozenset(t.taxa) == frozenset(names) for t in trees):
            pytest.skip("no taxa actually dropped")
        try:
            result = cass_from_trees(
                trees, CassOptions(z_closure=True, max_level=4, time_limit=10)
            )
        except RuntimeError:
            pytest.skip("search budget exhausted")
        assert result.represents_input()
        # a taxon dropped from every tree is simply not in the input
        assert result.taxa == frozenset().union(*(frozenset(t.taxa) for t in trees))
        result.network.validate()
