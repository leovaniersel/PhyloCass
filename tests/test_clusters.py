import pytest

from phylocass.clusters import (
    block_partition,
    collapse,
    incompatibility_graph,
    is_compatible,
    is_separated,
    is_st_set,
    maximal_st_sets,
    maximal_unseparated_sets,
    nontrivial_components,
    restrict,
    remove_taxa,
    round_up,
    unseparated_closure,
)

from conftest import fs, shown, singletons


class TestCompatibility:
    def test_disjoint_is_compatible(self):
        assert is_compatible(fs("ab"), fs("cd"))

    def test_nested_is_compatible(self):
        assert is_compatible(fs("ab"), fs("abc"))
        assert is_compatible(fs("abc"), fs("ab"))

    def test_equal_is_compatible(self):
        assert is_compatible(fs("ab"), fs("ab"))

    def test_overlap_is_incompatible(self):
        assert not is_compatible(fs("ab"), fs("bc"))

    def test_separated(self):
        assert is_separated(fs("ab"), [fs("bc")])
        assert not is_separated(fs("ab"), [fs("abc"), fs("de"), fs("a")])


class TestIncompatibilityGraph:
    def test_edges(self):
        clusters = [fs("ab"), fs("bc"), fs("de")]
        adj = incompatibility_graph(clusters)
        assert adj[0] == {1}
        assert adj[1] == {0}
        assert adj[2] == set()

    def test_nontrivial_components(self):
        clusters = [fs("ab"), fs("bc"), fs("de"), fs("ef"), fs("xy")]
        comps = nontrivial_components(clusters)
        assert shown([frozenset().union(*c) for c in comps]) == [
            ("a", "b", "c"),
            ("d", "e", "f"),
        ]

    def test_compatible_input_has_no_components(self):
        assert nontrivial_components([fs("ab"), fs("abc"), fs("de")]) == []


class TestRestrictAndRemove:
    def test_restrict(self):
        out = restrict([fs("abc"), fs("cd"), fs("ef")], fs("abcd"))
        assert shown(out) == [("a", "b", "c"), ("c", "d")]

    def test_restrict_drops_the_whole_set(self):
        assert restrict([fs("abc")], fs("abc")) == set()

    def test_remove_taxa(self):
        out = remove_taxa([fs("abc"), fs("cd")], fs("c"), fs("abde"))
        assert shown(out) == [("a", "b"), ("d",)]

    def test_remove_taxa_drops_ground(self):
        assert remove_taxa([fs("abc")], fs("c"), fs("ab")) == set()

    def test_round_up(self):
        blocks = [fs("ab"), fs("cd"), fs("e")]
        assert round_up(fs("a"), blocks) == fs("ab")
        assert round_up(fs("ac"), blocks) == fs("abcd")


class TestSTSets:
    def test_paper_collapse_example(self):
        """Section 4 of the paper: C = {{1,2},{2,3,4},{3,4}} collapses {3,4}."""
        clusters = [fs("12"), fs("234"), fs("34")]
        blocks = singletons("1234")
        assert shown(maximal_st_sets(clusters, blocks)) == [
            ("1",),
            ("2",),
            ("3", "4"),
        ]

        new_blocks, mapping, new_clusters = collapse(clusters, blocks)
        assert shown(new_blocks) == [("1",), ("2",), ("3", "4")]
        assert shown(new_clusters) == [("1", "2"), ("2", "3", "4")]
        assert shown(mapping[fs("34")]) == [("3",), ("4",)]

    def test_ground_set_is_not_an_st_set(self):
        assert not is_st_set(fs("abc"), [fs("ab")], fs("abc"))

    def test_st_set_needs_internal_compatibility(self):
        # {a,b,c} is unseparated but C|{a,b,c} = {{a,b},{b,c}} is not compatible
        clusters = [fs("ab"), fs("bc")]
        assert not is_separated(fs("abc"), clusters)
        assert not is_st_set(fs("abc"), clusters, fs("abcd"))

    def test_maximal_st_sets_partition_the_taxa(self):
        clusters = [fs("abc"), fs("cde"), fs("de")]
        parts = maximal_st_sets(clusters, singletons("abcde"))
        assert frozenset().union(*parts) == fs("abcde")
        assert sum(len(p) for p in parts) == 5

    def test_compatible_clusters_collapse_to_two_blocks(self):
        clusters = [fs("ab"), fs("abc")]
        parts = maximal_st_sets(clusters, singletons("abcd"))
        assert shown(parts) == [("a", "b", "c"), ("d",)]


class TestUnseparatedSets:
    def test_closure_absorbs_conflicting_clusters(self):
        # {a,b} conflicts with {b,c}, so the closure must swallow c too
        assert unseparated_closure(fs("ab"), [fs("bc")], fs("abcde")) == fs("abc")

    def test_closure_returns_none_when_it_reaches_the_ground_set(self):
        assert unseparated_closure(fs("ab"), [fs("bc")], fs("abc")) is None

    def test_greedy_pairwise_merging_would_stall_here(self):
        """Every pair inside {a,b,c} is separated, yet {a,b,c} is unseparated.

        A naive "merge any pair whose union is unseparated" loop returns four
        singletons; the closure-based merge finds {a,b,c}.
        """
        clusters = [fs("ab"), fs("ac")]
        ground = fs("abcd")
        for pair in (fs("ab"), fs("ac"), fs("bc")):
            assert is_separated(pair, clusters)
        assert not is_separated(fs("abc"), clusters)
        parts = maximal_unseparated_sets(clusters, singletons("abcd"))
        assert shown(parts) == [("a", "b", "c"), ("d",)]
        assert frozenset().union(*parts) == ground

    def test_unseparated_sets_are_proper_subsets(self):
        clusters = [fs("ab"), fs("bc")]
        parts = maximal_unseparated_sets(clusters, singletons("abc"))
        assert all(p != fs("abc") for p in parts)
        assert shown(parts) == [("a",), ("b",), ("c",)]

    def test_partition_property(self):
        clusters = [fs("abc"), fs("cd"), fs("ef")]
        parts = maximal_unseparated_sets(clusters, singletons("abcdefg"))
        assert frozenset().union(*parts) == fs("abcdefg")
        assert sum(len(p) for p in parts) == 7


class TestPhyloZooInterop:
    def test_blocks_convert_to_a_phylozoo_partition(self):
        clusters = [fs("abc"), fs("cd")]
        parts = maximal_unseparated_sets(clusters, singletons("abcde"))
        partition = block_partition(parts)
        assert partition.elements == fs("abcde")
        assert partition.size() == 5  # size() counts elements, not parts
        assert len(partition.parts) == len(parts)
        assert {frozenset(p) for p in partition.parts} == set(parts)
