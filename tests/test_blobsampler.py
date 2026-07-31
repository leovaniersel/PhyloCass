"""The generator-based sampler, and Cass round-tripped through it.

`conftest.random_network` builds a network by hanging a fresh leaf below each
new reticulation: always one blob, and almost every reticulation directly above
a leaf, which is the shape Cass's own construction produces. This sampler
expands a tree of blobs built from PhyloZoo's level-k generators instead, so
networks have several blobs of mixed level, reticulations sit above whole
subtrees, and blobs nest below other blobs' reticulations.
"""

import random

import pytest
from phylozoo.core.network.dnetwork.classifications import level as pz_level
from phylozoo.core.network.dnetwork.generator.side import HybridSide
from phylozoo.core.primitives.d_multigraph.features import biconnected_components

from phylocass import CassOptions, cass

from blobsampler import (
    generator_catalogue,
    min_leaves,
    random_blob_network,
    sample_blob,
    sample_simple_network,
)


def blob_reticulation_numbers(net):
    """Reticulation number of every non-trivial biconnected component."""
    edges = net.edges()
    out = []
    for comp in biconnected_components(net.g):
        inside = sum(1 for u, v in edges if u in comp and v in comp)
        r = inside - len(comp) + 1
        if r > 0:
            out.append(r)
    return out


def blob_nodes(net):
    edges = net.edges()
    inside = set()
    for comp in biconnected_components(net.g):
        m = sum(1 for u, v in edges if u in comp and v in comp)
        if m - len(comp) + 1 > 0:
            inside |= set(comp)
    return inside


class TestGeneratorCatalogue:
    @pytest.mark.parametrize(
        "level, count", [(1, 1), (2, 4), (3, 65)]
    )
    def test_published_counts(self, level, count):
        """PhyloZoo enumerates the generators; these are the known counts."""
        assert len(generator_catalogue(level)) == count

    def test_level_zero_has_no_blob_generators(self):
        assert generator_catalogue(0) == ()

    def test_catalogue_order_is_stable(self):
        assert [str(g.sides) for g in generator_catalogue(2)] == [
            str(g.sides) for g in generator_catalogue(2)
        ]

    @pytest.mark.parametrize("level", [1, 2, 3])
    def test_min_leaves_counts_the_hybrid_sides(self, level):
        for generator in generator_catalogue(level):
            hybrids = sum(1 for s in generator.sides if isinstance(s, HybridSide))
            assert min_leaves(generator) == max(2, hybrids)


class TestSampleBlob:
    @pytest.mark.parametrize("level", [1, 2, 3])
    def test_blob_has_the_requested_level_and_leaf_count(self, level):
        rng = random.Random(0)
        for generator in generator_catalogue(level):
            wanted = min_leaves(generator) + 2
            blob = sample_blob(generator, wanted, rng)
            assert blob is not None, f"could not build a blob from {generator.sides}"
            assert pz_level(blob) == level
            assert len(blob.taxa) == wanted
            assert len(list(blob.edges)) == len(set(blob.edges)), "no parallel edges"
            blob.validate()

    def test_too_few_leaves_is_refused(self):
        (generator,) = generator_catalogue(1)
        assert sample_blob(generator, 1, random.Random(0)) is None


class TestSimpleNetworks:
    """One blob carrying every taxon -- the shape that is actually hard."""

    @pytest.mark.parametrize("level", [1, 2, 3])
    @pytest.mark.parametrize("n_taxa", [5, 12])
    def test_shape(self, level, n_taxa):
        net = sample_simple_network(random.Random(level * 10 + n_taxa), level, n_taxa)
        assert net is not None
        assert len(net.taxa) == n_taxa
        assert pz_level(net) == level
        net.validate()

    @pytest.mark.parametrize("level", [1, 2])
    def test_nothing_collapses_so_the_component_is_the_whole_dataset(self, level):
        """Contrast with the blob tree, whose components stay tiny."""
        from phylocass.clusters import (
            canonical,
            maximal_unseparated_sets,
            nontrivial_components,
        )
        from phylocass.io import softwired_clusters

        net = sample_simple_network(random.Random(3), level, 12)
        taxa = frozenset(net.taxa)
        clusters = {c for c in softwired_clusters(net) if 1 < len(c) < len(taxa)}
        biggest = 0
        for comp in nontrivial_components(canonical(clusters)):
            support = frozenset().union(*comp)
            blocks = maximal_unseparated_sets(
                comp, [frozenset({t}) for t in sorted(support, key=str)]
            )
            biggest = max(biggest, len(blocks))
        assert biggest >= 10, "a simple network should barely collapse"

    def test_round_trips_through_cass(self):
        from phylocass.io import softwired_clusters

        net = sample_simple_network(random.Random(0), 2, 12)
        taxa = frozenset(net.taxa)
        clusters = {c for c in softwired_clusters(net) if 1 < len(c) < len(taxa)}
        r = cass(clusters, taxa, CassOptions(max_level=3, time_limit=120))
        assert r.represents_input()
        assert r.level <= 2


class TestSampledNetworks:
    @pytest.mark.parametrize("seed", range(15))
    def test_shape_is_valid(self, seed):
        net = random_blob_network(random.Random(seed), 20, max_level=2)
        assert len(net.taxa()) == 20
        assert net.level() <= 2
        net.to_phylozoo().validate()

    def test_networks_have_several_blobs(self):
        counts = [
            len(blob_reticulation_numbers(random_blob_network(random.Random(s), 20)))
            for s in range(40)
        ]
        assert max(counts) >= 4
        assert sum(c > 1 for c in counts) > 30

    def test_blobs_are_level_one_and_level_two(self):
        seen = set()
        for s in range(40):
            seen.update(blob_reticulation_numbers(random_blob_network(random.Random(s), 20)))
        assert seen == {1, 2}

    def test_level_three_blobs_can_be_sampled(self):
        seen = set()
        for s in range(30):
            net = random_blob_network(random.Random(s), 24, max_level=3)
            seen.update(blob_reticulation_numbers(net))
            assert net.level() <= 3
        assert 3 in seen, "expected at least one level-3 blob"

    def test_reticulations_sit_above_subtrees_not_only_leaves(self):
        """The coverage gap in the tree-plus-reticulations sampler."""
        above_leaf = above_subtree = 0
        for s in range(40):
            net = random_blob_network(random.Random(s), 20)
            children, parents = net.adjacency()
            for v in net.nodes():
                if len(parents[v]) >= 2:
                    for c in children[v]:
                        if children[c]:
                            above_subtree += 1
                        else:
                            above_leaf += 1
        assert above_leaf > 0
        assert above_subtree > 0.2 * (above_leaf + above_subtree)

    def test_blobs_occur_directly_below_reticulations(self):
        """Nesting the tree-plus-reticulations sampler cannot produce at all."""
        with_nesting = 0
        for s in range(40):
            net = random_blob_network(random.Random(s), 20)
            children, parents = net.adjacency()
            inside = blob_nodes(net)
            if any(
                c in inside
                for v in net.nodes()
                if len(parents[v]) >= 2
                for c in children[v]
            ):
                with_nesting += 1
        assert with_nesting > 30, f"only {with_nesting}/40 nested a blob under a reticulation"


class TestRoundTrip:
    @pytest.mark.parametrize("seed", range(20))
    def test_cass_matches_or_beats_the_source_level(self, seed):
        net = random_blob_network(random.Random(seed), 20, max_level=2)
        source_level = net.level()
        taxa = net.taxa()
        clusters = {c for c in net.softwired_clusters() if 1 < len(c) < len(taxa)}
        if not clusters:
            pytest.skip("no non-trivial clusters")
        r = cass(clusters, taxa, CassOptions(max_level=4, time_limit=60))
        assert r.represents_input()
        assert r.level <= source_level
        assert frozenset(r.network.taxa) == taxa
        r.network.validate()

    @pytest.mark.parametrize("seed", range(10))
    def test_galled_networks_stay_level_one(self, seed):
        net = random_blob_network(random.Random(100 + seed), 20, max_level=1)
        taxa = net.taxa()
        clusters = {c for c in net.softwired_clusters() if 1 < len(c) < len(taxa)}
        if not clusters:
            pytest.skip("no non-trivial clusters")
        r = cass(clusters, taxa, CassOptions(max_level=4, time_limit=60))
        assert r.represents_input()
        assert r.level <= 1
