"""The blob-wise shortcut must agree with plain enumeration.

Softwired clusters are computed one biconnected component at a time rather than
by enumerating every switching of the whole network, which turns
``2 ** (total reticulations)`` into a sum of ``2 ** (reticulations in one
blob)``. That is only worth having if it gives exactly the same answer, so
these tests pin it against PhyloZoo's ``displayed_trees``.
"""

import random

import pytest

from phylocass.io import (
    displays_trees,
    displays_trees_via_displayed_trees,
    hardwired_clusters,
    softwired_clusters,
    softwired_clusters_via_displayed_trees,
)
from phylocass.workgraph import WorkGraph

from blobsampler import random_blob_network
from conftest import fs, random_network, random_tree

# keep the reference affordable: it enumerates 2 ** reticulations trees
REFERENCE_CAP = 10


def sample(seed, taxa=(6, 12), max_level=2):
    rng = random.Random(seed)
    if seed % 3 == 0:
        work = random_network(rng, rng.randint(3, 6), rng.randint(0, 2))
        work.clean()
    else:
        work = random_blob_network(rng, rng.randint(*taxa), max_level=max_level)
    return work


class TestSoftwiredClusters:
    @pytest.mark.parametrize("seed", range(45))
    def test_agrees_with_displayed_trees(self, seed):
        work = sample(seed)
        if work.reticulation_number() > REFERENCE_CAP:
            pytest.skip("reference would be too slow")
        network = work.to_phylozoo()
        assert softwired_clusters(network) == softwired_clusters_via_displayed_trees(
            network
        )

    @pytest.mark.parametrize("seed", range(30))
    def test_workgraph_agrees_with_its_own_switching_enumeration(self, seed):
        work = sample(seed, max_level=3)
        if work.reticulation_number() > REFERENCE_CAP:
            pytest.skip("reference would be too slow")
        whole = set()
        for one in work.switching_cluster_sets():
            whole |= one
        assert work.softwired_clusters() == whole

    def test_a_tree_is_unaffected(self):
        work = random_tree(random.Random(0), list("abcdef"))
        assert work.softwired_clusters() == set(work.blob_cluster_sets()[0])

    def test_fixed_clusters_are_the_switching_invariant_ones(self):
        """Taxa below a cut-edge cannot be cut off by any switching."""
        work = random_blob_network(random.Random(1), 12)
        fixed, per_blob = work.blob_cluster_sets()
        for one in per_blob:
            for cluster_set in one:
                assert fixed <= work.softwired_clusters()
        assert fixed <= work.softwired_clusters()
        assert work.taxa() in fixed  # the root's cluster is always present


class TestDisplaysTrees:
    @pytest.mark.parametrize("seed", range(40))
    def test_agrees_with_displayed_trees(self, seed):
        work = sample(seed)
        if work.reticulation_number() > REFERENCE_CAP:
            pytest.skip("reference would be too slow")
        network = work.to_phylozoo()

        # ask about the network's own displayed trees, and about noise
        from phylozoo.core.network.dnetwork.derivations import displayed_trees

        groups = [hardwired_clusters(t) for t in displayed_trees(network)][:3]
        rng = random.Random(seed)
        taxa = sorted(work.taxa())
        if len(taxa) > 3:
            groups.append({frozenset(rng.sample(taxa, 2)), frozenset(rng.sample(taxa, 3))})

        for group in groups:
            assert displays_trees(network, [group]) == displays_trees_via_displayed_trees(
                network, [group]
            )

    def test_two_clusters_needing_different_switchings(self):
        g = WorkGraph()
        root, u, v, r = g.new_node(), g.new_node(), g.new_node(), g.new_node()
        g.add_edge(root, u)
        g.add_edge(root, v)
        g.add_edge(u, g.new_node(block=fs("a")))
        g.add_edge(v, g.new_node(block=fs("b")))
        g.add_edge(u, r)
        g.add_edge(v, r)
        g.add_edge(r, g.new_node(block=fs("c")))
        assert g.represents([fs("ac"), fs("bc")])
        assert not g.displays([[fs("ac"), fs("bc")]])
        assert g.displays([[fs("ac")], [fs("bc")]])

    def test_independent_blobs_may_use_different_switchings(self):
        """Two blobs are satisfied independently, which is the whole point."""
        from phylocass import CassOptions, cass_from_trees, read_trees

        trees = read_trees(
            "(((a,b),(c,d)),((e,f),(g,h))); (((a,c),(b,d)),((e,g),(f,h)));"
        )
        r = cass_from_trees(trees, CassOptions(display_trees=True))
        assert r.displays_input_trees() is True
        assert r.reticulation_number == 4


class TestScale:
    def test_a_network_far_past_enumeration(self):
        """2 ** 40-ish switchings, answered in milliseconds."""
        work = None
        for seed in range(200):
            candidate = random_blob_network(random.Random(seed), 120, max_level=2)
            if work is None or candidate.reticulation_number() > work.reticulation_number():
                work = candidate
            if work.reticulation_number() >= 30:
                break
        assert work.reticulation_number() >= 20, "expected a heavily reticulate sample"
        clusters = work.softwired_clusters()
        assert clusters
        assert work.taxa() in clusters
        network = work.to_phylozoo()
        assert softwired_clusters(network) == clusters
