"""End-to-end tests for the Cass algorithm."""

import random

import pytest
from phylozoo.core.network.dnetwork.classifications import level as pz_level

from phylocass import CassOptions, cass, cass_from_trees, read_trees
from phylocass.io import clusters_of_trees, hardwired_clusters, softwired_clusters

from conftest import fs, nontrivial, random_network, random_tree, shown


# ----------------------------------------------------------------------
# hand-picked inputs
# ----------------------------------------------------------------------


class TestTrivialInputs:
    def test_identical_trees_give_a_tree(self):
        r = cass_from_trees(read_trees("((a,b),(c,d)); ((a,b),(c,d));"))
        assert r.level == 0
        assert r.reticulation_number == 0
        assert r.network.is_tree()
        assert r.represents_input()

    def test_single_tree_reproduces_itself(self):
        r = cass_from_trees(read_trees("(((a,b),c),(d,e));"))
        assert r.level == 0
        assert nontrivial(hardwired_clusters(r.network), r.taxa) == r.clusters

    def test_nested_clusters_are_not_a_conflict(self):
        r = cass(({fs("ab"), fs("abc")}), fs("abcd"))
        assert r.level == 0
        assert r.represents_input()


class TestLevelOne:
    def test_single_leaf_move(self):
        r = cass_from_trees(read_trees("(((a,b),c),d); (((a,c),b),d);"))
        assert r.level == 1
        assert r.reticulation_number == 1
        assert r.represents_input()

    def test_galled_tree_clusters(self):
        r = cass({fs("ab"), fs("bc"), fs("abc")}, fs("abcd"))
        assert r.level == 1
        assert r.represents_input()


class TestLevelTwo:
    def test_four_taxa_double_conflict_needs_level_two(self):
        """((a,b),(c,d)) versus ((a,c),(b,d)) cannot be one reticulation.

        A level-1 network displays two trees differing by relocating a single
        leaf; these two differ by swapping a pair, so level 2 is optimal.
        """
        r = cass_from_trees(read_trees("((a,b),(c,d)); ((a,c),(b,d));"))
        assert r.level == 2
        assert r.reticulation_number == 2
        assert r.represents_input()

    def test_paper_figure_one(self):
        """The nine-taxon example from Figure 1 of the paper.

        The paper reports that a level-2 network with two reticulations exists
        here and that Cass finds it, where the galled-network algorithm needs
        four reticulations.
        """
        clusters = {
            fs(s)
            for s in (
                "abfgi", "abcfgi", "abfi", "bcfi", "cdeh", "deh",
                "bcfhi", "bcdfhi", "bci", "ag", "bi", "ci", "dh",
            )
        }
        r = cass(clusters, fs("abcdefghi"))
        assert r.level == 2
        assert r.reticulation_number == 2
        assert r.represents_input()


class TestDecomposition:
    def test_independent_conflicts_stay_in_separate_blobs(self):
        r = cass_from_trees(
            read_trees("(((a,b),(c,d)),((e,f),(g,h))); (((a,c),(b,d)),((e,g),(f,h)));")
        )
        assert r.level == 2
        assert r.reticulation_number == 4  # two independent level-2 blobs
        assert r.represents_input()

    def test_conflict_plus_clean_subtree(self):
        r = cass_from_trees(read_trees("(((a,b),c),(x,y)); (((a,c),b),(x,y));"))
        assert r.level == 1
        assert r.represents_input()
        assert fs("xy") in softwired_clusters(r.network)

    def test_taxa_are_preserved_exactly_once(self):
        r = cass_from_trees(read_trees("(((a,b),c),d); (((a,c),b),d);"))
        assert frozenset(r.network.taxa) == fs("abcd")
        assert len(r.network.leaves) == 4


# ----------------------------------------------------------------------
# the output is a well-formed PhyloZoo network
# ----------------------------------------------------------------------


class TestPhyloZooOutput:
    @pytest.mark.parametrize(
        "text",
        [
            "((a,b),(c,d)); ((a,c),(b,d));",
            "(((a,b),c),d); (((a,c),b),d);",
            "(((a,b),(c,d)),((e,f),(g,h))); (((a,c),(b,d)),((e,g),(f,h)));",
        ],
    )
    def test_output_validates(self, text):
        r = cass_from_trees(read_trees(text))
        r.network.validate()

    def test_reported_level_is_phylozoos_level(self):
        r = cass_from_trees(read_trees("((a,b),(c,d)); ((a,c),(b,d));"))
        assert r.level == pz_level(r.network)

    def test_enewick_round_trips_through_phylozoo(self):
        from phylozoo import DirectedPhyNetwork

        r = cass_from_trees(read_trees("((a,b),(c,d)); ((a,c),(b,d));"))
        again = DirectedPhyNetwork.from_string(r.to_enewick())
        assert softwired_clusters(again) == softwired_clusters(r.network)

    def test_no_edge_joins_two_reticulations(self):
        r = cass_from_trees(read_trees("((a,b),(c,d)); ((a,c),(b,d));"))
        hybrids = set(r.network.hybrid_nodes)
        for u, v in r.network.edges:
            assert not (u in hybrids and v in hybrids)


# ----------------------------------------------------------------------
# round-trip: rediscover a network from the clusters it displays
# ----------------------------------------------------------------------


@pytest.mark.parametrize("seed", range(40))
def test_roundtrip_level_one(seed):
    """Cass must find a level-<=1 network for clusters that came from one."""
    rng = random.Random(seed)
    work = random_network(rng, rng.randint(3, 5), 1)
    if work.level() != 1:
        pytest.skip("random network degenerated")
    clusters = nontrivial(work.softwired_clusters(), work.taxa())
    if not clusters:
        pytest.skip("no non-trivial clusters")
    r = cass(clusters, work.taxa(), CassOptions(max_level=3))
    assert r.represents_input()
    assert r.level <= 1, f"level {r.level} for clusters from a level-1 network"


@pytest.mark.parametrize("seed", range(25))
def test_roundtrip_level_two(seed):
    """The paper's guarantee: exact whenever a level-<=2 network exists."""
    rng = random.Random(1000 + seed)
    work = random_network(rng, rng.randint(3, 4), 2)
    if work.level() != 2:
        pytest.skip("random network degenerated")
    clusters = nontrivial(work.softwired_clusters(), work.taxa())
    if not clusters:
        pytest.skip("no non-trivial clusters")
    r = cass(clusters, work.taxa(), CassOptions(max_level=4))
    assert r.represents_input()
    assert r.level <= 2, f"level {r.level} for clusters from a level-2 network"


@pytest.mark.parametrize("seed", range(30))
def test_random_trees_always_produce_a_valid_network(seed):
    """Whatever the input trees, the output must display every input cluster."""
    rng = random.Random(500 + seed)
    names = [chr(ord("a") + i) for i in range(rng.randint(4, 6))]
    trees = [random_tree(rng, names).to_phylozoo() for _ in range(rng.randint(2, 3))]
    clusters, taxa = clusters_of_trees(trees)
    try:
        r = cass(clusters, taxa, CassOptions(max_level=4, time_limit=5.0))
    except RuntimeError:
        pytest.skip("search budget exhausted")
    assert r.represents_input()
    r.network.validate()
    assert frozenset(r.network.taxa) == taxa


# ----------------------------------------------------------------------
# budgets
# ----------------------------------------------------------------------


class TestDeterminism:
    """Several valid networks usually exist; the choice must not depend on hashing."""

    @pytest.mark.parametrize(
        "text",
        [
            "(((a,b),c),d); (((a,c),b),d);",
            "((a,b),(c,d)); ((a,c),(b,d));",
            "(((a,b),(c,d)),((e,f),(g,h))); (((a,c),(b,d)),((e,g),(f,h)));",
        ],
    )
    def test_repeated_runs_agree(self, text):
        answers = {cass_from_trees(read_trees(text)).to_enewick() for _ in range(5)}
        assert len(answers) == 1

    def test_shuffling_the_input_clusters_changes_nothing(self):
        clusters = [fs(s) for s in ("abfgi", "abfi", "bcfi", "bci", "ag", "bi", "ci")]
        taxa = fs("abcfgi")
        rng = random.Random(0)
        answers = set()
        for _ in range(5):
            shuffled = clusters[:]
            rng.shuffle(shuffled)
            answers.add(cass(set(shuffled), taxa).to_enewick())
        assert len(answers) == 1

    def test_subprocesses_with_different_hash_seeds_agree(self):
        """The real check: PYTHONHASHSEED only varies across processes."""
        import subprocess
        import sys

        code = (
            "from phylocass import read_trees, cass_from_trees;"
            "print(cass_from_trees(read_trees('((a,b),(c,d)); ((a,c),(b,d));'))"
            ".to_enewick())"
        )
        outputs = set()
        for seed in ("1", "2", "3"):
            out = subprocess.run(
                [sys.executable, "-c", code],
                capture_output=True,
                text=True,
                env={"PYTHONHASHSEED": seed, "PATH": ""},
                check=True,
            )
            outputs.add(out.stdout.strip())
        assert len(outputs) == 1


class TestBudgets:
    def test_max_level_zero_conflicts_are_reported(self):
        with pytest.raises(RuntimeError, match="gave up"):
            cass({fs("ab"), fs("bc")}, fs("abc"), CassOptions(max_level=0))

    def test_giveup_message_names_the_knobs(self):
        with pytest.raises(RuntimeError) as excinfo:
            cass({fs("ab"), fs("bc")}, fs("abc"), CassOptions(max_level=0))
        assert "max_level" in str(excinfo.value)

    def test_time_limit_is_one_budget_for_the_whole_component(self):
        """The budget must not be re-granted at every level the search climbs."""
        import time as _time

        clusters = {
            fs(s) for s in ("ab", "bc", "cd", "de", "ea", "ac", "bd", "ce", "da", "eb")
        }
        t0 = _time.monotonic()
        try:
            cass(clusters, fs("abcdef"), CassOptions(max_level=6, time_limit=1.0))
        except RuntimeError:
            pass
        elapsed = _time.monotonic() - t0
        assert elapsed < 4.0, f"took {elapsed:.1f}s; budget looks per-level"
