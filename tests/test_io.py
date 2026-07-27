"""Reading trees and extracting clusters (PhyloZoo does the parsing)."""

import pytest
from phylozoo import DirectedPhyNetwork

from phylocass.io import (
    clusters_of_trees,
    hardwired_clusters,
    read_cluster_file,
    read_tree_file,
    read_trees,
    softwired_clusters,
)

from conftest import fs, shown


class TestReadTrees:
    def test_single_tree(self):
        trees = read_trees("((a,b),(c,d));")
        assert len(trees) == 1
        assert sorted(trees[0].taxa) == ["a", "b", "c", "d"]

    def test_multiple_trees_on_one_line(self):
        assert len(read_trees("((a,b),c); ((a,c),b);")) == 2

    def test_multiple_trees_on_separate_lines(self):
        assert len(read_trees("((a,b),c);\n((a,c),b);\n")) == 2

    def test_comment_and_blank_lines_skipped(self):
        text = "# gene 1\n((a,b),c);\n\n# gene 2\n((a,c),b);\n"
        assert len(read_trees(text)) == 2

    def test_tree_split_across_lines(self):
        trees = read_trees("((a,b),\n(c,d));\n")
        assert len(trees) == 1
        assert len(trees[0].taxa) == 4

    def test_branch_lengths_and_internal_labels(self):
        trees = read_trees("((a:0.1,b:0.2)Node:0.3,c);")
        assert sorted(trees[0].taxa) == ["a", "b", "c"]

    def test_quoted_labels(self):
        trees = read_trees("(('Homo sapiens',Pan),Gorilla);")
        assert "Homo sapiens" in trees[0].taxa

    def test_malformed_input_raises(self):
        with pytest.raises(Exception):
            read_trees("((a,b,c);")

    def test_empty_input_gives_no_trees(self):
        assert read_trees("\n# nothing here\n\n") == []


class TestClusterExtraction:
    def test_hardwired_clusters_of_a_tree(self):
        tree = DirectedPhyNetwork.from_string("((a,b),(c,d));")
        assert shown(hardwired_clusters(tree)) == [
            ("a",), ("a", "b"), ("a", "b", "c", "d"), ("b",), ("c",), ("c", "d"), ("d",)
        ]

    def test_softwired_clusters_of_a_network(self):
        net = DirectedPhyNetwork.from_string("(((a,(b)#H1),c),(d,#H1));")
        soft = softwired_clusters(net)
        # b can sit with a or with d
        assert fs("ab") in soft
        assert fs("bd") in soft

    def test_softwired_superset_of_hardwired(self):
        net = DirectedPhyNetwork.from_string("(((a,(b)#H1),c),(d,#H1));")
        assert hardwired_clusters(net) <= softwired_clusters(net)

    def test_clusters_of_trees_drops_trivial_ones(self):
        trees = read_trees("((a,b),(c,d)); ((a,c),(b,d));")
        clusters, taxa = clusters_of_trees(trees)
        assert taxa == fs("abcd")
        assert shown(clusters) == [("a", "b"), ("a", "c"), ("b", "d"), ("c", "d")]
        assert all(1 < len(c) < 4 for c in clusters)

    def test_clusters_of_trees_unions_over_trees(self):
        trees = read_trees("(((a,b),c),d); ((a,b),(c,d));")
        clusters, _ = clusters_of_trees(trees)
        assert fs("ab") in clusters
        assert fs("abc") in clusters
        assert fs("cd") in clusters


class TestByteOrderMarks:
    """Windows sprinkles BOMs: Notepad writes one, PowerShell pipes one in.

    Left in place a BOM parses as a taxon, and it also stops a leading ``#``
    comment line from being recognised as one.
    """

    BOM = "﻿"

    def test_bom_at_the_start_of_a_string(self):
        trees = read_trees(self.BOM + "((a,b),c);")
        assert sorted(trees[0].taxa) == ["a", "b", "c"]

    def test_bom_before_a_comment_line(self):
        text = self.BOM + "# a comment\n((a,b),c);\n((a,c),b);\n"
        assert len(read_trees(text)) == 2

    def test_bom_written_by_notepad_style_utf8_sig(self, tmp_path):
        p = tmp_path / "t.newick"
        p.write_text("# trees\n((a,b),c);\n((a,c),b);\n", encoding="utf-8-sig")
        assert p.read_bytes().startswith(b"\xef\xbb\xbf")
        assert len(read_tree_file(p)) == 2

    def test_bom_in_a_cluster_file(self, tmp_path):
        p = tmp_path / "c.txt"
        p.write_text("# clusters\na b\nc d\n", encoding="utf-8-sig")
        clusters, taxa = read_cluster_file(p)
        assert shown(clusters) == [("a", "b"), ("c", "d")]
        assert taxa == fs("abcd")


class TestFiles:
    def test_read_tree_file(self, tmp_path):
        p = tmp_path / "t.newick"
        p.write_text("# two trees\n((a,b),c);\n((a,c),b);\n", encoding="utf-8")
        assert len(read_tree_file(p)) == 2

    def test_read_cluster_file(self, tmp_path):
        p = tmp_path / "c.txt"
        p.write_text("# clusters\na b\nc, d\n\n", encoding="utf-8")
        clusters, taxa = read_cluster_file(p)
        assert shown(clusters) == [("a", "b"), ("c", "d")]
        assert taxa == fs("abcd")
