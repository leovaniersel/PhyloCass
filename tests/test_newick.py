import pytest

from phylocass.newick import NewickError, parse_newick, read_trees


def fs(s):
    return frozenset(s)


def shown(sets):
    return sorted(tuple(sorted(s)) for s in sets)


class TestParsing:
    def test_simple_tree(self):
        t = parse_newick("((a,b),(c,d));")
        assert shown(t.taxon_labels()) == [("a",), ("b",), ("c",), ("d",)]
        assert fs("ab") in t.clusters()
        assert fs("cd") in t.clusters()

    def test_branch_lengths_are_ignored(self):
        t = parse_newick("((a:0.1,b:0.2):0.5,c:1.0);")
        assert shown(t.taxon_labels()) == [("a",), ("b",), ("c",)]

    def test_internal_labels_are_ignored(self):
        t = parse_newick("((a,b)Node1:0.4,c)Root;")
        assert shown(t.taxon_labels()) == [("a",), ("b",), ("c",)]
        assert fs("ab") in t.clusters()

    def test_comments_are_skipped(self):
        t = parse_newick("((a,b)[an internal comment],c);")
        assert shown(t.taxon_labels()) == [("a",), ("b",), ("c",)]

    def test_quoted_labels(self):
        t = parse_newick("(('Homo sapiens','Pan troglodytes'),Gorilla);")
        assert ("Homo sapiens",) in shown(t.taxon_labels())

    def test_underscores_become_spaces(self):
        t = parse_newick("((Homo_sapiens,Pan),Gorilla);")
        assert ("Homo sapiens",) in shown(t.taxon_labels())

    def test_missing_semicolon_is_tolerated(self):
        assert len(parse_newick("(a,b)").leaves()) == 2

    def test_duplicate_labels_rejected(self):
        with pytest.raises(NewickError):
            parse_newick("((a,b),a);")

    def test_unbalanced_parentheses_rejected(self):
        with pytest.raises(NewickError):
            parse_newick("((a,b,c);")

    def test_trailing_junk_rejected(self):
        with pytest.raises(NewickError):
            parse_newick("(a,b); (c,d)")


class TestReadTrees:
    def test_multiple_trees_on_one_line(self):
        trees = read_trees("((a,b),c); ((a,c),b);")
        assert len(trees) == 2

    def test_multiple_trees_on_separate_lines(self):
        trees = read_trees("((a,b),c);\n((a,c),b);\n")
        assert len(trees) == 2

    def test_comment_lines_and_blanks_skipped(self):
        trees = read_trees("# gene 1\n((a,b),c);\n\n# gene 2\n((a,c),b);\n")
        assert len(trees) == 2

    def test_tree_split_across_lines(self):
        trees = read_trees("((a,b),\n(c,d));\n")
        assert len(trees) == 1
        assert len(trees[0].leaves()) == 4
