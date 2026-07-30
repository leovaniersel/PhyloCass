"""The command line, including how it behaves on bad input."""

import subprocess
import sys
from pathlib import Path

import pytest

from phylocass.cli import main

EXAMPLES = Path(__file__).resolve().parent.parent / "examples"


def run(args, stdin=None):
    # encoding is pinned so the BOM cases below survive a cp1252 locale
    return subprocess.run(
        [sys.executable, "-m", "phylocass", *args],
        input=stdin,
        capture_output=True,
        text=True,
        encoding="utf-8",
    )


class TestNormalUse:
    def test_newick_file(self):
        proc = run([str(EXAMPLES / "conflicting_trees.newick")])
        assert proc.returncode == 0
        assert proc.stdout.strip().endswith(";")
        assert "level=1" in proc.stderr

    def test_cluster_file(self):
        proc = run(["--format", "clusters", str(EXAMPLES / "figure1_clusters.txt")])
        assert proc.returncode == 0
        assert "level=2" in proc.stderr
        assert "reticulations=2" in proc.stderr

    def test_stdin(self):
        proc = run(["-q"], stdin="((a,b),(c,d)); ((a,c),(b,d));")
        assert proc.returncode == 0
        assert proc.stdout.count("#H") == 4  # two reticulations, each written twice

    def test_quiet_prints_only_the_network(self):
        proc = run(["-q", str(EXAMPLES / "conflicting_trees.newick")])
        assert proc.returncode == 0
        assert proc.stderr == ""
        assert len(proc.stdout.strip().splitlines()) == 1

    def test_show_clusters(self):
        proc = run(["--show-clusters", str(EXAMPLES / "conflicting_trees.newick")])
        assert proc.returncode == 0
        assert "clusters on" in proc.stderr

    def test_out_writes_a_file(self, tmp_path):
        target = tmp_path / "net.enewick"
        proc = run(["--out", str(target), str(EXAMPLES / "double_conflict.newick")])
        assert proc.returncode == 0
        assert target.read_text().strip().endswith(";")


class TestBadInput:
    def test_missing_file_is_reported_cleanly(self, tmp_path):
        proc = run([str(tmp_path / "nope.newick")])
        assert proc.returncode == 2
        assert "no such file" in proc.stderr
        assert "Traceback" not in proc.stderr

    def test_malformed_newick_is_reported_cleanly(self, tmp_path):
        bad = tmp_path / "bad.newick"
        bad.write_text("((a,b,c);\n", encoding="utf-8")
        proc = run([str(bad)])
        assert proc.returncode == 2
        assert proc.stderr.startswith("error:")
        assert "Traceback" not in proc.stderr

    def test_empty_input_is_reported(self, tmp_path):
        empty = tmp_path / "empty.newick"
        empty.write_text("# nothing here\n", encoding="utf-8")
        proc = run([str(empty)])
        assert proc.returncode == 2
        assert "no trees found" in proc.stderr

    def test_impossible_level_is_reported_cleanly(self):
        proc = run(["--max-level", "0", str(EXAMPLES / "double_conflict.newick")])
        assert proc.returncode == 1
        assert "gave up" in proc.stderr
        assert "Traceback" not in proc.stderr


class TestDisplayTrees:
    THREE = "(((a,b),c),(d,e)); (((a,c),b),(d,e)); (((a,d),b),(c,e));"

    def test_flag_changes_the_answer(self):
        plain = run(["-"], stdin=self.THREE)
        display = run(["--display-trees", "-"], stdin=self.THREE)
        assert plain.returncode == 0 and display.returncode == 0
        assert "displays-trees=False" in plain.stderr
        assert "displays-trees=True" in display.stderr

    def test_summary_reports_the_display_check(self):
        proc = run([str(EXAMPLES / "conflicting_trees.newick")])
        assert "displays-trees=" in proc.stderr

    def test_rejected_for_cluster_input(self):
        proc = run(
            ["--display-trees", "--format", "clusters", str(EXAMPLES / "figure1_clusters.txt")]
        )
        assert proc.returncode == 2
        assert "needs tree input" in proc.stderr
        assert "Traceback" not in proc.stderr


class TestBomOnStdin:
    """PowerShell 5.1 prepends a UTF-8 BOM when piping into a program."""

    def test_bom_prefixed_stdin(self):
        proc = run(["-q"], stdin="﻿((a,b),(c,d)); ((a,c),(b,d));")
        assert proc.returncode == 0, proc.stderr
        assert proc.stdout.strip().endswith(";")

    def test_bom_prefixed_comment_line_on_stdin(self):
        proc = run(["-q"], stdin="﻿# a comment\n((a,b),c);\n((a,c),b);\n")
        assert proc.returncode == 0, proc.stderr
        assert proc.stdout.strip().endswith(";")
