"""Drawing the result. Skipped when the optional viz extra is absent."""

import subprocess
import sys
from pathlib import Path

import pytest

from phylocass import CassOptions, cass, cass_from_trees, read_trees

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from phylocass.plot import plot_result, save_plot  # noqa: E402

from conftest import fs  # noqa: E402

EXAMPLES = Path(__file__).resolve().parent.parent / "examples"


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


def run_cli(args):
    return subprocess.run(
        [sys.executable, "-m", "phylocass", *args],
        capture_output=True,
        text=True,
        encoding="utf-8",
    )


class TestPlotResult:
    def test_returns_a_figure_with_a_panel_per_tree_plus_the_network(self):
        trees = read_trees("(((a,b),c),d); (((a,c),b),d);")
        fig = plot_result(cass_from_trees(trees), trees)
        assert len(fig.axes) == len(trees) + 1

    def test_network_only_when_there_are_no_trees(self):
        result = cass({fs("ab"), fs("bc")}, fs("abcd"))
        fig = plot_result(result)
        assert len(fig.axes) == 1

    def test_many_trees_wrap_into_a_grid(self):
        trees = read_trees(
            "((a,b),(c,d)); ((a,c),(b,d)); ((a,d),(b,c)); "
            "(((a,b),c),d); (((a,c),b),d);"
        )
        result = cass_from_trees(trees, CassOptions(max_level=5))
        fig = plot_result(result, trees, max_columns=3)
        assert len(fig.axes) == 6

    def test_title_mentions_level_and_reticulations(self):
        trees = read_trees("(((a,b),c),d); (((a,c),b),d);")
        fig = plot_result(cass_from_trees(trees), trees)
        title = fig.axes[-1].get_title()
        assert "level 1" in title
        assert "1 reticulation" in title

    def test_title_flags_when_the_trees_are_not_displayed(self):
        trees = read_trees("(((a,b),c),(d,e)); (((a,c),b),(d,e)); (((a,d),b),(c,e));")
        fig = plot_result(cass_from_trees(trees), trees)
        assert "does not display the trees" in fig.axes[-1].get_title()

    def test_partial_trees_are_labelled_with_their_taxon_count(self):
        trees = read_trees("(((a,b),c),d); (((a,c),b),e);")
        result = cass_from_trees(trees, CassOptions(z_closure=True))
        fig = plot_result(result, trees)
        assert "missing" in fig.axes[0].get_title()


class TestSavePlot:
    @pytest.mark.parametrize("suffix", [".png", ".pdf", ".svg"])
    def test_writes_the_requested_format(self, tmp_path, suffix):
        trees = read_trees("(((a,b),c),d); (((a,c),b),d);")
        out = save_plot(cass_from_trees(trees), tmp_path / f"net{suffix}", trees=trees)
        assert out.exists() and out.stat().st_size > 0

    def test_creates_missing_directories(self, tmp_path):
        trees = read_trees("((a,b),c);")
        out = save_plot(cass_from_trees(trees), tmp_path / "deep" / "net.png", trees=trees)
        assert out.exists()

    def test_returns_a_resolved_path(self, tmp_path):
        trees = read_trees("((a,b),c);")
        out = save_plot(cass_from_trees(trees), tmp_path / "net.png", trees=trees)
        assert out.is_absolute()


class TestCli:
    def test_plot_flag_writes_an_image(self, tmp_path):
        out = tmp_path / "net.png"
        proc = run_cli([str(EXAMPLES / "conflicting_trees.newick"), "--plot", str(out)])
        assert proc.returncode == 0, proc.stderr
        assert out.stat().st_size > 0
        assert "wrote" in proc.stderr

    def test_plot_still_prints_the_network(self, tmp_path):
        out = tmp_path / "net.png"
        proc = run_cli([str(EXAMPLES / "conflicting_trees.newick"), "--plot", str(out)])
        assert proc.stdout.strip().endswith(";")

    def test_quiet_suppresses_the_wrote_line(self, tmp_path):
        out = tmp_path / "net.png"
        proc = run_cli([str(EXAMPLES / "conflicting_trees.newick"), "--plot", str(out), "-q"])
        assert proc.returncode == 0
        assert proc.stderr == ""
        assert out.exists()

    def test_plot_works_for_cluster_input(self, tmp_path):
        out = tmp_path / "net.png"
        proc = run_cli(
            ["--format", "clusters", str(EXAMPLES / "figure1_clusters.txt"),
             "--plot", str(out)]
        )
        assert proc.returncode == 0, proc.stderr
        assert out.stat().st_size > 0

    def test_plot_works_with_z_closure(self, tmp_path):
        out = tmp_path / "net.png"
        proc = run_cli(
            [str(EXAMPLES / "partial_trees.newick"), "--z-closure", "--plot", str(out)]
        )
        assert proc.returncode == 0, proc.stderr
        assert out.stat().st_size > 0

    def test_plot_dpi_changes_the_file(self, tmp_path):
        small, large = tmp_path / "s.png", tmp_path / "l.png"
        for target, dpi in ((small, "50"), (large, "200")):
            run_cli(
                [str(EXAMPLES / "conflicting_trees.newick"),
                 "--plot", str(target), "--plot-dpi", dpi]
            )
        assert large.stat().st_size > small.stat().st_size


class TestExampleScript:
    def test_script_runs(self, tmp_path):
        out = tmp_path / "net.png"
        proc = subprocess.run(
            [sys.executable, str(EXAMPLES / "plot_network.py"), "-o", str(out)],
            capture_output=True,
            text=True,
        )
        assert proc.returncode == 0, proc.stderr
        assert out.stat().st_size > 0
