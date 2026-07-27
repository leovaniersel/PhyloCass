"""The documented example must keep working.

Plotting is optional, so these skip when PhyloZoo's viz extra is absent.
"""

import subprocess
import sys
from pathlib import Path

import pytest

from phylocass import cass_from_trees, read_trees
from phylocass.io import hardwired_clusters

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")

EXAMPLES = Path(__file__).resolve().parent.parent / "examples"

README_TREES = """
    (((a,b),c),(d,e));
    (((a,c),b),(d,e));
"""


class TestReadmeExample:
    def test_numbers_quoted_in_the_readme(self):
        result = cass_from_trees(read_trees(README_TREES))
        assert result.level == 1
        assert result.reticulation_number == 1
        assert result.to_enewick() == "((d,e),((b,(a)#H1),(#H1,c)));"
        assert result.represents_input()

    def test_the_two_switchings_are_the_two_input_trees(self):
        from phylozoo.core.network.dnetwork.derivations import displayed_trees

        trees = read_trees(README_TREES)
        result = cass_from_trees(trees)
        displayed = [hardwired_clusters(t) for t in displayed_trees(result.network)]
        for tree in trees:
            assert hardwired_clusters(tree) in displayed
        assert len(displayed) == 2

    def test_plot_returns_axes(self):
        from phylozoo.viz import plot

        result = cass_from_trees(read_trees(README_TREES))
        ax = plot(result.network)
        ax.set_title("smoke")
        assert ax.figure is not None
        matplotlib.pyplot.close("all")

    def test_style_and_shared_axes(self):
        from phylozoo.viz import plot
        from phylozoo.viz.dnetwork import DNetStyle

        trees = read_trees(README_TREES)
        result = cass_from_trees(trees)
        fig, (left, right) = matplotlib.pyplot.subplots(1, 2)
        plot(trees[0], ax=left)
        plot(
            result.network,
            ax=right,
            style=DNetStyle(hybrid_edge_color="crimson", node_size=300),
        )
        matplotlib.pyplot.close("all")

    @pytest.mark.parametrize("layout", ["spring", "circular"])
    def test_networkx_layouts_work_without_graphviz(self, layout):
        from phylozoo.viz import plot

        result = cass_from_trees(read_trees(README_TREES))
        plot(result.network, layout=layout)
        matplotlib.pyplot.close("all")

    def test_saving_to_enewick_and_dot(self, tmp_path):
        result = cass_from_trees(read_trees(README_TREES))
        for name in ("network.enewick", "network.dot"):
            target = tmp_path / name
            result.network.save(target)
            assert target.stat().st_size > 0


class TestPlotScript:
    def test_runs_and_writes_an_image(self, tmp_path):
        out = tmp_path / "network.png"
        proc = subprocess.run(
            [sys.executable, str(EXAMPLES / "plot_network.py"), "-o", str(out)],
            capture_output=True,
            text=True,
        )
        assert proc.returncode == 0, proc.stderr
        assert out.stat().st_size > 0
        assert "level              1" in proc.stdout

    def test_accepts_an_input_file(self, tmp_path):
        out = tmp_path / "network.png"
        proc = subprocess.run(
            [
                sys.executable,
                str(EXAMPLES / "plot_network.py"),
                str(EXAMPLES / "double_conflict.newick"),
                "-o",
                str(out),
            ],
            capture_output=True,
            text=True,
        )
        assert proc.returncode == 0, proc.stderr
        assert out.stat().st_size > 0
        assert "level              2" in proc.stdout

    def test_committed_figure_is_present(self):
        figure = EXAMPLES.parent / "docs" / "example-network.png"
        assert figure.exists(), "README references docs/example-network.png"
        assert figure.stat().st_size > 0
