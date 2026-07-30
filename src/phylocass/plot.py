"""Drawing a result: the input trees beside the network Cass built from them.

Uses PhyloZoo's matplotlib backend, which is an optional dependency::

    pip install "phylocass[viz]"

Nothing here is imported at package import time, so the core package keeps
working without matplotlib installed.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Sequence

__all__ = ["PlotUnavailable", "plot_result", "save_plot"]


class PlotUnavailable(RuntimeError):
    """Raised when the optional plotting dependencies are missing."""


def _matplotlib(headless: bool):
    try:
        import matplotlib
    except ImportError as exc:  # pragma: no cover - depends on the environment
        raise PlotUnavailable(
            "drawing needs matplotlib, which comes with PhyloZoo's viz extra: "
            'pip install "phylocass[viz]"'
        ) from exc
    if headless:
        # a file is all that is wanted, so do not go looking for a display
        matplotlib.use("Agg", force=False)
    import matplotlib.pyplot as plt

    return plt


def _plot():
    try:
        from phylozoo.viz import plot
    except ImportError as exc:  # pragma: no cover - depends on the environment
        raise PlotUnavailable(
            "drawing needs PhyloZoo's viz extra: pip install \"phylocass[viz]\""
        ) from exc
    return plot


def _tree_title(tree, index: int, all_taxa: frozenset) -> str:
    taxa = frozenset(tree.taxa)
    if taxa == all_taxa:
        return f"input tree {index}"
    missing = len(all_taxa) - len(taxa)
    return f"input tree {index}  ({len(taxa)} taxa, {missing} missing)"


def plot_result(
    result,
    trees: Sequence = (),
    layout: str = "auto",
    style=None,
    max_columns: int = 4,
    panel_size: float = 3.2,
    title: str | None = None,
    headless: bool = True,
):
    """Draw ``result``'s network, with each input tree beside it.

    Returns the matplotlib ``Figure``.  The input trees are laid out in a grid
    above, the network gets a full-width panel below, and reticulation edges
    are red (PhyloZoo's default).

    ``trees`` may be empty -- for a bare cluster set there is nothing to draw
    alongside, and only the network is plotted.
    """
    plt = _matplotlib(headless)
    plot = _plot()

    trees = list(trees)
    all_taxa = frozenset(result.taxa)

    if not trees:
        fig, ax = plt.subplots(figsize=(max(6.0, panel_size * 2), 5.0))
        plot(result.network, ax=ax, layout=layout, style=style)
        ax.set_axis_off()
        ax.set_title(_network_title(result), fontsize=11)
        if title:
            fig.suptitle(title, fontsize=13)
        fig.tight_layout()
        return fig

    columns = max(1, min(max_columns, len(trees)))
    rows = math.ceil(len(trees) / columns)

    fig = plt.figure(figsize=(panel_size * columns + 1.5, panel_size * rows + 4.2))
    grid = fig.add_gridspec(
        rows + 1,
        columns,
        height_ratios=[1] * rows + [1.6],
        hspace=0.28,
        wspace=0.08,
        top=0.9 if title else 0.94,
        bottom=0.04,
    )

    for index, tree in enumerate(trees):
        ax = fig.add_subplot(grid[index // columns, index % columns])
        plot(tree, ax=ax, layout=layout, style=style)
        ax.set_axis_off()
        ax.set_title(_tree_title(tree, index + 1, all_taxa), fontsize=10)

    ax = fig.add_subplot(grid[rows, :])
    plot(result.network, ax=ax, layout=layout, style=style)
    ax.set_axis_off()
    ax.set_title(_network_title(result), fontsize=11)

    if title:
        fig.suptitle(title, fontsize=13)
    return fig


def _network_title(result) -> str:
    plural = "" if result.reticulation_number == 1 else "s"
    headline = (
        f"Cass network — level {result.level}, "
        f"{result.reticulation_number} reticulation{plural}"
    )

    notes = []
    if result.display_trees:
        notes.append("displays the input trees")
    elif result.displays_input_trees() is False:
        notes.append("represents the clusters, but does not display the trees")
    if result.z_closure is not None and result.z_closure.dropped:
        notes.append(
            f"{result.z_closure.dropped} partial cluster"
            f"{'' if result.z_closure.dropped == 1 else 's'} dropped by Z-closure"
        )

    return headline + ("\n" + "; ".join(notes) if notes else "")


def save_plot(
    result,
    path: str | Path,
    trees: Sequence = (),
    dpi: int = 150,
    show: bool = False,
    **kwargs,
) -> Path:
    """Draw the result and write it to ``path``.

    The format follows the file extension, so ``.png``, ``.pdf`` and ``.svg``
    all work.  Returns the resolved path.
    """
    path = Path(path)
    fig = plot_result(result, trees, headless=not show, **kwargs)
    if path.parent and not path.parent.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    if show:  # pragma: no cover - interactive
        import matplotlib.pyplot as plt

        plt.show()
    else:
        import matplotlib.pyplot as plt

        plt.close(fig)
    return path.resolve()
