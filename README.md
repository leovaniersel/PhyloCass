# PhyloCass

[![PyPI](https://img.shields.io/pypi/v/phylocass.svg)](https://pypi.org/project/phylocass/)
[![Python versions](https://img.shields.io/pypi/pyversions/phylocass.svg)](https://pypi.org/project/phylocass/)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

The **Cass** algorithm, implemented on top of [PhyloZoo](https://github.com/nholtgrefe/phylozoo).

> L. van Iersel, S. Kelk, R. Rupp, D. Huson.
> *Phylogenetic Networks Do not Need to Be Complex: Using Fewer Reticulations to Represent Conflicting Clusters.*
> Bioinformatics **26**(12):i124–i131, 2010 ([arXiv:0910.3082](https://arxiv.org/abs/0910.3082))

Cass takes a set of rooted phylogenetic trees, collects every cluster those
trees display, and builds a single rooted phylogenetic network that represents
**all** of those clusters in the softwired sense — while keeping the *level* of
the network as low as it can.

The level is the largest number of reticulations inside any single biconnected
component. Minimising it, rather than the total reticulation number, is what
lets Cass produce networks that stay readable even when the input trees
disagree a lot.

Results come back as PhyloZoo `DirectedPhyNetwork` objects, so they drop
straight into the rest of that ecosystem — plotting, format conversion,
displayed trees, quartets, and so on.

## Install

```bash
pip install phylocass            # core
pip install "phylocass[viz]"     # with plotting
```

Pulls in `phylozoo>=0.2.6`. Python 3.10+. Installing puts a `phylocass`
command on your path; `python -m phylocass` does the same thing.

From a checkout instead, for development:

```bash
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\Activate.ps1
pip install -e ".[dev]"
```

<details>
<summary>Windows: "running scripts is disabled on this system"</summary>

PowerShell refuses to run `Activate.ps1` when the execution policy is
`Restricted`, which is the default on Windows client editions. Check with
`Get-ExecutionPolicy -List` — if every scope says `Undefined`, that default is
what applies.

You do not have to activate anything; the launchers in the environment work
directly, and need no policy change at all:

```powershell
.\.venv\Scripts\phylocass.exe examples\conflicting_trees.newick
.\.venv\Scripts\python.exe -m pytest
```

To activate anyway, allow local scripts for this window:

```powershell
Set-ExecutionPolicy -Scope Process RemoteSigned
```

or once and for all for your user (no administrator rights needed):

```powershell
Set-ExecutionPolicy -Scope CurrentUser RemoteSigned
```

`RemoteSigned` is enough — it only blocks *downloaded* unsigned scripts, and
the one in your environment is local. `cmd.exe` is unaffected either way:
`.venv\Scripts\activate.bat`.

</details>

## Command line

```bash
# from a file of Newick trees, one tree per ';'
phylocass examples/conflicting_trees.newick

# from an explicit cluster list, one cluster per line
phylocass --format clusters examples/figure1_clusters.txt

# from stdin
echo '((a,b),(c,d)); ((a,c),(b,d));' | phylocass

# write the network out through PhyloZoo (.enewick/.nwk, or .dot)
phylocass examples/double_conflict.newick --out network.enewick
```

The network goes to stdout and a summary to stderr, so
`phylocass in.newick > out.enewick` gives you just the network:

```
# taxa=4 clusters=4 level=2 reticulations=2 verified=True
(((d,(b)#H2),((a)#H1,c)),(#H1,#H2));
```

Output is deterministic: the same input always gives the same network, whatever
the cluster ordering or Python's hash seed.

Exit codes are `0` on success, `2` for unusable input (missing file, malformed
Newick, nothing to do) and `1` when the search gives up under `--max-level` or
`--time-limit`. Input is always read as UTF-8, and a byte-order mark is
ignored, so piping from PowerShell and files saved by Notepad both work.

## Library

```python
from phylocass import read_trees, cass_from_trees

trees = read_trees("""
    (((a,b),c),d);
    (((a,c),b),d);
""")

result = cass_from_trees(trees)
result.level                 # 1
result.reticulation_number   # 1
result.to_enewick()          # '(d,((b,(a)#H1),(c,#H1)));'
result.represents_input()    # True

result.network               # a phylozoo DirectedPhyNetwork
result.network.validate()
sorted(result.network.taxa)  # ['a', 'b', 'c', 'd']
```

Or start from clusters directly:

```python
from phylocass import cass

result = cass({frozenset("ab"), frozenset("bc")}, taxa=frozenset("abcd"))
```

Because the result is an ordinary PhyloZoo network, everything PhyloZoo offers
applies to it:

```python
from phylozoo.core.network.dnetwork.derivations import displayed_trees
from phylozoo.core.network.dnetwork.classifications import level

[t.to_string() for t in displayed_trees(result.network)]
level(result.network)
result.network.save("out.dot")
```

`result.represents_input()` re-checks the finished network against the input
clusters through PhyloZoo's `displayed_trees` — a different code path from the
one the search uses, so it is a genuine verification rather than a restatement.

## Drawing it

`--plot` writes the input trees and the resulting network to one image:

```bash
phylocass examples/three_trees.newick --plot network.png
phylocass examples/three_trees.newick --plot network.pdf --plot-dpi 300
phylocass my_trees.newick --show                    # open a window instead
```

The format follows the extension — `.png`, `.pdf` and `.svg` all work — and
`--show` opens a window (on its own, or alongside `--plot`). Reticulation edges
are red. The network's caption reports the level and reticulation count, and
flags the things worth noticing: whether the trees themselves are displayed,
and how many partial clusters Z-closure had to drop.

Needs PhyloZoo's plotting extra:

```bash
pip install "phylocass[viz]"
```

From Python, the same thing:

```python
from phylocass import cass_from_trees, read_tree_file
from phylocass.plot import plot_result, save_plot

trees = read_tree_file("examples/three_trees.newick")
result = cass_from_trees(trees)

save_plot(result, "network.png", trees=trees, dpi=300)

fig = plot_result(result, trees)     # or take the matplotlib Figure and adjust
fig.axes[-1].set_title("my caption")
```

`plot_result` accepts `layout=` and `style=` and passes them to PhyloZoo, lays
the trees out in a grid (`max_columns`), and works with no trees at all — for a
bare cluster set it just draws the network.

## Worked example: trees in, network out, drawn

Two gene trees that agree that `d` and `e` are sisters but disagree about where
`a` belongs — one groups it with `b`, the other with `c`. No single tree
displays both `{a,b}` and `{a,c}`, so Cass resolves the conflict with one
reticulation.

```python
import matplotlib.pyplot as plt
from phylozoo.viz import plot

from phylocass import cass_from_trees, read_trees

trees = read_trees("""
    (((a,b),c),(d,e));
    (((a,c),b),(d,e));
""")

result = cass_from_trees(trees)

print(result.level)                # 1
print(result.reticulation_number)  # 1
print(result.to_enewick())         # ((d,e),((b,(a)#H1),(#H1,c)));
print(result.represents_input())   # True

# result.network is a PhyloZoo DirectedPhyNetwork, so PhyloZoo draws it
ax = plot(result.network)
ax.set_title(f"level {result.level}, {result.reticulation_number} reticulation")
plt.savefig("network.png", dpi=150, bbox_inches="tight")
plt.show()
```

That draws the network on its own. To get the input trees beside it, as below,
use `--plot` (or `phylocass.plot.save_plot`):

```bash
phylocass examples/conflicting_trees.newick --plot network.png
```

[`examples/plot_network.py`](examples/plot_network.py) is the same thing as a
script you can edit.

![Two conflicting input trees and the level-1 network Cass builds from their clusters](docs/example-network.png)

Reticulation edges are red by default. The one reticulation has two incoming
edges, and switching on one or the other recovers exactly the two input trees:

```python
from phylozoo.core.network.dnetwork.derivations import displayed_trees

[t.to_string() for t in displayed_trees(result.network)]
# ['((d,e),(b,(a,c)));', '((d,e),(c,(a,b)));']
```

Those are the input trees written with a different child order — the same two
cluster sets. Nothing was lost and nothing spurious was added.

`plot` takes a matplotlib `ax`, a layout name and a style object, so it composes
normally:

```python
from phylozoo.viz.dnetwork import DNetStyle

fig, (left, right) = plt.subplots(1, 2, figsize=(11, 4))
plot(trees[0], ax=left)
plot(result.network, ax=right, style=DNetStyle(hybrid_edge_color="crimson", node_size=300))
for ax in (left, right):
    ax.set_axis_off()
```

Layouts other than the default `'pz-dag'` are available: any NetworkX layout
(`'spring'`, `'circular'`, `'kamada_kawai'`, …) works out of the box, and the
Graphviz layouts (`'dot'`, `'neato'`, …) work if `pygraphviz` is installed —
without it, `plot(net, layout='dot')` raises `PhyloZooImportError`. For a
network with several tangled reticulations, `layout='dot'` is often easier to
read than the default.

If you would rather draw it elsewhere, save in a format another tool
understands:

```python
result.network.save("network.enewick")   # eNewick, e.g. for Dendroscope
result.network.save("network.dot")       # DOT, e.g. for Graphviz or Gephi
```

## Displaying the trees, not just their clusters

By default Cass does what the paper describes: every input cluster must be
represented by **some** switching of the network. Two clusters of the same
input tree are free to need *different* switchings — so the network represents
all the clusters without necessarily displaying any of the trees they came
from.

Setting `display_trees` strengthens the requirement: all clusters of one input
tree must appear in **one and the same** switching. That switching's tree then
displays the input tree, so the network displays the trees themselves. Its
reticulation number is therefore an upper bound on the **hybridization number**
of the input trees, which turns Cass into a heuristic for Hybridization Number
on multiple trees.

```bash
phylocass examples/three_trees.newick                  # level 2, 2 reticulations
phylocass examples/three_trees.newick --display-trees  # level 3, 3 reticulations
```

```python
from phylocass import CassOptions, cass_from_trees, read_tree_file

trees = read_tree_file("examples/three_trees.newick")

plain = cass_from_trees(trees)
plain.reticulation_number       # 2
plain.represents_input()        # True  - all clusters are there
plain.displays_input_trees()    # False - but not one switching per tree

hybrid = cass_from_trees(trees, CassOptions(display_trees=True))
hybrid.reticulation_number      # 3  <- heuristic bound on the hybridization number
hybrid.displays_input_trees()   # True
```

Those three trees are the smallest case in `examples/` where the modes differ:
they agree that `d` and `e` are sisters but each place `a` somewhere else.

`displays_input_trees()` is checked through PhyloZoo's `displayed_trees`, so it
verifies the finished network independently of the search — and it is worth
calling in either mode, since cluster mode often happens to display the trees
anyway. It returns `None` when Cass was handed a bare cluster set, because
then nothing records which cluster came from which tree.

### Two trees are a special case

For **two binary trees on the same taxon set** the minimum is the same whichever
model you use — displaying the trees, or representing their clusters in the
softwired sense — both for the reticulation number and for the level:

> L. van Iersel, S. Kelk. *When two trees go to war.*
> Journal of Theoretical Biology **269**(1):245–255, 2011
> ([arXiv:1004.5332](https://arxiv.org/abs/1004.5332)),
> Corollary 1 (reticulation number) and Theorem 2 (level).

Over 600 random pairs of binary trees, display mode indeed cost **nothing**: the
same level and the same reticulation number in all 600 cases.

That is a statement about the *minimum*, though, not about any particular
network — and the paper says so explicitly, giving a two-tree instance whose
minimum cluster network does not display both trees. Cass runs into this
routinely: in **268 of those 600 cases (45%) the cluster-mode network did not
display the two input trees**, even though a network with the same number of
reticulations that does display them exists.

```bash
phylocass examples/two_trees.newick                  # 2 reticulations, displays-trees=False
phylocass examples/two_trees.newick --display-trees  # 2 reticulations, displays-trees=True
```

So for two binary trees, switch the option on: it is free, and without it you
often do not get the trees. The result needs both trees binary and on the same
taxon set — it fails for three or more trees, and for two non-binary trees.

For three and four trees the models genuinely diverge, and then the option
costs something: one or two extra reticulations in 20 of 31 random cases.

Two caveats worth stating plainly:

- **It is a heuristic, with no optimality guarantee.** The paper's exactness
  results are about cluster mode. Worse, its decomposition theorem is about
  *level*, and the paper notes explicitly that the analogous statement for
  reticulation number does **not** hold — so splitting the problem along the
  incompatibility graph, which is what makes Cass fast, is itself an
  approximation when you are counting reticulations. Treat the number as an
  upper bound produced by a heuristic, not as the hybridization number.
- **It is slower, and climbs higher.** A stronger acceptance test rejects more
  candidates, so the search reaches further up the levels; expect to raise
  `max_level` and to want a `time_limit`.

If your input trees do not all have the same taxon set, the check still
compares cluster sets, which is the natural reading but is not the same as
displaying each tree on its own taxon set. Restrict the trees to their common
taxa first if that distinction matters to you.

## Trees with different taxon sets

Cass needs clusters of a single taxon set. A tree missing some taxa only gives
**partial** clusters: `{a,b}` seen on `{a,b,c}` says a and b group against c,
and says *nothing* about the taxa that tree lacks. Handing such a cluster to
Cass as if it were a full cluster asserts more than the tree does, so
`cass_from_trees` refuses input whose trees disagree about the taxa, and says
why.

Setting `z_closure` completes them first, using the Z-closure of Huson,
Dezulian, Klöpper and Steel (2004): partial splits are combined by a rule until
nothing new appears, and whichever have grown to cover the whole taxon set are
handed to the ordinary algorithm.

```bash
phylocass examples/partial_trees.newick              # refuses, and explains
phylocass examples/partial_trees.newick --z-closure  # completes, then builds
```

```python
from phylocass import CassOptions, cass_from_trees, read_trees

trees = read_trees("(((a,b),c),d); (((a,c),b),e);")   # neither has all 5 taxa
result = cass_from_trees(trees, CassOptions(z_closure=True))

result.level                      # 1
result.z_closure.dropped          # clusters that stayed partial and were lost
```

### The rule

For partial splits `S₁ = A₁|Ã₁` and `S₂ = A₂|Ã₂`, the Z-rule fires when exactly
one of the four intersections is empty — which is where the name comes from:

> if `A₁∩A₂ ≠ ∅`, `A₂∩Ã₁ ≠ ∅`, `Ã₁∩Ã₂ ≠ ∅` and **`A₁∩Ã₂ = ∅`**,
> produce `A₁|(Ã₁∪Ã₂)` and `(A₁∪A₂)|Ã₂`.

Rooted clusters need one adjustment. A split does not say which side is which,
but a cluster does — and every input tree shares the root, so conceptually each
tree gains a root taxon that always sits opposite the cluster. That pins the
orientation and makes `Ã₁∩Ã₂ ≠ ∅` automatic. Since `Aᵢ` may be *either* part of
split *i*, working through the combinations leaves two rules (one never fires,
because the root taxon is on both root sides; one is another with the pair
swapped). Writing `R = X \ C` for the root side of cluster `C` known on `X`:

| | condition | consequence |
| --- | --- | --- |
| **overlapping** | `C₁∩C₂ ≠ ∅`, `C₂∩R₁ ≠ ∅`, `C₁∩R₂ = ∅` | `C₁` is known on `X₁∪R₂`; `C₁∪C₂` is a cluster known on `C₁∪X₂` |
| **disjoint** | `C₁∩C₂ = ∅`, `C₁∩R₂ ≠ ∅`, `C₂∩R₁ ≠ ∅` | `C₁` is known on `X₁∪C₂` |

Widening the taxon set is what does the real work — it is how a partial cluster
becomes full. The disjoint rule matters more than it looks: most clusters of one
tree are disjoint from most clusters of another, and without it very little
completes.

PhyloCass keeps only the widest taxon set per cluster, so the state only grows
and the closure is a genuine fixed point — order-independent, and therefore
deterministic. That differs from the "replace the two inputs" reading of the
rule, which is order-dependent and is why implementations usually randomise the
order and repeat.

### What to expect

Z-closure is **incomplete**: it does not derive every cluster the input
implies, which is what motivated the later M- and Y-rules. Clusters that never
reach the full taxon set are dropped, so information is lost twice over.

Measured on random trees restricted to subsets of their taxa, then rebuilt —
checking each input tree against the *sub-network on that tree's taxa*, which
is the right question for partial input:

| input | input trees recovered | clusters dropped |
| --- | --- | --- |
| 6–8 taxa, 2–3 trees, 1 taxon missing | 78% | 39% |
| 8–10 taxa, 3–4 trees, 1 taxon missing | 93% | 40% |
| 6–8 taxa, 2–3 trees, 1–2 missing | 72% | 47% |
| 8–10 taxa, 4–5 trees, 1–3 missing | 89% | 59% |

More trees and more overlap give the rule more to work with. Sparse overlap
recovers less. Treat the result as a supernetwork-style summary, not as a
reconstruction.

`displays_partial_trees(network, trees)` is the check used above: it restricts
the network to each tree's taxa with PhyloZoo's `subnetwork` and then asks
whether that displays the tree. Use it rather than `displays_input_trees()`
when the taxon sets differ — comparing raw cluster sets asks something
stronger and slightly wrong.

> D. Huson, T. Dezulian, T. Klöpper, M. Steel. *Phylogenetic super-networks
> from partial trees.* IEEE/ACM TCBB **1**(4):151–158, 2004.

## What the guarantees are

| input | Cass gives you |
| --- | --- |
| clusters that fit a tree | that tree, level 0 |
| a level-1 network exists | a level-1 network (**proved**) |
| a level-2 network exists | a level-2 network (**proved**) |
| otherwise | a low-level network, in practice far below galled-network methods — but **not proved optimal** |
| `display_trees=True` | a network displaying the input trees — **heuristic**, no optimality claim |
| `z_closure=True` | trees on differing taxon sets accepted — **heuristic**, and lossy: what Z-closure cannot complete is dropped |

The paper conjectures Cass is optimal at every level and proves a decomposition
theorem supporting that, but only levels 1 and 2 are settled.
`CassOptions.max_level` (default 4) bounds how far the search climbs.

## Module map

| file | contents |
| --- | --- |
| [`clusters.py`](src/phylocass/clusters.py) | compatibility, incompatibility graph, ST-sets, `Collapse` |
| [`workgraph.py`](src/phylocass/workgraph.py) | the mutable search representation, and the hand-off to PhyloZoo |
| [`treebuild.py`](src/phylocass/treebuild.py) | the unique tree representing a compatible cluster set |
| [`io.py`](src/phylocass/io.py) | multi-tree reading, cluster extraction, verification |
| [`plot.py`](src/phylocass/plot.py) | drawing the input trees and the network together |
| [`zclosure.py`](src/phylocass/zclosure.py) | completing partial clusters from trees on differing taxon sets |
| [`cass.py`](src/phylocass/cass.py) | the algorithm |
| [`cli.py`](src/phylocass/cli.py) | command line |

## Performance

Cass runs in `O(|X|^(3k+2) · |C|)` time for fixed level *k* — polynomial, but
the exponent bites. `CassOptions` gives you three brakes:

```python
from phylocass import cass, CassOptions

result = cass(clusters, taxa, CassOptions(
    max_level=3,        # give up above this level
    time_limit=60.0,    # seconds per conflicting component, shared across levels
    max_networks=5000,  # intermediate networks kept per subproblem
))
```

`time_limit` is one budget for a whole component, not a fresh allowance at each
level the search climbs through. When a component exhausts its budget or
exceeds `max_level`, `cass` raises `RuntimeError` naming the component's size
and which limit it hit.

### Which datasets is this usable on?

The taxon count is almost irrelevant. What decides the running time is the
**largest conflicting component** — the incompatibility-graph component that is
left after unseparated sets have been collapsed — and its level.

Two regimes, measured by sampling networks, taking the clusters they display
and handing them back to Cass. First the pessimal one: a *simple* level-k
network, where every taxon hangs off a single blob so nothing collapses and the
component is the entire dataset.

Median seconds:

| taxa in one blob | level 1 | level 2 | level 3 | level 4 |
| ---: | ---: | ---: | ---: | ---: |
| 10 | 0.05 | 0.08 | 0.24 | 0.23 |
| 20 | 0.21 | 0.80 | 1.6 | 8.8 |
| 30 | 0.46 | 3.5 | 6.2 | |
| 40 | 1.0 | 7.3 | 50 | |
| 50 | 1.6 | 14 | 81 | |
| 60 | 2.6 | 30 | *gave up* | |
| 80 | 5.1 | 88 | | |
| 100 | 7.7 | 108 | | |
| 120 | 13.5 | 182 | | |

Every completed run returned a network of level ≤ *k* representing every input
cluster. A blank is untested, not a failure. Levels 1–2 and the 30-taxa-and-up
part of level 3 are medians of three samples at a 300 s budget; the 10- and
20-taxa entries for levels 3 and 4 are medians of five samples at 60 s, from a
second run that went after the small sizes specifically.

**Medians hide a long tail.** Occasional instances take a hundred times the
median: at level 3 on 10 taxa the median is 0.24 s and the slowest of five
samples hit the 60 s cap; at level 4 on 12 taxa, 0.40 s against 41 s. Set a
`time_limit` and treat the medians as typical, not as a bound.

Now the realistic one: a tree of blobs, i.e. many small conflicts scattered
across the taxa rather than one enormous tangle.

With blobs of up to level 3, the hardest of the three runs:

| taxa | reticulations | blobs | largest component | median |
| ---: | ---: | ---: | ---: | ---: |
| 100 | 61 | 32 | 4 blocks | 0.10 s |
| 200 | 132 | 69 | 4 blocks | 0.36 s |
| 400 | 271 | 132 | 4 blocks | 0.44 s |
| 800 | 543 | 275 | 4 blocks | 1.34 s |

**800 taxa with 543 reticulations in 1.3 s.** The component never exceeds four
blocks however large the dataset gets, because collapsing removes everything
that is not in genuine conflict.

Capping the blobs at level 1 or 2 instead gives fewer reticulations and slightly
lower times — 800 taxa comes to 267 reticulations in 0.88 s and 385 in 0.98 s.
The taxon count barely matters in any of the three; the component size is what
does, and it stays at four blocks throughout.

So, as a rule of thumb:

- **conflicts local** (each disagreement involves a handful of taxa) — hundreds
  of taxa are fine; measured with blobs up to level 3, and the level of an
  individual blob barely shows, because what the search sees is a component of
  a few blocks either way;
- **one big conflict, level 1** — 120 taxa in 13 s, and still climbing slowly;
- **one big conflict, level 2** — comfortable to around 100 taxa, minutes at 120;
- **one big conflict, level 3** — around 50 taxa;
- **one big conflict, level 4** — around 20 taxa.

Note how much better level 1 does than its own bound: `O(|X|^5)` would have
twelve times the taxa costing a quarter of a million times the time, and the
measured factor is 270 — roughly `n^2.3`. The bound counts an exhaustive
search, while Cass stops at the first network that works.

### What that means in practice

The `|X|` in the bound is misleading, because the decomposition means cost is
driven by the largest *conflicting component* after unseparated sets have been
collapsed — not by the taxon count. Sampling level-2 networks, taking the
clusters they display and handing them back to Cass:

| leaves | median time | hardest component | level ≤ 2 |
| ---: | ---: | ---: | :---: |
| 10 | 0.02 s | 7 blocks | 35/35 |
| 20 | 0.06 s | 11 blocks | 34/34 |
| 40 | 0.17 s | 15 blocks | 30/30 |
| 80 | 0.42 s | 19 blocks | 30/30 |
| 160 | 2.0 s | 26 blocks | 31/31 |

Sixteen times the taxa costs about thirty times the time, because the component
Cass actually searches grows far more slowly than the input. A dataset whose
disagreements are local stays fast however many taxa it has.

That effect is starker still on networks with many *small* conflicts. Sampling
multi-blob networks (a tree of blobs, each blob a level-≤2 generator) and
holding the reticulation count near ten while the taxa grow:

| taxa | 10 | 20 | 40 | 80 | 160 |
| --- | ---: | ---: | ---: | ---: | ---: |
| median time | 0.003 s | 0.006 s | 0.010 s | 0.011 s | 0.019 s |

Sixteen times the taxa for six times the time, because each blob collapses to a
component of at most four blocks no matter how large the dataset is.

### Verification scales too

Checking a *finished* network — `represents_input()`, `displays_input_trees()` —
means asking which clusters it displays, and the obvious way to do that
enumerates one tree per switching, so `2 ** (total reticulations)`. That is
fine for a handful of reticulations and hopeless past twenty.

PhyloCass does it one biconnected component at a time instead. Switchings in
different blobs are independent, and the taxa below a cut-edge are the same
under *every* switching — both parents of any reticulation down there lie in
that same blob, so nothing can cut it off. A node's descendant set therefore
varies only with its own blob's switching, and the cost becomes a sum of
`2 ** (reticulations in one blob)`:

| taxa | reticulations | blobs | blob-wise | trees the naive way needs |
| ---: | ---: | ---: | ---: | ---: |
| 20 | 14 | 8 | 0.002 s | 16,384 |
| 80 | 37 | 24 | 0.005 s | 1.4 × 10¹¹ |
| 160 | 78 | 51 | 0.013 s | 3.0 × 10²³ |
| 320 | 155 | 104 | 0.033 s | 4.6 × 10⁴⁶ |

`displays_input_trees()` decomposes the same way, which needs one more fact:
a cluster that is *not* switching-invariant belongs to exactly one blob, so the
blobs can be satisfied independently. `softwired_clusters_via_displayed_trees`
and `displays_trees_via_displayed_trees` are the naive versions, kept as the
references the fast ones are tested against.

The exponent does bite when the *conflict* is large rather than the dataset:
level 4 on 6–7 mutually conflicting taxa takes seconds to tens of seconds, and
that is inherent rather than an implementation problem.

## Tests

```bash
pip install -e ".[dev]"
pytest
```

Four kinds of test carry the weight:

- **Round-trips.** Generate a random network of level ≤ 2, take the clusters it
  displays, run Cass on them, and require a network of no higher level
  displaying every one of them. This tests the paper's exactness guarantee
  directly rather than comparing against hard-coded output. Run at 20 leaves as
  well as at small sizes; over 375 sampled 20-leaf level-2 networks, Cass
  returned a level-≤2 network representing every input cluster every time, and
  in three cases found a level-1 network — the sampled network's clusters did
  not need its second reticulation.

- **Generator-based sampling.** The round-trip above builds its networks by
  hanging a fresh leaf below each new reticulation, which always gives one blob
  and puts ~95% of reticulations directly above a leaf — the very shape Cass's
  own construction produces. [`tests/blobsampler.py`](tests/blobsampler.py)
  samples properly instead: it expands a random *tree of blobs*, each internal
  node becoming a simple level-k network built from one of PhyloZoo's level-k
  generators (`all_level_k_generators`, 1/1/4/65 for k = 0…3) with one
  attachment point per child.

  That covers what the simple sampler misses — networks come out with a median
  of 6 blobs and up to 12, 46% of reticulations sit above whole subtrees rather
  than leaves, and 97% of networks nest a blob directly below another blob's
  reticulation. Across 499 such networks on 20 taxa with blobs up to level 1, 2
  and 3, Cass represented every input cluster and never exceeded the source
  network's level.
- **Agreement with PhyloZoo.** The search's own softwired-cluster and level
  routines are checked against PhyloZoo's `displayed_trees` and `level` on
  random networks, and every finished network must pass `validate()`. The
  blob-wise shortcut described under *Performance* is pinned the same way:
  1481 random networks agreed on their cluster sets, and 9937 display queries
  over 1157 networks agreed exactly, with no mismatches either way.
- **The paper's example.** The nine-taxon cluster set from Figure 1 reproduces
  the level-2, two-reticulation network the paper reports, where the
  galled-network algorithm needs four reticulations.
- **Determinism.** Repeated runs, shuffled input, and subprocesses started
  under different `PYTHONHASHSEED` values must all agree.

The worked example above is tested too, down to the eNewick string quoted in
this README, so the documentation cannot drift from the code. Plotting tests
skip themselves if the `viz` extra is not installed.

## Licence

MIT. PhyloZoo is MIT-licensed as well; if you use this, cite both the Cass
paper above and the PhyloZoo preprint
([bioRxiv:10.64898/2026.06.09.731120](https://doi.org/10.64898/2026.06.09.731120)).
