# PhyloCass

A Python implementation of **Cass**, the algorithm from

> L. van Iersel, S. Kelk, R. Rupp, D. Huson.
> *Phylogenetic Networks Do not Need to Be Complex: Using Fewer Reticulations to Represent Conflicting Clusters.*
> Bioinformatics **26**(12):i124–i131, 2010 ([arXiv:0910.3082](https://arxiv.org/abs/0910.3082))

Cass takes a set of rooted phylogenetic trees, collects every cluster those
trees display, and builds a single rooted phylogenetic network that represents
**all** of those clusters in the softwired sense — while keeping the *level* of
the network as low as it can.

The level of a network is the largest number of reticulations inside any single
biconnected component. Minimising it, rather than the total reticulation
number, is what lets Cass produce networks that stay readable even when the
input trees disagree a lot.

## Install

```bash
pip install -e .
```

No dependencies beyond the standard library. Python 3.10+.

## Command line

```bash
# from a file of Newick trees, one tree per ';'
phylocass examples/conflicting_trees.newick

# from an explicit cluster list, one cluster per line
phylocass --format clusters examples/figure1_clusters.txt

# or from stdin
echo '((a,b),(c,d)); ((a,c),(b,d));' | phylocass
```

Output is the network in extended Newick, with reticulations written as
`#H1`, `#H2`, …:

```
# taxa=4 clusters=4 level=2 reticulations=2 verified=True
(((a,(b)#H1),(c)#H2),((d,#H1),#H2));
```

The summary line goes to stderr and the network to stdout, so
`phylocass in.newick > out.enewick` gives you just the network.

## Library

```python
from phylocass import read_trees, cass_from_trees

trees = read_trees("""
    (((a,b),c),d);
    (((a,c),b),d);
""")

result = cass_from_trees(trees)
print(result.level)                  # 1
print(result.reticulation_number)    # 1
print(result.to_enewick())           # (d,((c,(a,(b)#H1)),#H1));
print(result.represents_input())     # True
```

Or start from clusters directly:

```python
from phylocass import cass

clusters = {frozenset("ab"), frozenset("bc"), frozenset("abc")}
result = cass(clusters, taxa=frozenset("abcd"))
```

`result.represents_input()` re-checks the finished network against the input
clusters from scratch, independently of the search. It is worth calling.

## What the guarantees are

| input | Cass gives you |
| --- | --- |
| clusters that fit a tree | that tree, level 0 |
| a level-1 network exists | a level-1 network (**proved**) |
| a level-2 network exists | a level-2 network (**proved**) |
| otherwise | a low-level network, in practice much lower than galled-network methods — but **not proved optimal** |

The paper conjectures Cass is optimal for every level and proves a
decomposition theorem supporting that, but only levels 1 and 2 are settled.
`CassOptions.max_level` (default 4) bounds how far the search will climb.

## How it works

Two halves, matching the paper.

**Decomposition (Section 3, Steps 1–4)** — in [`cass.py`](src/phylocass/cass.py).
The incompatibility graph `IG(C)` has the clusters as nodes and joins two
clusters when they overlap without nesting. Each non-trivial connected
component becomes an independent subproblem: collapse the maximal unseparated
subsets of its support, solve it, build the tree on everything left over, then
splice each component's network into that tree at the LCA of its taxa. Theorem
1 of the paper is what makes this safe — if any level-*k* network exists, one
respecting this decomposition exists.

**Simple level-*k* networks (Section 4, Algorithm 1)** — also in `cass.py`, as
`_SimpleSearch`. Loop over the taxa; for each one, delete it from every cluster
and collapse the maximal ST-sets of what remains. Do that *k* times and only
two taxa are left. Then walk back out: decollapse each ST-set leaf into its
strict subtree, and hang the removed taxon below a new reticulation fed by
every possible pair of edges, keeping any network that represents the clusters
at that level.

Two details from the paper that are easy to miss and that PhyloCass
implements:

- a **dummy taxon** `d` is offered alongside the real taxa at every round, and
  when `d` is the one removed the collapse step is deliberately skipped. This
  is what lets Cass build reticulations of indegree 3: hang `d`, hang the real
  taxon, then delete `d` and contract the edge between the two reticulations.
- every tree built at the base of the recursion gets a **dummy root** above its
  real root, so that a reticulation edge can also be attached above the root.
  Both dummies are stripped on output.

### Module map

| file | contents |
| --- | --- |
| [`clusters.py`](src/phylocass/clusters.py) | compatibility, incompatibility graph, ST-sets, `Collapse` |
| [`network.py`](src/phylocass/network.py) | the DAG, softwired clusters, biconnected components, level, eNewick |
| [`treebuild.py`](src/phylocass/treebuild.py) | the unique tree representing a compatible cluster set |
| [`newick.py`](src/phylocass/newick.py) | Newick reader |
| [`cass.py`](src/phylocass/cass.py) | the algorithm |
| [`cli.py`](src/phylocass/cli.py) | command line |

### A note on the representation

Collapsing taxa into groups happens over and over in this algorithm. Rather
than minting composite taxon names, PhyloCass keeps every cluster written in
terms of the *original* taxa and tracks the current level of collapsing as a
separate partition into **blocks** — a block being a `frozenset` of original
taxa that currently acts as one taxon. Restriction, deletion and the
"does this network represent C?" test then need no translation layer at all.
Leaf labels are blocks; the empty block marks a dummy leaf, which therefore
contributes nothing to any descendant set for free.

## Performance

Cass runs in `O(|X|^(3k+2) · |C|)` time for fixed level *k* — polynomial, but
the exponent bites. Level 2 is comfortable; level 4 and up on more than a
handful of conflicting taxa can run for a long time. `CassOptions` gives you
three brakes:

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

Note that the decomposition means the cost is driven by the largest *conflicting
component*, not by the total number of taxa — a 200-taxon dataset whose
disagreements are local stays fast. On this machine, levels 0–2 come back in
well under a second for the sizes above; level 4 on 6–7 mutually conflicting
taxa takes seconds to tens of seconds, which is the `O(|X|^(3k+2))` exponent
showing up rather than anything avoidable.

## Tests

```bash
pip install -e ".[dev]"
pytest
```

Besides unit tests for each module, the suite includes **round-trip** tests:
generate a random network of level ≤ 2, take the clusters it displays, run Cass
on them, and require that Cass comes back with a network of no higher level
that displays every one of those clusters. That is a direct test of the paper's
exactness guarantee rather than a test against hard-coded output. The
nine-taxon example from Figure 1 of the paper is checked too — Cass finds the
level-2, two-reticulation network the paper reports, where the galled-network
algorithm needs four reticulations.

## Licence

MIT.
