"""PhyloCass -- the Cass algorithm, built on PhyloZoo.

Cass takes a set of rooted phylogenetic trees, collects every cluster they
display, and constructs a single rooted phylogenetic network representing all
of those clusters in the softwired sense, while keeping the *level* of the
network -- the largest number of reticulations inside any biconnected
component -- as small as it can.

Networks are PhyloZoo ``DirectedPhyNetwork`` objects, so the result drops
straight into the rest of the PhyloZoo ecosystem:

    >>> from phylocass import read_trees, cass_from_trees
    >>> trees = read_trees("(((a,b),c),d);  (((a,c),b),d);")
    >>> result = cass_from_trees(trees)
    >>> result.level
    1
    >>> result.network                      # doctest: +ELLIPSIS
    <phylozoo...DirectedPhyNetwork object at ...>

Reference: L. van Iersel, S. Kelk, R. Rupp, D. Huson, "Phylogenetic Networks
Do not Need to Be Complex: Using Fewer Reticulations to Represent Conflicting
Clusters", Bioinformatics 26(12):i124-i131, 2010.
"""

from .cass import CassOptions, CassResult, cass, cass_from_trees, cass_simple
from .clusters import (
    block_partition,
    collapse,
    incompatibility_graph,
    is_compatible,
    is_separated,
    is_st_set,
    maximal_st_sets,
    maximal_unseparated_sets,
    nontrivial_components,
    unseparated_closure,
)
from .io import (
    clusters_of_trees,
    hardwired_clusters,
    read_cluster_file,
    read_tree_file,
    read_trees,
    softwired_clusters,
)
from .treebuild import build_tree
from .workgraph import WorkGraph

__version__ = "0.2.0"

__all__ = [
    "CassOptions",
    "CassResult",
    "WorkGraph",
    "block_partition",
    "build_tree",
    "cass",
    "cass_from_trees",
    "cass_simple",
    "clusters_of_trees",
    "collapse",
    "hardwired_clusters",
    "incompatibility_graph",
    "is_compatible",
    "is_separated",
    "is_st_set",
    "maximal_st_sets",
    "maximal_unseparated_sets",
    "nontrivial_components",
    "read_cluster_file",
    "read_tree_file",
    "read_trees",
    "softwired_clusters",
    "unseparated_closure",
]
