"""PhyloCass -- the Cass algorithm for building low-level phylogenetic networks.

Cass takes a set of rooted phylogenetic trees, collects every cluster they
display, and constructs a single rooted phylogenetic network that represents
all of those clusters in the softwired sense, while keeping the *level* of the
network -- the largest number of reticulations inside any biconnected
component -- as small as it can.

    >>> from phylocass import read_trees, cass_from_trees
    >>> trees = read_trees("((a,b),(c,d));  ((a,c),(b,d));")
    >>> result = cass_from_trees(trees)
    >>> result.level
    1

Reference: L. van Iersel, S. Kelk, R. Rupp, D. Huson, "Phylogenetic Networks
Do not Need to Be Complex: Using Fewer Reticulations to Represent Conflicting
Clusters", Bioinformatics 26(12):i124-i131, 2010.
"""

from .cass import CassOptions, CassResult, cass, cass_from_trees, cass_simple
from .clusters import (
    clusters_of_trees,
    incompatibility_graph,
    is_compatible,
    is_separated,
    is_st_set,
    maximal_st_sets,
    maximal_unseparated_sets,
    nontrivial_components,
)
from .network import Network
from .newick import parse_newick, parse_newick_file, read_trees
from .treebuild import build_tree

__version__ = "0.1.0"

__all__ = [
    "CassOptions",
    "CassResult",
    "Network",
    "build_tree",
    "cass",
    "cass_from_trees",
    "cass_simple",
    "clusters_of_trees",
    "incompatibility_graph",
    "is_compatible",
    "is_separated",
    "is_st_set",
    "maximal_st_sets",
    "maximal_unseparated_sets",
    "nontrivial_components",
    "parse_newick",
    "parse_newick_file",
    "read_trees",
]
