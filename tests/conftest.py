"""Shared helpers: random trees and networks built as working graphs."""

import random

from phylocass.workgraph import WorkGraph


def fs(s):
    return frozenset(s)


def singletons(s):
    return [frozenset({t}) for t in s]


def shown(sets):
    """Sets of taxa as sorted tuples, for readable assertions."""
    return sorted(tuple(sorted(s)) for s in sets)


def taxon_names(n: int) -> list[str]:
    """``a``..``z``, then ``aa``, ``ab``, ... -- spreadsheet-column style.

    ``chr(ord("a") + i)`` only works up to 26 taxa; past that it walks into
    punctuation and control characters, which are still distinct strings but
    are not labels anyone would want to see in Newick.
    """
    names = []
    for i in range(n):
        name, i = "", i
        while True:
            name = chr(ord("a") + i % 26) + name
            i = i // 26 - 1
            if i < 0:
                break
        names.append(name)
    return names


def random_tree(rng: random.Random, names) -> WorkGraph:
    """A random binary rooted tree on ``names``."""
    tree = WorkGraph()
    nodes = [tree.new_node(block=frozenset({n})) for n in names]
    while len(nodes) > 1:
        i, j = rng.sample(range(len(nodes)), 2)
        parent = tree.new_node()
        tree.add_edge(parent, nodes[i])
        tree.add_edge(parent, nodes[j])
        nodes = [n for k, n in enumerate(nodes) if k not in (i, j)] + [parent]
    return tree


def random_network(rng: random.Random, n_tree_leaves: int, n_ret: int) -> WorkGraph:
    """A random tree with ``n_ret`` extra leaves hung below new reticulations.

    Note what this does and does not sample. Each reticulation is fed by two
    uniformly chosen edges and carries a *fresh leaf*, so nearly every
    reticulation ends up directly above a leaf -- which is the shape Cass's own
    construction produces. Level-2 networks whose reticulations sit above whole
    subtrees are under-represented, and the four level-2 generators are not
    sampled uniformly. Round-trip tests built on this are still sound, because
    the sampled network itself witnesses that a level-2 solution exists; it is
    the coverage that is narrow, not the check.
    """
    all_names = taxon_names(n_tree_leaves + n_ret)
    names = all_names[:n_tree_leaves]
    net = random_tree(rng, names)
    for i in range(n_ret):
        edges = net.edges()
        # allow attaching above the root as well
        old_root = net.root
        new_root = net.new_node()
        net.add_edge(new_root, old_root)
        edges.append((new_root, old_root))

        e1, e2 = rng.sample(edges, 2)
        w1 = net.subdivide(*e1)
        w2 = net.subdivide(*e2)
        reticulation = net.new_node()
        leaf = net.new_node(block=frozenset({all_names[n_tree_leaves + i]}))
        net.add_edge(reticulation, leaf)
        net.add_edge(w1, reticulation)
        net.add_edge(w2, reticulation)
        net.tidy()
    return net


def nontrivial(clusters, taxa):
    return {c for c in clusters if 1 < len(c) < len(taxa)}
