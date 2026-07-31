"""Sampling networks by expanding a tree of blobs, using PhyloZoo's generators.

``conftest.random_network`` builds a network by hanging a fresh leaf below each
new reticulation. That is easy but narrow: it always yields a single blob,
nearly every reticulation ends up directly above a leaf -- which is exactly the
shape Cass's own construction produces -- and the level-2 topologies are not
covered evenly.

This module samples in two stages instead:

1. sample a **tree of blobs**: a random rooted tree whose internal nodes are
   placeholders;
2. expand each internal node into a **blob** -- a simple level-k network built
   from a level-k *generator* -- with one attachment point per child, and
   splice the children in.

Each blob becomes a biconnected component of the result, so networks come out
with several blobs of mixed level. Because a child spliced onto a hybrid side
can itself be a blob, blobs nest below other blobs' reticulations, which the
tree-plus-reticulations sampler cannot produce at all.

PhyloZoo supplies the hard part. ``all_level_k_generators(k)`` enumerates the
level-k generators (1, 1, 4 and 65 for k = 0, 1, 2, 3) by the R1/R2 rules of
Gambette et al., and ``attach_leaves_to_generator`` turns a generator plus an
assignment of taxa to sides into a network. What is added here is the
randomisation and the blob-tree scaffolding; PhyloZoo has no random network
sampler of its own.
"""

from __future__ import annotations

import random
from functools import lru_cache

from phylozoo import DirectedPhyNetwork
from phylozoo.core.network.dnetwork.classifications import level as pz_level
from phylozoo.core.network.dnetwork.generator.attachment import (
    attach_leaves_to_generator,
)
from phylozoo.core.network.dnetwork.generator.construction import (
    all_level_k_generators,
)
from phylozoo.core.network.dnetwork.generator.side import HybridSide

from phylocass.workgraph import WorkGraph

__all__ = [
    "generator_catalogue",
    "min_leaves",
    "sample_blob",
    "random_blob_network",
]


@lru_cache(maxsize=None)
def generator_catalogue(level: int) -> tuple:
    """The level-``level`` generators, ordered so sampling is reproducible."""
    if level < 1:
        return ()
    return tuple(
        sorted(all_level_k_generators(level), key=lambda g: str(sorted(map(str, g.sides))))
    )


def min_leaves(generator) -> int:
    """Fewest attachment points this generator can carry.

    Every hybrid side needs exactly one taxon, and PhyloZoo refuses a network
    with fewer than two taxa in total.
    """
    return max(2, sum(1 for s in generator.sides if isinstance(s, HybridSide)))


def sample_blob(
    generator, n_leaves: int, rng: random.Random, attempts: int = 40
) -> DirectedPhyNetwork | None:
    """A blob from ``generator`` with exactly ``n_leaves`` attachment points.

    Leaves are labelled ``slot0``..``slotN``; the caller replaces them with
    whole subtrees. Returns ``None`` if this generator cannot carry that many.
    """
    hybrid_sides = [s for s in generator.sides if isinstance(s, HybridSide)]
    edge_sides = [s for s in generator.sides if not isinstance(s, HybridSide)]
    if n_leaves < min_leaves(generator) or not edge_sides and n_leaves != len(hybrid_sides):
        return None

    spare = n_leaves - len(hybrid_sides)
    if spare < 0:
        return None

    for _ in range(attempts):
        names = [f"slot{i}" for i in range(n_leaves)]
        rng.shuffle(names)
        pool = iter(names)

        side_taxa = {side: [next(pool)] for side in hybrid_sides}
        buckets: dict = {side: [] for side in edge_sides}
        for _ in range(spare):
            buckets[edge_sides[rng.randrange(len(edge_sides))]].append(next(pool))
        side_taxa.update(buckets)

        try:
            blob = attach_leaves_to_generator(generator, side_taxa)
        except Exception:
            continue
        # a bundle of parallel edges survives unless all but one is subdivided
        if len(list(blob.edges)) != len(set(blob.edges)):
            continue
        if pz_level(blob) != generator.level:
            continue
        return blob
    return None


def _random_blob_tree(rng: random.Random, n_taxa: int):
    """A random rooted tree on ``n_taxa`` leaves; internal nodes get expanded.

    Arity varies between two and four so blobs of different sizes fit.
    """
    children: dict[int, list[int]] = {leaf: [] for leaf in range(n_taxa)}
    nodes = list(range(n_taxa))
    next_id = n_taxa
    while len(nodes) > 1:
        take = min(len(nodes), rng.choice([2, 2, 2, 3, 3, 4]))
        picked = rng.sample(nodes, take)
        children[next_id] = picked
        nodes = [v for v in nodes if v not in picked] + [next_id]
        next_id += 1
    return nodes[0], children


def random_blob_network(
    rng: random.Random,
    n_taxa: int,
    max_level: int = 2,
    blob_probability: float = 0.55,
    names=None,
) -> WorkGraph:
    """Sample a network by expanding a tree of blobs.

    Every blob has level between 1 and ``max_level``, so the network's level is
    at most ``max_level`` by construction -- no filtering needed. Internal nodes
    that are not expanded stay plain split nodes.
    """
    from conftest import taxon_names

    names = names or taxon_names(n_taxa)
    root, children = _random_blob_tree(rng, n_taxa)
    catalogue = {k: generator_catalogue(k) for k in range(1, max_level + 1)}

    net = WorkGraph()

    def build(v: int) -> int:
        kids = children[v]
        if not kids:
            return net.new_node(block=frozenset({names[v]}))

        built = [build(k) for k in kids]

        if rng.random() < blob_probability:
            levels = [k for k in range(1, max_level + 1) if catalogue[k]]
            rng.shuffle(levels)
            for k in levels:
                options = [g for g in catalogue[k] if min_leaves(g) <= len(built)]
                rng.shuffle(options)
                for generator in options:
                    blob = sample_blob(generator, len(built), rng)
                    if blob is not None:
                        order = list(built)
                        rng.shuffle(order)
                        return _graft_blob(net, blob, order)

        parent = net.new_node()
        for child in built:
            net.add_edge(parent, child)
        return parent

    build(root)
    return net


def _graft_blob(net: WorkGraph, blob: DirectedPhyNetwork, subroots) -> int:
    """Copy ``blob`` into ``net``, its ``slotN`` leaves replaced by subtrees."""
    mapping: dict = {}
    for v in blob.nodes:
        label = blob.get_label(v)
        if label is not None and label.startswith("slot"):
            mapping[v] = subroots[int(label[4:])]
        else:
            mapping[v] = net.new_node()
    for u, v in blob.edges:
        net.add_edge(mapping[u], mapping[v])
    return mapping[blob.root_node]
