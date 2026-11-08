"""Check native length distances against an independent DendroPy oracle."""

import math

import dendropy
import pytest
from dendropy.calculate import treecompare

from treesignal import TreeSpace, _treesignalc


def oracle_distances(gene_newick, species_newick, normalize=False):
    namespace = dendropy.TaxonNamespace()
    trees = [
        dendropy.Tree.get(
            data=newick,
            schema="newick",
            rooting="force-unrooted",
            preserve_underscores=True,
            taxon_namespace=namespace,
        )
        for newick in (gene_newick, species_newick)
    ]
    common = sorted(
        set(node.taxon.label for node in trees[0].leaf_node_iter())
        & set(node.taxon.label for node in trees[1].leaf_node_iter())
    )
    for tree in trees:
        tree.retain_taxa_with_labels(common)
        tree.deroot()
        tree.suppress_unifurcations()
        tree.seed_node.edge_length = 0.0

    vectors = []
    for tree in trees:
        matrix = tree.phylogenetic_distance_matrix()
        vectors.append(
            [
                matrix.distance(namespace.get_taxon(a), namespace.get_taxon(b))
                for i, a in enumerate(common)
                for b in common[i + 1 :]
            ]
        )
    if normalize:
        vectors = [
            [value / (sum(vector) / len(vector)) for value in vector]
            for vector in vectors
        ]
        for tree in trees:
            total = sum(edge.length for edge in tree.postorder_edge_iter())
            for edge in tree.postorder_edge_iter():
                edge.length /= total
    rms = math.sqrt(
        sum((a - b) ** 2 for a, b in zip(*vectors)) / len(vectors[0])
    )
    if len(common) == 2:
        # There is a single unrooted edge after suppressing the root. DendroPy's
        # split dictionary keeps only one of the complementary terminal edges
        # in this degenerate case, so use the full path length as its weight.
        difference = abs(vectors[0][0] - vectors[1][0])
        return rms, difference, difference
    return (
        rms,
        treecompare.weighted_robinson_foulds_distance(*trees),
        treecompare.euclidean_distance(*trees),
    )


PAIRS = [
    (
        "((A:1,B:2):3,(C:4,D:5):6):99;",
        "(B:2,A:1,(D:5,C:4):9):1000;",
    ),
    (
        "((A:1,B:2):3,(C:4,D:5):6);",
        "((D:5,A:2):2,(C:4,B:2):4);",
    ),
    (
        "(((A:1,B:2):3,X:4):5,(C:6,D:7):8);",
        "((D:7,C:6):4,((B:2,Y:8):3,A:1):5);",
    ),
    (
        "(((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6):0.7,"
        "((E:0.8,F:0.9):1,G:1.1):1.2);",
        "(((G:0.5,C:0.4):0.3,(A:0.2,F:0.1):0.6):0.7,"
        "((B:1.1,E:1):0.9,D:0.8):1.2);",
    ),
    (
        "(((A:1,X:2):3,B:4):5,(Y:6,Z:7):8);",
        "(((A:2,U:3):4,B:5):6,(V:7,W:8):9);",
    ),
    (
        "((((A:1,X:2):3,B:4):5,C:6):7,(Y:8,Z:9):10);",
        "((A:2,B:3):4,(C:5,U:6):7);",
    ),
]


@pytest.mark.parametrize("gene,species", PAIRS)
@pytest.mark.parametrize("normalize", [False, True])
def test_native_metrics_match_dendropy(gene, species, normalize):
    expected = oracle_distances(gene, species, normalize)
    actual = _treesignalc.from_tree_spaces_lengths(
        TreeSpace(gene), TreeSpace(species), normalize=normalize
    )
    assert actual == pytest.approx(expected, abs=1e-12)


def test_terminal_edge_difference_has_known_distances():
    gene = TreeSpace("((A:1,B:1):1,(C:1,D:1):1);")
    species = TreeSpace("((A:2,B:1):1,(C:1,D:1):1);")

    assert _treesignalc.from_tree_spaces_lengths(gene, species) == pytest.approx(
        (math.sqrt(0.5), 1.0, 1.0)
    )


def test_zero_edges_are_valid_raw_distances():
    gene = TreeSpace("((A:0,B:0):0,(C:0,D:0):0);")
    species = TreeSpace("((A:0,C:0):0,(B:0,D:0):0);")

    assert _treesignalc.from_tree_spaces_lengths(gene, species) == (0.0, 0.0, 0.0)


@pytest.mark.parametrize("scale", [1e-150, 1e150])
def test_finite_extreme_scales_preserve_normalized_distances(scale):
    first, second = PAIRS[1]
    expected = oracle_distances(first, second, normalize=True)
    trees = []
    for newick in (first, second):
        tree = dendropy.Tree.get(data=newick, schema="newick")
        for edge in tree.postorder_edge_iter():
            if edge.length is not None:
                edge.length *= scale
        trees.append(TreeSpace(tree.as_string(schema="newick")))

    actual = _treesignalc.from_tree_spaces_lengths(*trees, normalize=True)
    assert all(math.isfinite(value) for value in actual)
    assert actual == pytest.approx(expected, abs=1e-12)


def test_zero_normalized_trees_remain_zero_vectors():
    zero = TreeSpace("((A:0,B:0):0,(C:0,D:0):0);")

    assert _treesignalc.from_tree_spaces_lengths(zero, zero, normalize=True) == (
        0.0, 0.0, 0.0
    )
