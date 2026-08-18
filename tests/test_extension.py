import dendropy
import treesignal
from treesignal import _treesignalc


def test_native_extension_is_package_local():
    assert _treesignalc.__name__ == "treesignal._treesignalc"


def test_legacy_high_level_api_is_still_available():
    assert treesignal.TreeSignal
    assert treesignal.lowlevel_calculate_spectrum_from_tree_strings


def test_generate_spr_trees_round_trip():
    tree_string = _treesignalc.generate_spr_trees(4, 3, 1)

    assert isinstance(tree_string, str)
    parsed = dendropy.TreeList.get(
        data=tree_string,
        schema="newick",
        preserve_underscores=True,
    )
    # The initial random tree plus three SPR-derived trees.
    assert len(parsed) == 4


def test_distances_from_tree_strings():
    gene_tree = "(((A,B),(C,D)),((E,F),(G,H)));"
    species_trees = (
        "(((A,B),(C,D)),((E,F),(G,H)));"
        "(((A,C),(B,D)),((E,G),(F,H)));"
    )

    distances = _treesignalc.fromtrees(gene_tree, species_trees)

    # Six distance measures for each reference species tree.
    assert len(distances) == 12
    assert all(isinstance(value, float) for value in distances)
    assert distances[:6] == (0.0,) * 6
