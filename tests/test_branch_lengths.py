import math

import pytest

from treesignal import TreeSpace, _treesignalc


TREE = "((A:1,B:1):1,(C:1,D:1):1);"
DOUBLED = "((A:2,B:2):2,(C:2,D:2):2);"


def test_length_variants_are_retained_in_order():
    references = TreeSpace(TREE + DOUBLED + TREE, rooted=True)
    assert references.ntrees == references.ndistinct == 3
    values = _treesignalc.from_tree_spaces_lengths(TreeSpace(TREE), references)
    assert values[:3] == (0., 0., 0.)
    assert values[3:6] == pytest.approx((math.sqrt(12), 6, math.sqrt(8)))
    assert values[6:] == (0., 0., 0.)


def test_serialization_preserves_matrix_and_native_label_order():
    original = TreeSpace(TREE)
    copied = TreeSpace(original.newick)
    matrix = copied.patristic_matrix()
    expected = original.patristic_matrix()
    for i, label in enumerate(original.taxa):
        for j, other in enumerate(original.taxa):
            assert matrix[copied.taxa.index(label)][copied.taxa.index(other)] == expected[i][j]


def test_scale_normalization_and_zero_lengths():
    assert _treesignalc.from_tree_spaces_lengths(TreeSpace(TREE), TreeSpace(DOUBLED), normalize=True) == (0., 0., 0.)
    zero = TreeSpace("((A:0,B:0):0,(C:0,D:0):0);")
    assert _treesignalc.from_tree_spaces_lengths(zero, zero, normalize=True) == (0., 0., 0.)


def test_missing_lengths_and_no_shared_taxa_raise():
    with pytest.raises(ValueError):
        _treesignalc.from_tree_spaces_lengths(TreeSpace(TREE), TreeSpace("((A,B),(C,D));"))
    with pytest.raises(ValueError):
        _treesignalc.from_tree_spaces_lengths(TreeSpace(TREE), TreeSpace("((E:1,F:1):1,(G:1,H:1):1);"))
    with pytest.raises(ValueError):
        TreeSpace("((A,B),(C,D));").patristic_matrix()


def test_mixed_length_collections_do_not_replace_missing_lengths():
    references = TreeSpace(TREE + "(A,(B,(C,D)));")
    assert references.ndistinct == 2
    assert not references.has_branch_lengths
    with pytest.raises(ValueError):
        references.patristic_matrix(index=1)


def test_generation_is_native_seeded_and_exact_count():
    a = _treesignalc.generate_tree_space(8, 7, seed=42)
    b = _treesignalc.generate_tree_space(n_species=8, chain_size=7, n_spr=1, seed=42)
    assert isinstance(a, TreeSpace)
    assert a.ntrees == a.ndistinct == 7
    assert a.has_branch_lengths
    assert a.newick == b.newick
    assert set(a.taxa) == {f"T{i}" for i in range(1, 9)}
    assert TreeSpace(a.newick).ntrees == 7
    with pytest.raises(OverflowError):
        _treesignalc.generate_tree_space(8, 2, seed=-1)
    with pytest.raises(ValueError):
        _treesignalc.generate_tree_space(8, 2, n_spr=0)


def test_overflow_cannot_silently_normalize_to_zero():
    enormous = TreeSpace("((A:1e308,B:1e308):1e308,(C:1e308,D:1e308):1e308);")
    with pytest.raises(ValueError):
        _treesignalc.from_tree_spaces_lengths(enormous, TreeSpace(TREE), normalize=True)
    with pytest.raises(ValueError):
        enormous.patristic_matrix()
