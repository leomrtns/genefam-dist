import gc
import subprocess
import sys

import pytest

from treesignal import (
    TreeSpace,
    from_tree_spaces,
    from_tree_spaces_pvalue,
    from_tree_spaces_rescale,
)
from treesignal import _treesignalc


GENE_NEWICK = "(((A,B),(C,D)),((E,F),(G,H)));"
SPECIES_NEWICK = (
    "(((A,B),(C,D)),((E,F),(G,H)));"
    "(((A,C),(B,D)),((E,G),(F,H)));"
)


def test_tree_space_parses_once_and_exposes_metadata():
    gene = TreeSpace(GENE_NEWICK)
    species = TreeSpace(SPECIES_NEWICK, rooted=True)

    assert gene.ntrees == gene.ndistinct == 1
    assert gene.rooted is False
    assert species.ntrees == species.ndistinct == 2
    assert species.rooted is True
    assert "trees=2" in repr(species)


def test_native_and_legacy_string_apis_match():
    gene = TreeSpace(GENE_NEWICK)
    species = TreeSpace(SPECIES_NEWICK, rooted=True)

    assert from_tree_spaces(gene, species) == _treesignalc.fromtrees(
        GENE_NEWICK, SPECIES_NEWICK
    )
    assert from_tree_spaces_rescale(
        gene, species
    ) == _treesignalc.fromtrees_rescale(GENE_NEWICK, SPECIES_NEWICK)


def test_native_objects_remain_unchanged_after_mutating_algorithms():
    gene = TreeSpace(GENE_NEWICK)
    species = TreeSpace(SPECIES_NEWICK, rooted=True)
    expected = from_tree_spaces(gene, species)

    pvalues = from_tree_spaces_pvalue(gene, species, 10)

    assert len(pvalues) == 24
    assert from_tree_spaces(gene, species) == expected
    assert gene.ntrees == 1
    assert species.ntrees == 2


def test_native_functions_reject_non_tree_space_arguments():
    species = TreeSpace(SPECIES_NEWICK, rooted=True)

    with pytest.raises(TypeError):
        from_tree_spaces(GENE_NEWICK, species)


def test_tree_space_owns_and_releases_native_memory():
    for _ in range(100):
        tree_space = TreeSpace(GENE_NEWICK)
        assert tree_space.ntrees == 1
        del tree_space
    gc.collect()


def test_importing_native_api_does_not_import_optional_dependencies():
    result = subprocess.run(
        [
            sys.executable,
            "-c",
            "import sys, treesignal; "
            "assert 'dendropy' not in sys.modules; "
            "assert 'numpy' not in sys.modules; "
            "assert treesignal.TreeSpace",
        ],
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr
