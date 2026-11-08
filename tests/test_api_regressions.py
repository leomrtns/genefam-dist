"""Input failures must raise exceptions, never terminate the Python process."""

import subprocess
import sys

import pytest

from treesignal import TreeSpace, from_tree_spaces


def rejected_in_subprocess(expression):
    # Keep fatal-error and hang regressions isolated from the rest of pytest.
    program = (
        "from treesignal import TreeSpace, from_tree_spaces\n"
        "from treesignal import _treesignalc\n"
        "try:\n"
        f"    {expression}\n"
        "except ValueError:\n"
        "    pass\n"
        "else:\n"
        "    raise AssertionError('invalid input was accepted')\n"
    )
    result = subprocess.run(
        [sys.executable, "-c", program],
        capture_output=True,
        text=True,
        timeout=10,
        check=False,
    )
    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize(
    "newick",
    [
        "",
        "not a tree",
        "((A,B),(C,D);",
        "(A,B,C,D);",
        "((A,A),(B,C));",
        "((A,B),(C,D));((A,B),(C,X));",
    ],
)
def test_invalid_newick_raises_without_exiting(newick):
    rejected_in_subprocess(f"TreeSpace({newick!r})")


@pytest.mark.parametrize(
    "function",
    ["from_tree_spaces", "from_tree_spaces_rescale", "from_tree_spaces_pvalue"],
)
def test_unknown_gene_species_raises_without_exiting(function):
    rejected_in_subprocess(
        f"_treesignalc.{function}("
        "TreeSpace('((A,B),(C,X));'), TreeSpace('((A,B),(C,D));', rooted=True))"
    )


def test_negative_spr_count_raises_instead_of_hanging():
    rejected_in_subprocess("_treesignalc.generate_spr_trees(4, 3, -1)")


def test_reinitialization_keeps_native_objects_usable():
    gene = TreeSpace("((A,B),(C,D));")
    species = TreeSpace("((A,B),(C,D));", rooted=True)
    expected = from_tree_spaces(gene, species)
    gene.__init__("((D,C),(B,A));")
    assert from_tree_spaces(gene, species) == expected


def test_failed_reinitialization_preserves_previous_tree():
    gene = TreeSpace("((A,B),(C,D));")
    species = TreeSpace("((A,B),(C,D));", rooted=True)
    expected = from_tree_spaces(gene, species)
    with pytest.raises(ValueError):
        gene.__init__("not a tree")
    assert from_tree_spaces(gene, species) == expected


def test_a_tree_space_can_be_its_own_reference():
    tree = TreeSpace("(((A,B),(C,D)),((E,F),(G,H)));", rooted=True)
    assert from_tree_spaces(tree, tree) == (0.0,) * 6


def test_reference_reinitialization_does_not_reuse_stale_patristic_matrix():
    from treesignal import TreeAnalysis

    newick = "((A:1,B:1):1,(C:1,D:1):1);"
    references = TreeSpace(newick, rooted=True)
    analysis = TreeAnalysis(references)
    assert analysis.compare(newick, metrics="patristic")["patristic"] == (0.0,)
    references.__init__("((A:2,B:2):2,(C:2,D:2):2);", rooted=True)
    assert analysis.compare(newick, metrics="patristic")["patristic"][0] > 0.0
