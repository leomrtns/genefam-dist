import subprocess
import sys

import pytest

from treesignal import TreeSpace, _treesignalc


@pytest.mark.parametrize(
    "newick",
    [
        "",
        "not a tree",
        "(A,B);",
        "((A,B),C",
        "((A,B),);",
        "((A,A),C);",
        "((A,B),C);((A,B),D);",
        "((A,B),C);((A,B),C,D);",
        "((A,B,C),D);",
        "((A),B,C);",
        "((A:1,B:1):1,C);",
        "((A:1,B:1):1,C:nan);",
        "((A:1,B:1):1,C:inf);",
        "((A:1,B:1):1,C:-1);",
        "((A:1,B:1):1,C:1e10000);",
        "((A:1x,B:1):1,C:1);",
        "(('A name',B),C);",
        "((A name,B),C);",
        "((A[comment],B),C);",
        "((A,B),C); garbage",
    ],
)
def test_bad_newick_raises_value_error(newick):
    with pytest.raises(ValueError):
        TreeSpace(newick)


@pytest.mark.parametrize(
    "newick",
    [
        "((A,B),C);",
        "((A,B),C)",
        " ( ( A , B ) , C ) ; ",
        "((A:0,B:2e-3)95:1,C:2);",
        "(A:1,B:2,C:3);",
        "((A:1,B:2):3,C:4):5;",
        "((A,B),C);(C,(B,A));",
    ],
)
def test_supported_newick_is_accepted(newick):
    tree = TreeSpace(newick)
    assert tree.ntrees >= 1


def test_invalid_legacy_input_raises_without_exiting_python():
    for function in (
        _treesignalc.fromtrees,
        _treesignalc.fromtrees_rescale,
        _treesignalc.fromtrees_pvalue,
    ):
        with pytest.raises(ValueError):
            function("invalid", "((A,B),C);")


def test_invalid_input_preserves_existing_tree_space():
    tree = TreeSpace("((A,B),C);")
    with pytest.raises(ValueError):
        tree.__init__("invalid")
    assert tree.ntrees == 1


def test_bad_input_does_not_terminate_interpreter():
    result = subprocess.run(
        [sys.executable, "-c", "from treesignal import TreeSpace; "
         "\ntry: TreeSpace('((A,B),C);((A,B),D);')"
         "\nexcept ValueError: print('caught')"],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "caught"
