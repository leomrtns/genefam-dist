"""Native and legacy interfaces for tree-signal calculations.

The native API has no DendroPy or NumPy dependency. Legacy high-level objects
and helpers are imported lazily so existing callers retain their public names.
"""

from ._treesignalc import (
    TreeSpace,
    from_tree_spaces,
    from_tree_spaces_pvalue,
    from_tree_spaces_rescale,
)


_LEGACY_UTILS = {
    "lowlevel_create_random_tree_string",
    "lowlevel_generate_spr_trees_string",
    "lowlevel_randomise_trees_with_spr_string",
    "lowlevel_calculate_spectrum_from_tree_strings",
    "lowlevel_calculate_spectrum_from_tree_strings_rescale",
    "lowlevel_calculate_spectrum_from_tree_strings_pvalue",
}

__all__ = [
    "TreeSpace",
    "from_tree_spaces",
    "from_tree_spaces_rescale",
    "from_tree_spaces_pvalue",
    "TreeSignal",
    *_LEGACY_UTILS,
]


def __getattr__(name):
    """Load the DendroPy/NumPy-based legacy API only when requested."""
    try:
        if name == "TreeSignal":
            from .treesignal import TreeSignal

            globals()[name] = TreeSignal
            return TreeSignal
        if name in _LEGACY_UTILS:
            from . import utils

            value = getattr(utils, name)
            globals()[name] = value
            return value
    except ModuleNotFoundError as error:
        if error.name in {"dendropy", "numpy"}:
            raise ModuleNotFoundError(
                "The legacy TreeSignal API requires optional dependencies; "
                "install them with `python -m pip install 'treesignal[legacy]'`."
            ) from error
        raise
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
