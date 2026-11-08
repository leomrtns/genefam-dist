# Python API and command line

## Recommended interface

```python
from treesignal import TreeSpace, TreeAnalysis, analyze_files, generate_references
```

The native extension is compiled into the Python package. No external executable,
DendroPy object, or NumPy array is required by this interface.

### `TreeSpace(newick, rooted=False)`

Parse a string containing one or more Newick trees and own their native C
structures until the Python object is destroyed. All trees in one space must
use the same unique leaf labels. Gene copies should use distinct labels and an
explicit gene-to-species mapping for the analysis interface.

`ntrees` counts input trees; `ndistinct` describes stored native alternatives.
Branch-length inputs retain their alternatives in input order, including trees
with identical topology and different lengths. Historical topology-only spaces
may compact repeated topologies. `rooted` controls the topology-only comparison
of alternatives; the branch-length split metrics are unrooted.

`taxa` exposes the native leaf-label order. `newick` serializes the stored panel,
and `has_branch_lengths` indicates whether lengths were supplied. Use
`patristic_matrix(index=0)` for a matrix of raw path sums for one tree in the
native taxon order. The diagonal is zero. A stem above the root is excluded
from these sums.

Lengths must be finite and nonnegative. If any edge is annotated with a length,
supply lengths for all non-stem edges; incomplete lengths raise `ValueError`.
The supported Newick subset has at least three leaves, binary internal branching,
and an optional trifurcating root for an unrooted tree. Taxon labels must be
unique, simple, unquoted labels without whitespace. Internal support labels are
accepted and ignored. Quotes, square-bracket comments, and Nexus input are not
supported; remove such annotations before passing output from tree-inference
programs to `TreeSpace`. Excessive nesting (over 512 levels) is rejected.

### `TreeAnalysis(references)`

Accept a `TreeSpace` or Newick string and retain the reference panel for reuse.

```python
analysis.compare(
    gene,
    metrics=("patristic", "wrf", "kf"),
    normalize=False,
    gene_to_species=None,
    aggregation="min",
)
```

`gene` is a native space or Newick string containing one gene-family tree.
Returns a dictionary keyed by metric name; each value is an ordered tuple of
reference scores. `gene_to_species` maps unique gene-copy labels to reference
species labels. With no mapping the label itself identifies the species.

| Metric | Definition | Applicable data |
| --- | --- | --- |
| `patristic` | RMS discrepancy between species-pair path matrices | Single-copy or mapped paralogous families |
| `wrf` | L1 difference between unrooted split lengths | Single-copy species labels |
| `kf` | L2 (Kuhner–Felsenstein) split-length difference | Single-copy species labels |
| `topology` | Original six topological scores | Existing legacy gene-label mapping semantics |

The `topology` result is flattened in reference-major order, with six entries
per reference: duplications, losses, deep coalescences, RF, Hdist, approximate
SPR. `normalize=False` returns raw counts; `normalize=True` uses the existing
theoretical upper bounds. These scores intentionally ignore branch lengths. The explicit mapping
argument is for the new patristic calculation; it does not change the legacy
reconciliation algorithm's species matching.

Comparisons use species present in both trees; at least two are required for
the length-aware scores. Ordinary split metrics are not defined by this API for
multiple gene copies per species. Specify only `patristic` for that case.
`aggregation="min"` uses the minimum cross-copy path distance;
`aggregation="mean"` averages all cross-copy distances for the species pair.

### Normalization

`normalize=False` retains input branch-length units. With `normalize=True`:

- `patristic`: divide each paired species matrix by its own mean upper-triangle
  entry after projecting to the shared species. Compute the RMS of these scaled
  entries. An all-zero matrix remains zero.
- `wrf` and `kf`: divide each induced tree's edge lengths by its own total
  length before comparing its split vector. An all-zero tree remains zero.

These operations remove overall scale independently in the two inputs. They
are **not** `phangorn`'s sum-of-lengths normalization of the raw wRF score.
Both raw and normalized scores can be useful; record which was used. Neither
normalization corrects for inference bias, taxon sampling, or nonuniform rates.

### `analyze_files(paths, references, output=None, ...)`

Read exactly one Newick tree from each path, reuse one reference object, and
return rows with fields `family`, `reference`, `metric`, and `value`. Supply
`output` to write a TSV. Reference indices are 1-based. The `topology` metric is
expanded into named rows `duplication`, `loss`, `deep_coalescence`, `rf`, `hdist`,
and `spr`. The metric, normalization, mapping, and aggregation
arguments follow `TreeAnalysis.compare`. Errors identify the relevant input file.
Use deterministic path order and retain reference/input files alongside output.

### `generate_references(species, count=20, n_spr=1, seed=0, branch_length=None)`

Create a native reference panel with the supplied species labels and SPR
generation parameters, without spawning a standalone executable. A positive
seed is reproducible; zero chooses native random initialization. Generated trees
have unit unrooted edge lengths by default; the two edges incident to the
artificial degree-two root are each half a unit. Supplying `branch_length`
scales all these synthetic edge lengths. Lengths are not estimated from data
by reference generation.

## Native numerical interface

For callers working directly with native objects:

```python
from treesignal import from_tree_spaces_lengths

values = from_tree_spaces_lengths(gene, references, normalize=False)
```

Returns reference-major triples `(patristic_rms, weighted_rf, kuhner_felsenstein)`.
The native routine uses unique exact taxon labels and induced common-species
trees; it does not aggregate paralogs. Its `normalize=True` follows the same
conventions as `TreeAnalysis`: paired-matrix mean scaling for patristic RMS,
and induced total-tree-length scaling for wRF/KF.

`from_tree_spaces`, `from_tree_spaces_rescale`, and
`from_tree_spaces_pvalue` preserve the established topology-only interfaces.
Their input is native `TreeSpace` objects. `fromtrees*` in
`treesignal._treesignalc` accepts strings for older callers. `TreeSignal` and
the low-level legacy Python helpers remain available after installing
`treesignal[legacy]`.

## Command line

The two operations in `python -m treesignal` cover reference generation and
analysis. Run `python -m treesignal --help`, `analyze --help`, or
`references --help` for the accepted options. The [tutorial](tutorial.md) gives
complete file-based examples. The older C programs and OMA preparation scripts
are retained for existing users; they are independent of the recommended
analysis interface.
