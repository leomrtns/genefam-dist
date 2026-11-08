# A single Python workflow

Start with inferred Newick trees containing finite, nonnegative branch lengths.
The reference species trees must also contain branch lengths for the new
length-aware metrics. Gene and species branches should represent comparable
units, or you should deliberately choose normalized comparisons.

Use binary Newick trees (an unrooted root may have three children), with unique
unquoted leaf labels. Remove square-bracket comments and quoted labels from
tree-inference output before loading it. Internal support labels are accepted
but ignored. If lengths occur, every non-stem edge must have one. See the
[supported input format](api.md) for details.

## Compare small trees

```python
from treesignal import TreeAnalysis, TreeSpace

gene = TreeSpace("((A:1,B:1):1,(C:1,D:1):1);")
references = TreeSpace(
    "((A:1,B:1):1,(C:1,D:1):1);"
    "((A:2,B:1):1,(C:1,D:1):1);",
    rooted=True,
)
analysis = TreeAnalysis(references)
scores = analysis.compare(gene)
print(scores)

assert scores["patristic"][0] == 0.0
assert scores["wrf"] == (0.0, 1.0)
assert scores["kf"] == (0.0, 1.0)
```

The two references share a topology but differ in the pendant length for `A`.
They remain separate reference alternatives. The second patristic score is
`sqrt(3/6)`: three of the six species pairs differ by one unit. The wRF and KF
scores are both one because exactly one split length differs by one unit.

Reuse `analysis` for each new gene family. `TreeSpace` owns native structures;
the object calls pass opaque handles and do not repeatedly parse the reference
Newick string.

## Paralogs and absent species

Use unique gene-copy labels and an explicit dictionary identifying their species:

```python
gene = TreeSpace("(((A_1:0.5,A_2:0.5):0.5,B_1:1):1,(C_1:1,D_1:1):1);")
mapping = {"A_1": "A", "A_2": "A", "B_1": "B", "C_1": "C", "D_1": "D"}
scores = analysis.compare(
    gene,
    metrics=("patristic",),
    gene_to_species=mapping,
    aggregation="min",
    normalize=True,
)
```

`min` uses the closest cross-species copy pair for each species pair; `mean`
uses their average. Neither is an automatic orthology assignment. Missing
reference species are excluded from the comparison, with at least two species
needed. The ordinary wRF/KF scores are restricted to single-copy labels, so
select `patristic` for paralogous data. See the [technical report](technical-report.md)
for the assumptions behind each reduction.

## One gene-family tree per file

Assume you have already inferred branch-length trees from your homologous gene
sets. Put exactly one tree in each family file and your ordered reference panel
in `references.nwk`. In Python:

```python
from pathlib import Path
from treesignal import analyze_files

rows = analyze_files(
    sorted(Path("families").glob("*.nwk")),
    references=Path("references.nwk").read_text(),
    output="signatures.tsv",
    metrics=("patristic", "wrf", "kf"),
    normalize=True,
)
```

The equivalent command line uses the same analysis code:

```bash
python -m treesignal analyze \
  --references references.nwk \
  --genes 'families/*.nwk' \
  --output signatures.tsv \
  --metrics patristic wrf kf \
  --normalize
```

Output is a tidy table with columns `family`, `reference`, `metric`, and `value`.
Retain the reference file with the results so reference indices can be traced to
the original input order. Invalid input is reported with the family filename.

For gene-copy labels, create a two-column TSV mapping without a header (gene
label, then species label), and use only the patristic score:

```text
A_1	A
A_2	A
B_1	B
C_1	C
D_1	D
```

```bash
python -m treesignal analyze \
  --references references.nwk --genes 'families/*.nwk' \
  --mapping mapping.tsv --metrics patristic --aggregation min \
  --normalize --output signatures.tsv
```

The existing OMA download and alignment scripts are data preparation steps.
Their FoldMason guide trees are not automatically validated evolutionary trees:
infer and assess gene-family trees with a suitable model before this analysis.
Large OMA data acquisition is not required to execute the example.

## Generate a reference panel from Python

For topology exploration, reference generation is part of the same package:

```python
from treesignal import generate_references, TreeAnalysis

references = generate_references(
    ["A", "B", "C", "D"], count=20, n_spr=1, seed=42, branch_length=1.0
)
analysis = TreeAnalysis(references)
```

Or save the generated panel for reuse:

```bash
python -m treesignal references --species species.txt \
  --count 20 --spr 1 --seed 42 --branch-length 1 \
  --output references.nwk
```

`species.txt` contains one species label per line. SPR neighbours can revisit
the same topology; `count` is a requested panel size, not a promise of unique
trees. A positive seed makes generation reproducible; seed zero requests the
native random initialization. Generated trees have unit unrooted edge lengths
by default. `branch_length=1.0` explicitly requests that default; another value
scales those synthetic edges. Generation does not estimate rates, times, or a biological
species tree. For biological length comparisons, prefer inferred reference
lengths from an appropriate data set and model.

The optional archived C executables remain available, but they are no longer
needed to generate references or compare a collection of gene trees through
this workflow.

## Executable example

[`notebooks/010.branch_lengths.ipynb`](https://github.com/leomrtns/genefam-dist/blob/native-tree-api/notebooks/010.branch_lengths.ipynb)
contains a small self-contained version of the examples, including temporary
family files and a result table. It uses only the installed `treesignal` package
and the Python standard library; the notebook frontend itself requires Jupyter.
