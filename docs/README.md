# Tree signatures with branch lengths

`genefam-dist` describes a gene family by its distances to a common panel of
reference species trees. The `treesignal` Python package is the recommended
entry point: it keeps parsed native C trees alive between calculations and
does not require a separate C library installation or subprocesses.

Use [the tutorial](tutorial.md) to compare trees with branch lengths and to run
a directory containing one inferred gene-family tree per file. The
[API reference](api.md) describes reusable native objects and the analysis
interface. The [technical report](technical-report.md) defines the distances,
their treatment of paralogs and missing species, and the limits of interpreting
branch lengths as reconciliation costs.

The original notebooks, manuscript, figures, and poster have been preserved in
[`notebooks/`](https://github.com/leomrtns/genefam-dist/tree/native-tree-api/notebooks).
They document the historical interface and may require its optional dependencies
and external programs. The new tutorial uses the current Python interface.

## Installation

From a checkout, with Python 3.11 or later, a C compiler, and Python headers:

```bash
python -m pip install .
```

For development and the test suite:

```bash
python -m pip install --editable '.[test]'
python -m pytest
```

If Python was built with a compiler absent from your machine, select the
available compiler explicitly, for example `CC=gcc python -m pip install .`.
The native API requires neither DendroPy nor NumPy. Install `.[legacy]` for the
historical `TreeSignal` class, and the notebook dependencies separately when
executing the archived notebooks.

## Reading these pages locally

Docsify renders these Markdown pages in a browser. From the repository root:

```bash
python -m http.server 8000 --directory docs
```

Open the server in your own browser. The Docsify theme and scripts load from
`cdn.jsdelivr.net`, so the rendered site needs access to that host. The Markdown
files can also be read directly without JavaScript. For GitHub Pages, select
the branch containing this work and `/docs` as the source directory.
