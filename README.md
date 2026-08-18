# genefam_dist

The main goal of *genefam_dist* is to describe the "tree signal" of a gene family, which is a vector of 'scores' describing how well
supported this gene family is by a set of reference species trees. You can download [here (PDF)](docs/Evolution2018_LeoMartinsPoster.pdf) 
a low-resolution version of the poster presented at [Evolution 2018](https://www.evolutionmontpellier2018.org/).

The software is composed of a low-level C library to handle mul-trees, a Python module, and standalone programs using the library.
The Python module is called `treesignal` and is the standard way of creating tree signatures. It uses
[DendroPy](https://www.dendropy.org/) for tree manipulation and NumPy for its result arrays. The standalone programs are described below.

The library is written in C, with functions to calculate distances between mul-trees (gene families) and species trees. 
It uses a modified internal biomcmc library, part of the [guenomu software](https://bitbucket.org/leomrtns/guenomu/). 

## Installation

### Python package

The Python package uses the native [CPython C API](https://docs.python.org/3/c-api/) directly. Its build compiles the
`genefam_dist` C sources into `treesignal._treesignalc`, so no separate library installation, hard-coded paths, or 
runtime linker configuration is needed. 
A C compiler and Python development headers are required. 

Install the dependency-free native API from the repository root:

```console
python -m pip install .
```

The reusable native API does not require DendroPy or NumPy:

```python
from treesignal import TreeSpace, from_tree_spaces

gene = TreeSpace("(((A,B),(C,D)),((E,F),(G,H)));")
species = TreeSpace(
    "(((A,B),(C,D)),((E,F),(G,H)));"
    "(((A,C),(B,D)),((E,G),(F,H)));",
    rooted=True,
)
distances = from_tree_spaces(gene, species)
```

`TreeSpace` owns the parsed C `topology_space` and can be reused across calls.
The original string-based functions remain available in `treesignal._treesignalc`
for compatibility, but parse their Newick arguments on every call. Install the
DendroPy/NumPy-based `TreeSignal` API with:

```console
python -m pip install '.[legacy]'
```

For development, use an editable installation and run the tests:

```console
python -m pip install --editable '.[test]'
python -m pytest
```

To use conda, the included environment file creates an isolated environment and installs the extension in editable mode:

```console
conda env create --file environment.yml
conda activate treesignal
python -m pytest
```

You can then use the package with `import treesignal`. The extension source is in `python/treesignal/_treesignalc.c`; new low-level wrappers can be added
there and the corresponding C implementation can be added under `lib/`.

### C library and standalone programs

The standalone programs use the Autotools build system. A clean out-of-tree build can be installed with:

```console
mkdir build
cd build
../configure --prefix="$HOME/.local"
make
make install
```

This standalone installation is independent of the Python package. See the notebooks in the [docs directory](docs) for examples.


## Programs
Here is a list of standalone programs that will be compiled and installed with autotools. Run them without arguments for a brief explanation
(they are simple but require a rigid structure for the arguments).

- ***gf_distmatrix_genetree_sptree*** given two files, one with a list of gene trees and one with a list of species trees,
  calculates a set of distances between all tree pairs

- ***gf_distsignal_genetree_sptree*** given two files, like the distmatrix program above, but all distances for one gene
  tree per line, as in the tree signature vector

- ***gf_generate_spr_trees*** generates a series of random trees with a given number of SPR operations apart from the next
  one.

- ***gf_spr_distance*** calculates the approximate unrooted SPR distance between consecutive pairs of trees in a nexus file
  (this program assumes single-labelled trees, that is, no paralogs etc. are allowed)

- ***gf_concatenate_trees*** takes two nexus trees (from e.g. [MrBayes](http://nbisweden.github.io/MrBayes) analyses) and merge them by
  excluding first trees and thinning, if needed. It also excludes less frequent trees and assumes unrooted topologies (i.e. it neglects root
  information), important since returns a *trprobs* nexus file. 

- ***gf_find_best_trees*** takes a (large) number of gene families and estimates patristic distance-based species trees, as well as tries to
  find species trees minimising all/some distances. This is an experimental program, intended as a replacement for the `bmc2_maxtree` program from 
  [guenomu](https://bitbucket.org/leomrtns/guenomu/).


## More info
More details to follow, 
- For running the programs, you can have a brief description by running them with no parameters.
- For more details about the calculations and theory, please check the [guenomu software](https://bitbucket.org/leomrtns/guenomu/) and the accompanying publication
  - A Bayesian Supertree Model for Genome-Wide Species Tree Reconstruction, Systematic Biology, 2016 [DOI: 10.1093/sysbio/syu082](http://dx.doi.org/10.1093/sysbio/syu082)
- If you have any question, comment or request please do not hesitate to contact me at leomrtns@gmail.com
- For an (not extensive) introduction on tree distances that I will be focusing on, please have a look at [this talk at UNIL](http://www.slideshare.net/leomrtns/comparing-phylogenetic-trees-20160616).

## Performance of the SPR approximation
The `dSPR` is an approximation to the unrooted SPR which was used on the [recombination detection model](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0002651)
(but faster than the original version and capable of working with arbitrarily large trees). 
The `Hdist` is related to the [clade alignment score](http://bioinformatics.oxfordjournals.org/content/22/1/117), that I rediscovered (with modifications) when speeding up the 
dSPR algortihm.

Notice that our current implementations of these distances can work on mul-trees on a very experimental setting (using the *extended species tree* concept from the 
[mulRF algorithm](https://www.ncbi.nlm.nih.gov/pubmed/25273112)), but the results below assume the same leaf set between trees. 

![](docs/performance_hdist.png)
![](docs/performance_spr_rf.png)

## License 
Copyright (C) 2016-today  [Leonardo de Oliveira Martins](https://github.com/leomrtns)

genefam-dist is free software; you can redistribute it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later
version (http://www.gnu.org/copyleft/gpl.html).
