# Branch-length-aware gene-family signatures

## Scope

A tree signature compares each gene-family tree with the same ordered panel of
candidate species trees. It is a feature representation of agreement with the
panel, not a posterior probability that a candidate is the species tree. Trees
with different taxon coverage or evolutionary rates need explicit handling
before their signatures can be compared across families.

The extension adds branch-length information through two complementary
representations: patristic matrices and vectors of weighted splits. The original
duplication, loss, deep-coalescence, RF, Hdist, and approximate SPR calculations
retain their topological meanings. Changing the interpretation of those
existing columns silently would invalidate historical analyses.

## Patristic matrix discrepancy

For leaves `i` and `j` of a gene tree `G`, let `p_G(i,j)` be the sum of edge
lengths on the unique path joining them. If `m(i)` maps gene copies to species,
define the species-level gene distance by

```text
D_G(a,b) = min { p_G(i,j) : m(i)=a, m(j)=b },  a != b.
```

This is the minimum cross-copy aggregation used by one family of methods in
[super_distance](https://github.com/leomrtns/super_distance). It selects the
closest cross-species copies, which can reduce the contribution of ancient
paralog divergences. It is a descriptive reduction: hidden paralogy, horizontal
transfer, erroneous short edges, or loss of the closest copies can still bias
it. Pairwise minima need not correspond to a single ortholog selection or to
an additive tree metric. In particular, these reductions can violate the
triangle inequality. This implementation does not claim statistical consistency
under a duplication/loss model. It is a discrepancy between reduced
representations, rather than a metric on the space of all paralogous trees:
different gene trees can reduce to the same matrix.

For the shared species set `X` of size `n`, compare this matrix with the reference
species-tree patristic matrix `D_S` using the upper-triangle root mean square:

```text
patristic_rms(G,S) = sqrt( sum_{a<b in X} (D_G(a,b)-D_S(a,b))^2 / choose(n,2) ).
```

Counting each unordered pair once avoids duplicating the off-diagonal entries;
the diagonal is identically zero and is excluded. Dividing by the pair count
makes the scale interpretable as a typical pairwise discrepancy. It does not
remove sampling bias from incomplete families. The unnormalised value has the
same units as the input branch lengths, usually substitutions per site.

The API also offers the mean over all cross-copy path distances as an explicit
alternative. The minimum and mean describe different aspects of a family and
should not be mixed between runs without recording the choice.

Missing species contribute no pair. A gene label without a species assignment
must be corrected rather than assigned a zero distance. At least two represented
species are necessary for a patristic comparison. References should use the
same species labels throughout a panel.

## Weighted split distances

For a single-copy gene tree and a reference on the common species set, let
`B` be the union of their unrooted bipartitions. Each split `b` has an edge length
`l_T(b)` in tree `T`, or zero if absent. Pendant edges are included. After pruning
unused species, degree-two vertices are suppressed by summing their incident
path lengths. The two child edges of an artificial degree-two root therefore
describe the same unrooted split and their lengths are added. A stem above the
root is not part of any path between leaves and is excluded.

The two complementary distances are

```text
weighted_rf(G,S) = sum_{b in B} abs(l_G(b) - l_S(b))
kuhner_felsenstein(G,S) = sqrt(sum_{b in B} (l_G(b) - l_S(b))^2).
```

These are respectively the L1 and L2 distances between split-length vectors.
The latter is the Kuhner–Felsenstein branch-score distance. Their definitions
follow `wRF.dist` and `KF.dist` in
[phangorn](https://github.com/KlausVigo/phangorn/blob/master/R/treedist.R).
The L2 score gives proportionally more weight to large individual edge
discrepancies. Even identical topologies can have nonzero scores; incompatible
splits contribute their full lengths against zero. A zero-length resolution is
equivalent to a contracted edge for these scores.

For paralogous families, the species-level matrix remains well defined, but
there is no unique ordinary split-length vector on species. Mapping two copies
to the same species can put that species on both sides of a gene split. The
ordinary wRF/KF formulas cannot resolve this. Report these metrics as unavailable
for such families, and use the patristic score and the existing mul-tree
topological distances. Do not silently select an arbitrary copy or average
incompatible splits.

## Rate and scale

Raw values are meaningful only when gene and reference branches have comparable
units and calibrations. A gene evolving twice as fast as another can have a
different patristic signature even if both have the same topology. A scale-free
comparison deliberately removes overall length information and asks about
relative distances or relative edge allocation instead.

The precise normalization used by the API is described in [the API
reference](api.md). It must not be confused with `phangorn::wRF.dist(normalize=TRUE)`,
which divides the *raw L1 score* by the sum of the two trees' total lengths.
These are different operations. A tree or matrix with zero total length also
requires an explicit convention; dividing by zero is not acceptable.

## Why branch lengths are not generic reconciliation weights

The existing reconciliation scores count inferred duplications, losses, or extra
lineages over possible rootings of the unrooted gene tree. These are discrete
events or lineage counts. A weighted objective would need to be optimized again
over candidate rootings; multiplying an already minimized count need not choose
the rooting minimizing the weighted score.
Multiplying a duplication count by the adjacent gene branch length has no
general likelihood justification. A branch normally measures substitutions per
site, while duplication and loss intensities are rates per unit time. A long
gene edge can indicate time, a high substitution rate, or estimation error.
It is not automatically a costly duplication or evidence of a loss.

A defensible extension is possible when the trees carry the extra biological
assumptions needed for it. For a dated, rooted species tree and a temporally
valid reconciliation `R`, let `k_R(e,t)` be the number of gene lineages occupying
species edge `e` at time `t`. The *extra-lineage exposure* is

```text
C_exposure(R,S) = sum_e integral_{t in e} max(k_R(e,t)-1, 0) dt.
```

When the lineage number is constant along each edge, this becomes the weighted
extra-lineage count `sum_e length(e) * max(k_e-1,0)`. It measures excess lineage
time, with units of time, and can distinguish a short from a long period of
deep coalescence. It is an explicit descriptive objective, not the complete
multispecies-coalescent likelihood. Under Kingman's coalescent model (reference
7), the waiting-time contribution is instead proportional to

```text
sum_e integral choose(k_R(e,t),2) / (2 N_e(t)) dt
```

for a diploid population-size convention, together with event-density terms.
Lineage counts change at embedded coalescences, so edge duration alone is not
enough to calculate the integral in general.

Likewise, a duplication–loss birth–death model combines event log-rates and
integrated lineage hazards, schematically
`-N_dup log(lambda) - N_loss log(mu) + integral (lambda+mu) k(t) dt`, with further
conditioning and observation terms depending on the model. Arbitrarily weighting
event counts by substitution lengths is not this likelihood. Dated gene/species
trees, event placements, sampling assumptions, and rate parameters would be
required for a valid implementation.

Consequently, this release implements branch-length scores where they have
standard definitions and retains reconciliation as counts. It does not label
an ad hoc product a weighted reconciliation distance. The exposure objective
above provides a concrete target for a later dated-reconciliation extension.

## Hdist and SPR

Hdist compares splits through an assignment cost, and the current SPR score is
an approximation to the minimum number of prune-and-regraft moves. Neither
becomes a standard branch-length metric by multiplying its final integer by
tree length. A weighted SPR objective needs an explicit cost for each move and
a rule for redistributing lengths after pruning and regrafting. A length-aware
split assignment likewise needs a specified matching cost. The weighted split
and patristic scores avoid inventing these semantics while preserving the
established topology-only outputs.

## Representation and reproducibility

`TreeSpace` is a Python object owning a native C `topology_space` pointer. The
pointer is opaque to callers and released with the object's lifetime. Reusing
the object avoids repeated Newick parsing. Native calculations copy mutable
reconciliation state where required by the existing algorithms.

Branch-length reference alternatives must not be collapsed merely because their
topologies are identical. Input order and the association between a reference
and its output columns are part of the analysis. Generated references are a
panel of candidate topologies, not inferred species trees; their branch-length
assignment is a modeling choice. Reference generation and all analysis parameters
should be recorded with the results.

The executable notebook and tests cover small examples with known distances.
They demonstrate behavior, not validation on a large OMA data set. OMA download
and alignment scripts remain upstream data-preparation tools; branch-length tree
inference and biological quality control remain the user's responsibility.

## References and implementation sources

1. Robinson, D. F. and Foulds, L. R. (1981). Comparison of phylogenetic trees.
   *Mathematical Biosciences* 53:131–147.
   [DOI](https://doi.org/10.1016/0025-5564(81)90043-2).
2. Kuhner, M. K. and Felsenstein, J. (1994). A simulation comparison of phylogeny
   algorithms under equal and unequal evolutionary rates. *Molecular Biology and
   Evolution* 11:459–468.
   [DOI](https://doi.org/10.1093/oxfordjournals.molbev.a040126).
3. Steel, M. A. and Penny, D. (1993). Distributions of tree comparison metrics—some
   new results. *Systematic Biology* 42:126–141.
   [DOI](https://doi.org/10.1093/sysbio/42.2.126).
4. Schliep, K. P. (2011). phangorn: phylogenetic analysis in R.
   *Bioinformatics* 27:592–593.
   [DOI](https://doi.org/10.1093/bioinformatics/btq706).
   [Distance documentation](https://klausvigo.github.io/phangorn/reference/treedist.html)
   and [source](https://github.com/KlausVigo/phangorn/blob/master/R/treedist.R).
5. de Oliveira Martins, L. et al. (2016). A Bayesian Supertree Model for Genome-Wide
   Species Tree Reconstruction. *Systematic Biology* 65:397–416.
   [DOI](https://doi.org/10.1093/sysbio/syu082).
6. [super_distance source and method descriptions](https://github.com/leomrtns/super_distance),
   especially its within-family minimum/mean patristic aggregation and scaling
   conventions. This release uses the matrix representation as inspiration; it
   does not embed or execute that program.
7. Kingman, J. F. C. (1982). The coalescent. *Stochastic Processes and their
   Applications* 13:235–248.
   [DOI](https://doi.org/10.1016/0304-4149(82)90011-4).
8. Szöllősi, G. J. et al. (2015). The inference of gene trees with species trees.
   *Systematic Biology* 64:e42–e62.
   [DOI](https://doi.org/10.1093/sysbio/syu048). This reviews model-based
   reconciliation beyond unweighted event-count objectives.
