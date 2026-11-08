/* 
 * This file is part of genefam-dist, a library for calculating distances between gene families and species trees. 
 * Copyright (C) 2016  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 *
 * genefam-dist is free software; you can redistribute it and/or modify it under the terms of the GNU General Public 
 * License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more 
 * details (file "COPYING" or http://www.gnu.org/copyleft/gpl.html).
 */

#include "genefam_dist.h"
#define NDISTS 6

/*! \brief calculate distances between gene tree and collection of sptrees, and storing into newly allocated distances[] */ 
bool generate_output_distances (topology_space gtree, topology_space stree, double **distances, int *n, bool rescale);
/*! \brief calculate distances between gene tree and random sptrees, accumulating into previously allocated pvalues[] */ 
int update_randomized_distances (topology gtree, topology stree, int n_obs, double *obs_dist, int *pvalues, int n_reps, double *maxdist);
/*! \brief defines upper bound for all distances based on theory, not empirical as p-value above */
bool calculate_max_distances (topology gt, double *maxd, splitset split);
/*! \brief copy a topology for calculations which mutate reconciliation or tree state */
static topology copy_topology_for_calculation (topology source);
/*! \brief wrap a copied topology in a minimal owning topology space */
static topology_space new_calculation_tree_space (topology source);

/* A split-length coordinate after pruning to a shared leaf set. Complementary
 * root edges and chains left by pruning are summed into the same coordinate. */
typedef struct {
  bipartition split;
  double length;
} length_split;

static int
induced_length_splits (topology tree, const int *leaf_ids, int nleaves, length_split **result, double *total)
{
  int i, j, k, count = 0;
  bipartition split = new_bipartition (nleaves);
  *result = (length_split*) biomcmc_malloc (tree->nnodes * sizeof (length_split));
  *total = 0.;
  if (!tree->traversal_updated) update_topology_traversal (tree);
  for (i = 0; i < tree->nnodes; i++) {
    if (tree->nodelist[i] == tree->root) continue; /* no root stem in an unrooted metric */
    bipartition_zero (split);
    for (j = 0; j < nleaves; j++)
      if (bipartition_is_bit_set (tree->nodelist[i]->split, leaf_ids[j])) bipartition_set (split, j);
    bipartition_count_n_ones (split);
    if (!split->n_ones || split->n_ones == nleaves) continue;
    bipartition_flip_to_smaller_set (split);
    for (k = 0; k < count; k++) if (bipartition_is_equal ((*result)[k].split, split)) break;
    if (k == count) {
      (*result)[count].split = new_bipartition_copy_from (split);
      (*result)[count++].length = 0.;
    }
    (*result)[k].length += tree->blength[i];
    *total += tree->blength[i];
  }
  del_bipartition (split);
  return count;
}

static void
free_length_splits (length_split *splits, int count)
{
  int i;
  for (i = 0; i < count; i++) del_bipartition (splits[i].split);
  free (splits);
}

int
genefam_module_treesignal_lengths (topology_space gene, topology_space species, bool normalize, double **output)
{
  int i, j, k, s, ncommon = 0, ng, ns, *gid, *sid;
  double gt, st, delta, rms, wrf, kf, gscale, sscale, gd, sd, gmean, smean, pairs;
  length_split *gs, *ss;
  distance_matrix gm, sm;
  topology g;

  if (!gene || !species || !gene->ndistinct || !species->ndistinct) return 0;
  g = gene->distinct[0];
  if (!g->blength) return 0;
  for (s = 0; s < species->ndistinct; s++) if (!species->distinct[s]->blength) return 0;
  gid = (int*) biomcmc_malloc (g->nleaves * sizeof (int));
  sid = (int*) biomcmc_malloc (g->nleaves * sizeof (int));
  for (i = 0; i < g->nleaves; i++) {
    j = lookup_hashtable (species->taxlabel_hash, g->taxlabel->string[i]);
    if (j >= 0) { gid[ncommon] = i; sid[ncommon++] = j; }
  }
  if (ncommon < 2) { free (gid); free (sid); return 0; }
  ng = induced_length_splits (g, gid, ncommon, &gs, &gt);
  if (!isfinite (gt)) { free_length_splits (gs, ng); free (gid); free (sid); return -1; }
  gm = new_distance_matrix_for_topology (g->nleaves);
  fill_distance_matrix_from_topology (gm, g, g->blength, true);
  *output = (double*) biomcmc_malloc (3 * species->ndistinct * sizeof (double));
  gscale = normalize && gt > 0. ? gt : 1.;
  for (s = 0; s < species->ndistinct; s++) {
    topology sp = species->distinct[s];
    ns = induced_length_splits (sp, sid, ncommon, &ss, &st);
    sm = new_distance_matrix_for_topology (sp->nleaves);
    fill_distance_matrix_from_topology (sm, sp, sp->blength, true);
    if (!isfinite (st)) {
      free_length_splits (ss, ns); del_distance_matrix (sm);
      free_length_splits (gs, ng); del_distance_matrix (gm); free (gid); free (sid);
      free (*output); *output = NULL; return -1;
    }
    sscale = normalize && st > 0. ? st : 1.;
    pairs = (double)ncommon * (ncommon - 1) / 2.;
    gmean = smean = 0.;
    if (normalize) for (i = 0; i < ncommon; i++) for (j = i + 1; j < ncommon; j++) {
      gd = gid[i] < gid[j] ? gm->d[gid[i]][gid[j]] : gm->d[gid[j]][gid[i]];
      sd = sid[i] < sid[j] ? sm->d[sid[i]][sid[j]] : sm->d[sid[j]][sid[i]];
      if (!isfinite (gd) || !isfinite (sd)) {
        free_length_splits (ss, ns); del_distance_matrix (sm);
        free_length_splits (gs, ng); del_distance_matrix (gm); free (gid); free (sid);
        free (*output); *output = NULL; return -1;
      }
      gmean += gd / pairs; smean += sd / pairs;
    }
    if (gmean <= 0.) gmean = 1.;
    if (smean <= 0.) smean = 1.;
    rms = wrf = kf = 0.;
    for (i = 0; i < ncommon; i++) for (j = i + 1; j < ncommon; j++) {
      gd = gid[i] < gid[j] ? gm->d[gid[i]][gid[j]] : gm->d[gid[j]][gid[i]];
      sd = sid[i] < sid[j] ? sm->d[sid[i]][sid[j]] : sm->d[sid[j]][sid[i]];
      delta = gd / gmean - sd / smean;
      rms = hypot (rms, delta);
    }
    for (i = 0; i < ng; i++) {
      for (j = 0; j < ns; j++) if (bipartition_is_equal (gs[i].split, ss[j].split)) break;
      delta = gs[i].length / gscale - (j < ns ? ss[j].length / sscale : 0.);
      wrf += fabs (delta); kf = hypot (kf, delta);
    }
    for (j = 0; j < ns; j++) {
      for (k = 0; k < ng; k++) if (bipartition_is_equal (gs[k].split, ss[j].split)) break;
      if (k == ng) { delta = ss[j].length / sscale; wrf += delta; kf = hypot (kf, delta); }
    }
    (*output)[3*s] = rms / sqrt (pairs);
    (*output)[3*s + 1] = wrf;
    (*output)[3*s + 2] = kf;
    free_length_splits (ss, ns);
    del_distance_matrix (sm);
  }
  free_length_splits (gs, ng);
  del_distance_matrix (gm);
  free (gid); free (sid);
  return 3 * species->ndistinct;
}

topology_space
genefam_module_generate_tree_space (int nleaves, int count, int nspr, unsigned long long seed)
{
  int i, j;
  char label[40];
  topology current, copy;
  topology_space space = new_topology_space ();

  space->taxlabel = new_char_vector (nleaves);
  space->taxlabel_hash = new_hashtable (nleaves);
  for (i = 0; i < nleaves; i++) {
    snprintf (label, sizeof (label), "T%d", i + 1);
    char_vector_add_string_at_position (space->taxlabel, label, i);
    insert_hashtable (space->taxlabel_hash, space->taxlabel->string[i], i);
  }
  biomcmc_random_number_init (seed);
  current = new_topology (nleaves);
  randomize_topology (current);
  for (i = 0; i < count; i++) {
    if (i) for (j = 0; j < nspr; j++) topology_apply_spr_unrooted (current, false);
    copy = copy_topology_for_calculation (current);
    copy->taxlabel = space->taxlabel;
    space->taxlabel->ref_counter++;
    topology_malloc_blength (copy);
    for (j = 0; j < copy->nnodes; j++) copy->blength[j] = 1.;
    copy->blength[copy->root->id] = 0.;
    copy->blength[copy->root->left->id] = copy->blength[copy->root->right->id] = 0.5;
    add_topology_to_topology_space_if_distinct (copy, space, true);
  }
  space->has_branch_lengths = true;
  del_topology (current);
  biomcmc_random_number_finalize ();
  return space;
}

topology_space
genefam_module_parse_newick_trees (const char *treelist_string, bool use_root_location)
{
  size_t str_size;
  char *str_start, *str_end;
  topology_space ts = NULL;

  str_start = strchr (treelist_string, '(');
  if (!str_start) biomcmc_error ("Couldn't read first tree (single or from list) in newick format");
  while (str_start) {
    str_end = strchr (str_start, ';');
    if (str_end == NULL) str_size = strlen (str_start);
    else str_size = str_end - str_start;
    add_string_with_size_to_topology_space (&ts, str_start, str_size, use_root_location);
    if (str_end) str_start = strchr (str_end, '(');
    else str_start = NULL;
  }
  return ts;
}

static topology
copy_topology_for_calculation (topology source)
{
  topology copy = new_topology (source->nleaves);
  copy_topology_from_topology (copy, source);
  if (source->taxlabel) {
    copy->taxlabel = source->taxlabel;
    copy->taxlabel->ref_counter++;
  }
  return copy;
}

static topology_space
new_calculation_tree_space (topology source)
{
  topology_space working = new_topology_space ();
  working->ndistinct = working->ntrees = 1;
  working->distinct = (topology*) biomcmc_malloc (sizeof (topology));
  working->tree = (topology*) biomcmc_malloc (sizeof (topology));
  working->distinct[0] = copy_topology_for_calculation (source);
  working->tree[0] = working->distinct[0];
  return working;
}

int
genefam_module_treesignal_from_topology_spaces (topology_space genetree, topology_space sptree, double **output_distances)
{
  int n_output = 0;
  topology_space working_gene;

  if ((!genetree) || (!sptree) || (genetree->ndistinct < 1) || (sptree->ndistinct < 1)) return 0;
  working_gene = new_calculation_tree_space (genetree->distinct[0]);
  generate_output_distances (working_gene, sptree, output_distances, &n_output, false);
  del_topology_space (working_gene);
  return n_output;
}

int
genefam_module_treesignal_from_topology_spaces_rescale (topology_space genetree, topology_space sptree, double **output_distances)
{
  int n_output = 0;
  topology_space working_gene;

  if ((!genetree) || (!sptree) || (genetree->ndistinct < 1) || (sptree->ndistinct < 1)) return 0;
  working_gene = new_calculation_tree_space (genetree->distinct[0]);
  generate_output_distances (working_gene, sptree, output_distances, &n_output, true);
  del_topology_space (working_gene);
  return n_output;
}

int
genefam_module_treesignal_fromtrees (const char *gtree_str, const char *splist_str, double **output_distances)
{
  int n_output = 0;
  topology_space genetree = NULL, sptree = NULL;

  genetree = genefam_module_parse_newick_trees (gtree_str, false);
  sptree = genefam_module_parse_newick_trees (splist_str, true);
  // TODO: reduce sptrees to genetree species here (and notice that afterwards some sptrees may be same)
  /* 
  int i;
  char *s; s = topology_to_string_by_id (genetree->distinct[0], NULL);  printf ("DEBUG2 g: %s\n", s); free(s); 
  for (i=0; i < sptree->ndistinct; i++) {
    s = topology_to_string_by_id (sptree->distinct[i], NULL);    printf ("DEBUG2 s[%d]: %s\n", i, s); free(s); 
  } */

  n_output = genefam_module_treesignal_from_topology_spaces (genetree, sptree, output_distances);
  del_topology_space (genetree);
  del_topology_space (sptree);
  return n_output;
}

int
genefam_module_treesignal_fromtrees_rescale (const char *gtree_str, const char *splist_str, double **output_distances)
{
  int n_output = 0;
  topology_space genetree = NULL, sptree = NULL;


  sptree = genefam_module_parse_newick_trees (splist_str, true);
  genetree = genefam_module_parse_newick_trees (gtree_str, false);

  n_output = genefam_module_treesignal_from_topology_spaces_rescale (genetree, sptree, output_distances);
  del_topology_space (genetree);
  del_topology_space (sptree);
  return n_output;
}

int
genefam_module_treesignal_from_topology_spaces_pvalue (topology_space genetree, topology_space sptree, int n_reps, double **output_distances)
{
  int n_output = 0, i, j, k, *pvalues, n_samples = 0;
  double *obs_distances = NULL;
  double dist_bounds[2 * NDISTS]; // min and max values (first and last values, respect.)
  bool valid_trees = true;
  topology_space working_gene;
  topology working_species;

  if ((!genetree) || (!sptree) || (genetree->ndistinct < 1) || (sptree->ndistinct < 1)) return 0;
  working_gene = new_calculation_tree_space (genetree->distinct[0]);
  working_species = copy_topology_for_calculation (sptree->distinct[0]);

  valid_trees = generate_output_distances (working_gene, sptree, &obs_distances, &n_output, false);
  *output_distances = (double*) biomcmc_malloc (sizeof (double) * 2 * n_output); /* pointer used (and freed) by calling function */
  if (! valid_trees) { // exit now, returning bogus vector 
    for (i=0; i < 2 * n_output; i++) (*output_distances)[i] = -1.;
    if (obs_distances) free (obs_distances);
    del_topology (working_species);
    del_topology_space (working_gene);
    return 2 * n_output;
  }

  pvalues = (int*) biomcmc_malloc (sizeof (int) * n_output); /* for openMP this should be a 2D matrix */
  biomcmc_random_number_init(0ULL);

  for (i=0; i < n_output; i++) pvalues[i] = 0;
  for (i=0; i < sptree->ndistinct; i++) for (j=0; j < sptree->ndistinct; j++) for (k=0; k < NDISTS; k++) 
    if ((i != j ) && (obs_distances[NDISTS*j + k] >= obs_distances[NDISTS*i + k])) pvalues[NDISTS*j + k]++;

  for (k=0; k < NDISTS; k++) dist_bounds[k] = dist_bounds[k + NDISTS] = obs_distances[k]; // initial values for min and max
  for (i=1; i < sptree->ndistinct; i++) for (k=0; k < NDISTS; k++) {
    if (dist_bounds[k] > obs_distances[NDISTS*i + k]) dist_bounds[k] = obs_distances[NDISTS*i + k]; // min 
    if (dist_bounds[NDISTS+k] < obs_distances[NDISTS*i + k]) dist_bounds[NDISTS+k] = obs_distances[NDISTS*i + k]; // max
  }

  n_samples = sptree->ndistinct - 1; /* number of comparisons per tree, for calculating p-value */
  n_samples += update_randomized_distances (working_gene->distinct[0], working_species, sptree->ndistinct, obs_distances, pvalues, n_reps, dist_bounds);
  for (i=0; i < n_output; i++) (*output_distances)[i] = ((double)pvalues[i])/(double)(n_samples); /* first trees X NDISTS are p-values */

  /* second half of output vector has tree-distance pairs, scaled to min and max observed values */
  for (k=0; k < NDISTS; k++) if (dist_bounds[NDISTS+k] - dist_bounds[k] < 0.00001) dist_bounds[k] = -1.; // min==max, then make all values = 1
  //for (k=0; k < NDISTS; k++) { printf ("%3.0lf %4.0lf | ", dist_bounds[k], dist_bounds[NDISTS+k]); } printf (" bounds::DEBUG \n"); 
  for (i=0; i < sptree->ndistinct; i++) for (k=0; k < NDISTS; k++) // y = (x-min)/(max-min)
    (*output_distances)[n_output + NDISTS * i + k] = (obs_distances[NDISTS * i + k] - dist_bounds[k])/(dist_bounds[NDISTS+k] - dist_bounds[k]);

  del_topology (working_species);
  del_topology_space (working_gene);
  if (pvalues) free (pvalues);
  if (obs_distances) free (obs_distances);
  biomcmc_random_number_finalize(); /* free the global variable */
  return 2 * n_output;
}

int
genefam_module_treesignal_fromtrees_pvalue (const char *gtree_str, const char *splist_str, int n_reps, double **output_distances)
{
  int n_output;
  topology_space genetree = genefam_module_parse_newick_trees (gtree_str, false);
  topology_space sptree = genefam_module_parse_newick_trees (splist_str, true);

  n_output = genefam_module_treesignal_from_topology_spaces_pvalue (genetree, sptree, n_reps, output_distances);
  del_topology_space (genetree);
  del_topology_space (sptree);
  return n_output;
}

bool
generate_output_distances (topology_space gtree, topology_space stree, double **distances, int *n, bool rescale)
{
  int i, j=0;
  double maxdistance[NDISTS];
  splitset split;

  *n = stree->ndistinct * NDISTS;
  (*distances) = (double*) biomcmc_malloc (sizeof (double) * *n);

  split = create_splitset_dSPR_genespecies (gtree->distinct[0], stree->distinct[0]);
  init_tree_recon_from_species_topology (gtree->distinct[0], stree->distinct[0]);
  if (! calculate_max_distances (gtree->distinct[0], maxdistance, split)) { // one of the trees is too small 
    for (j = 0; j < *n; j++) (*distances)[j] = -1.;
    fprintf (stderr, "WARNING: one of the trees is too small, or they do not have enough species in common; setting whole spectrum to '-1'\n");
    del_splitset (split);
    return false;
  }
  //for (j = 0; j < NDISTS; j++) { printf ("%.1lf | ", maxdistance[j]); } printf (" max::DEBUG \n"); 

  for (i=0; i < stree->ndistinct; i++) {
    //printf ("%d ", i); fflush(stdout);  // ::DEBUG::
    gene_tree_reconcile (gtree->distinct[0], stree->distinct[i]);
    dSPR_gene_species (gtree->distinct[0], stree->distinct[i], split);
    (*distances)[NDISTS * i + 0] = (double) gtree->distinct[0]->rec->ndups; 
    (*distances)[NDISTS * i + 1] = (double) gtree->distinct[0]->rec->nloss; 
    (*distances)[NDISTS * i + 2] = (double) gtree->distinct[0]->rec->ndcos;
    (*distances)[NDISTS * i + 3] = (double) split->rf;
    (*distances)[NDISTS * i + 4] = (double) split->hdist;
    (*distances)[NDISTS * i + 5] = (double) split->spr + split->spr_extra;   // order changed in 2018.08.15 (used to be spr, rf, hdist)
    if (rescale) for (j = 0; j < NDISTS; j++) (*distances)[NDISTS * i + j] /= maxdistance[j]; 
  }
  //printf ("sptrees::DEBUG \n"); 
  del_splitset (split);
  // TODO: If sptrees are repeated, we must transform distances[distinct] into distances_new[tree] and then re-map later to avoid repeated computations
  return true;
}

int
update_randomized_distances (topology gtree, topology stree, int n_obs, double *obs_dist, int *pvalues, int n_reps, double *bounds)
{
  int i, j, k;
  splitset split;
  double tmp_dist[NDISTS];

  split = create_splitset_dSPR_genespecies (gtree, stree);
  init_tree_recon_from_species_topology (gtree, stree);

  for (i=0; i < n_reps; i++) {
    randomize_topology (gtree); 
    randomize_topology (stree); 
    gene_tree_reconcile (gtree, stree);
    dSPR_gene_species (gtree, stree, split);
    tmp_dist[0] = (double) gtree->rec->ndups; 
    tmp_dist[1] = (double) gtree->rec->nloss; 
    tmp_dist[2] = (double) gtree->rec->ndcos;
    tmp_dist[3] = (double) split->rf;
    tmp_dist[4] = (double) split->hdist;
    tmp_dist[5] = (double) split->spr + split->spr_extra;
    /* pvalues[] is freq of distances as low as observed */
    for (j=0; j < n_obs; j++) for (k=0; k < NDISTS; k++) if (obs_dist[NDISTS*j + k] >= tmp_dist[k]) pvalues[NDISTS*j + k]++;
    for (k=0; k < NDISTS; k++) {
      if (bounds[k]        > tmp_dist[k]) bounds[k]        = tmp_dist[k]; 
      if (bounds[NDISTS+k] < tmp_dist[k]) bounds[NDISTS+k] = tmp_dist[k];
    }
  }

  del_splitset (split);
  return n_reps;
}

bool
calculate_max_distances (topology gt, double *maxd, splitset split)
{
  if ((gt->nleaves < 4) || (gt->rec->sp_size < 4) || (split->n_g < 2) || (split->n_s < 2)) return false;
  // max dups: all gene nodes (n-1) point to same (root) sptree node = n-2 dups 
  maxd[0] = (double) (gt->nleaves - 2);
  // max losses: for each dup, (n-1) losses leaving only one leaf = n_dups*(n-1), where n is number of species
  maxd[1] = (double) ((gt->nleaves - 2) * gt->rec->sp_size);
  // max deepcoal: dcos=loss-2 dups + 2 |diff in gt and st|; this bound is based on simulations as well
  maxd[2] = (double) ((gt->nleaves - 2) * (gt->rec->sp_size - 1) + (2 * gt->rec->size_diff));
  // max mulRF: total number of bipartitions
  maxd[3] = (double) (split->n_g + split->n_s);
  // max hdist = when disagreement large (n/2) for all edge pairs, w/ n=gene tree leaves
  maxd[4] = (double) (split->n_g * split->n_g);
  // max uSPR: n-3
  maxd[5] = (double) (((split->n_g) < (split->n_s)) ? (split->n_g) : (split->n_s));
  return true;
}

char*
genefam_module_randomise_trees_with_spr (const char *splist_str, int n_copies, int n_spr)
{
  topology_space strees = NULL;
  topology randtree;
  int n_strees, i, j, k;
  size_t  new_str_size;
  char *s, *output_tree_string;

  biomcmc_random_number_init(0ULL);
  strees = genefam_module_parse_newick_trees (splist_str, true);
  n_strees = strees->ndistinct;

  for (i=0; i < n_strees; i++) for (j = 0; j < n_copies; j++) {
    randtree = new_topology (strees->distinct[0]->nleaves);
    copy_topology_from_topology (randtree, strees->distinct[i]);
    randtree->taxlabel = strees->taxlabel; /* taxlabel is shared among all topologies */
    strees->taxlabel->ref_counter++;    /* since it is shared, it cannot be deleted if still in use */

    for (k = 0; k < n_spr; k++) topology_apply_spr_unrooted (randtree, false); /* SPR that does NOT preserve root location */
    add_topology_to_topology_space_if_distinct (randtree, strees, true); /* 'true' means that different rootings are treated as distinct trees */
  }

  s = topology_to_string_by_name (strees->distinct[0], NULL); /* second parameter is vector with branch lengths */
  new_str_size = strlen (s) + 2; /* semicolon and terminating NUL */
  output_tree_string = (char*) biomcmc_malloc (new_str_size);
  snprintf (output_tree_string, new_str_size, "%s;", s);
  free (s);
  for (i=1; i < strees->ndistinct; i++) {
    size_t current_size = strlen (output_tree_string);
    s = topology_to_string_by_name (strees->distinct[i], NULL);
    new_str_size = current_size + strlen (s) + 3; /* space, semicolon, NUL */
    output_tree_string = (char*) biomcmc_realloc (output_tree_string, new_str_size);
    snprintf (output_tree_string + current_size, new_str_size - current_size, " %s;", s);
    free (s);
  }

  biomcmc_random_number_finalize(); /* free the global variable */
  del_topology_space (strees);

  return output_tree_string;
}

char*
genefam_module_generate_spr_trees (int n_leaves, int n_iter, int n_spr)
{
  topology origtree, randtree; /* these trees don't have taxa labels */
  char *s, *output_tree_string;
  int i, j;
  size_t  new_str_size;
  splitset split;

  if (n_leaves < 4) n_leaves = 4; 
  if (n_iter < 0) n_iter = 0; 
  if (!n_spr) n_spr = n_leaves - 3; 
  biomcmc_random_number_init(0ULL);
  split = new_splitset_dSPR (n_leaves);
  origtree = new_topology (n_leaves);
  randtree = new_topology (n_leaves);

  randomize_topology (randtree);
  s = topology_to_string_create_name (randtree, NULL); /* second parameter is vector with branch lengths */
  new_str_size = strlen (s) + 2; /* semicolon and terminating NUL */
  output_tree_string = (char*) biomcmc_malloc (new_str_size);
  snprintf (output_tree_string, new_str_size, "%s;", s);
  free (s);

  for (i=0; i < n_iter; i++) {
    copy_topology_from_topology (origtree, randtree);
    do {
      for (j = 0; j < n_spr; j++) topology_apply_spr_unrooted (randtree, false);
    } while (topology_is_equal_unrooted (origtree, randtree, split, false));
    /* add random tree to list (string) */
    size_t current_size = strlen (output_tree_string);
    s = topology_to_string_create_name (randtree, NULL);
    new_str_size = current_size + strlen (s) + 3; /* space, semicolon, NUL */
    output_tree_string = (char*) biomcmc_realloc (output_tree_string, new_str_size);
    snprintf (output_tree_string + current_size, new_str_size - current_size, " %s;", s);
    free (s);
  }

  biomcmc_random_number_finalize(); /* free the global variable */
  del_topology (randtree);
  del_topology (origtree);
  del_splitset (split);

  return output_tree_string;
}
