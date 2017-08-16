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

/*! \brief create a new topology_space from string with list of trees in newick format */
topology_space topology_space_from_newick_string (const char *treelist_string);
/*! \brief calculate distances between gene tree and collection of sptrees, and storing into newly allocated distances[] */ 
bool generate_output_distances (topology_space gtree, topology_space stree, double **distances, int *n, bool rescale);
/*! \brief calculate distances between gene tree and random sptrees, accumulating into previously allocated pvalues[] */ 
int update_randomized_distances (topology gtree, topology stree, int n_obs, double *obs_dist, int *pvalues, int n_reps, double *maxdist);
/*! \brief defines upper bound for all distances based on theory, not empirical as p-value above */
bool calculate_max_distances (topology gt, double *maxd, splitset split);

int
genefam_module_treesignal_fromtrees (const char *gtree_str, const char *splist_str, double **output_distances)
{
  int n_output = 0;
  topology_space genetree = NULL, sptree = NULL;

  genetree = topology_space_from_newick_string (gtree_str);
  sptree = topology_space_from_newick_string (splist_str);
  // TODO: reduce sptrees to genetree species here (and notice that afterwards some sptrees may be same)
  /* 
  int i;
  char *s; s = topology_to_string_by_id (genetree->distinct[0], NULL);  printf ("DEBUG2 g: %s\n", s); free(s); 
  for (i=0; i < sptree->ndistinct; i++) {
    s = topology_to_string_by_id (sptree->distinct[i], NULL);    printf ("DEBUG2 s[%d]: %s\n", i, s); free(s); 
  } */

  generate_output_distances (genetree, sptree, output_distances, &n_output, false);
  del_topology_space (genetree);
  del_topology_space (sptree);
  return n_output;
}

int
genefam_module_treesignal_fromtrees_rescale (const char *gtree_str, const char *splist_str, double **output_distances)
{
  int n_output = 0;
  topology_space genetree = NULL, sptree = NULL;

  //printf("genetree::DEBUG %s", gtree_str);

  genetree = topology_space_from_newick_string (gtree_str);
  sptree = topology_space_from_newick_string (splist_str);

  generate_output_distances (genetree, sptree, output_distances, &n_output, true);
  del_topology_space (genetree);
  del_topology_space (sptree);
  return n_output;
}

int
genefam_module_treesignal_fromtrees_pvalue (const char *gtree_str, const char *splist_str, int n_reps, double **output_distances)
{
  int n_output = 0, i, j, k, *pvalues, n_samples = 0;
  double *obs_distances;
  double dist_bounds[2 * NDISTS]; // min and max values (first and last values, respect.)
  bool valid_trees = true;
  topology_space genetree = NULL, sptree = NULL;

  genetree = topology_space_from_newick_string (gtree_str);
  sptree = topology_space_from_newick_string (splist_str);

  valid_trees = generate_output_distances (genetree, sptree, &obs_distances, &n_output, false);
  *output_distances = (double*) biomcmc_malloc (sizeof (double) * 2 * n_output); /* pointer used (and freed) by calling function */
  if (! valid_trees) { // exit now, returning bogus vector 
    for (i=0; i < 2 * n_output; i++) (*output_distances)[i] = -1.;
    if (obs_distances) free (obs_distances);
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
  n_samples += update_randomized_distances (genetree->distinct[0], sptree->distinct[0], sptree->ndistinct, obs_distances, pvalues, n_reps, dist_bounds);
  for (i=0; i < n_output; i++) (*output_distances)[i] = ((double)pvalues[i])/(double)(n_samples); /* first trees X NDISTS are p-values */

  /* second half of output vector has tree-distance pairs, scaled to min and max observed values */
  for (k=0; k < NDISTS; k++) if (dist_bounds[NDISTS+k] - dist_bounds[k] < 0.00001) dist_bounds[k] = -1.; // min==max, then make all values = 1
  //for (k=0; k < NDISTS; k++) { printf ("%3.0lf %4.0lf | ", dist_bounds[k], dist_bounds[NDISTS+k]); } printf (" bounds::DEBUG \n"); 
  for (i=0; i < sptree->ndistinct; i++) for (k=0; k < NDISTS; k++) // y = (x-min)/(max-min)
    (*output_distances)[n_output + NDISTS * i + k] = (obs_distances[NDISTS * i + k] - dist_bounds[k])/(dist_bounds[NDISTS+k] - dist_bounds[k]);

  del_topology_space (genetree);
  del_topology_space (sptree);
  if (pvalues) free (pvalues);
  if (obs_distances) free (obs_distances);
  biomcmc_random_number_finalize(); /* free the global variable */
  return 2 * n_output;
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
    (*distances)[NDISTS * i + 3] = (double) split->spr + split->spr_extra;
    (*distances)[NDISTS * i + 4] = (double) split->rf;
    (*distances)[NDISTS * i + 5] = (double) split->hdist;
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
    tmp_dist[3] = (double) split->spr + split->spr_extra;
    tmp_dist[4] = (double) split->rf;
    tmp_dist[5] = (double) split->hdist;
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
  // max uSPR: n-3
  maxd[3] = (double) (((split->n_g) < (split->n_s)) ? (split->n_g) : (split->n_s));
  // max mulRF: total number of bipartitions
  maxd[4] = (double) (split->n_g + split->n_s);
  // max hdist = when disagreement large (n/2) for all edge pairs, w/ n=gene tree leaves
  maxd[5] = (double) (split->n_g * split->n_g);
  return true;
}

char*
genefam_module_randomise_trees_with_spr (const char *splist_str, int n_copies, int n_spr)
{
  topology_space strees = NULL;
  topology randtree;
  int n_strees, i, j, k;
  size_t  new_str_size;
  char *s, *output_tree_string, *tmp_string;

  biomcmc_random_number_init(0ULL);
  strees = topology_space_from_newick_string (splist_str);
  n_strees = strees->ndistinct;

  for (i=0; i < n_strees; i++) for (j = 0; j < n_copies; j++) {
    randtree = new_topology (strees->distinct[0]->nleaves);
    copy_topology_from_topology (randtree, strees->distinct[i]);
    randtree->taxlabel = strees->taxlabel; /* taxlabel is shared among all topologies */
    strees->taxlabel->ref_counter++;    /* since it is shared, it cannot be deleted if still in use */

    for (k = 0; k < n_spr; k++) topology_apply_spr_unrooted (randtree, false);
    add_topology_to_topology_space_if_distinct (randtree, strees);
  }

  s = topology_to_string_by_name (strees->distinct[0], NULL); /* second parameter is vector with branch lengths */
  new_str_size = strlen (s);
  output_tree_string = (char*) biomcmc_malloc (sizeof(char) * strlen(s));
  sprintf (output_tree_string, "%s;", s); free (s); /* s[] is recycled (malloc'ed again when new topol string is created) */
  for (i=1; i < strees->ndistinct; i++) {
    new_str_size = strrchr (output_tree_string, ';') - output_tree_string + 1; // size of output_tree_string[] with extra '\0'
    tmp_string = (char*) biomcmc_malloc (sizeof(char) * new_str_size + 1);
    strncpy (tmp_string, output_tree_string, new_str_size); /* adds '\0' only if output_tree_string is smaller!! */
    tmp_string[new_str_size] = '\0'; /* one after the semicolon */
    s = topology_to_string_by_name (strees->distinct[i], NULL); /* second parameter is vector with branch lengths */
    new_str_size += strlen (s) + 3;
    output_tree_string = (char*) biomcmc_realloc ((char*) output_tree_string, sizeof(char) * new_str_size);
    sprintf (output_tree_string, "%s %s;", tmp_string, s);  if (s) free (s); if (tmp_string) free (tmp_string); 
  }

  biomcmc_random_number_finalize(); /* free the global variable */
  del_topology (randtree);
  del_topology_space (strees);

  return output_tree_string;
}

char*
genefam_module_generate_spr_trees (int n_leaves, int n_iter, int n_spr)
{
  topology origtree, randtree; /* these trees don't have taxa labels */
  char *s, *output_tree_string, *tmp_string;
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
  new_str_size = strlen (s);
  output_tree_string = (char*) biomcmc_malloc (sizeof(char) * new_str_size);
  sprintf (output_tree_string, "%s;", s); free (s); /* s[] is recycled (malloc'ed again when new topol string is created) */

  for (i=0; i < n_iter; i++) {
    copy_topology_from_topology (origtree, randtree);
    do {
      for (j = 0; j < n_spr; j++) topology_apply_spr_unrooted (randtree, false);
    } while (topology_is_equal_unrooted (origtree, randtree, split, false));
    /* temporary string with copy of current trees */
    new_str_size = strrchr (output_tree_string, ';') - output_tree_string + 1;
    tmp_string = (char*) biomcmc_malloc (sizeof(char) * new_str_size + 1);
    strncpy (tmp_string, output_tree_string, new_str_size); /* adds '\0' only if output_tree_string is smaller!! */
    tmp_string[new_str_size] = '\0'; /* one after the semicolon */
    /* add random tree to list (string) */
    s = topology_to_string_create_name (randtree, NULL); /* second parameter is vector with branch lengths */
    new_str_size += strlen (s) + 3;
    output_tree_string = (char*) biomcmc_realloc ((char*) output_tree_string, sizeof(char) * new_str_size);
    sprintf (output_tree_string, "%s %s;", tmp_string, s);  if (s) free (s); /* s[] is recycled (malloc'ed again when new topol string is created) */
    if (tmp_string) free (tmp_string); 
  }

  biomcmc_random_number_finalize(); /* free the global variable */
  del_topology (randtree);
  del_topology (origtree);
  del_splitset (split);

  return output_tree_string;
}

topology_space
topology_space_from_newick_string (const char *treelist_string)
{
  size_t str_size;
  char *str_start, *str_end;
  topology_space ts = NULL;

  str_start = strchr (treelist_string, '(');
  if (!str_start) biomcmc_error ("Couldn't read first tree (single or from list) in newick format");
  while (str_start) {
    str_end = strchr (str_start, ';');
    if (str_end == NULL) str_size = strlen (str_start); /* function strchrnul() does this, but may not be portable? */
    else str_size = str_end - str_start; 
    add_string_with_size_to_topology_space (&ts, str_start, str_size);
    if (str_end) str_start = strchr (str_end, '('); 
    else str_start = NULL;
  }
  return ts;
}
