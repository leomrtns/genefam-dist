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

/*! \brief create a new topology_space from string with list of trees in newick format */
topology_space topology_space_from_newick_string (const char *treelist_string);
/*! \brief calculate distances between gene tree and collection of sptrees, and storing into newly allocated distances[] */ 
void generate_output_distances (topology_space gtree, topology_space stree, double **distances, int *n, bool rescale);
/*! \brief calculate distances between gene tree and random sptrees, accumulating into previously allocated pvalues[] */ 
int update_randomized_distances (topology gtree, topology stree, int n_obs, double *obs_dist, int *pvalues, int n_reps);

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
  topology_space genetree = NULL, sptree = NULL;

  genetree = topology_space_from_newick_string (gtree_str);
  sptree = topology_space_from_newick_string (splist_str);
  biomcmc_random_number_init(0ULL);

  generate_output_distances (genetree, sptree, &obs_distances, &n_output, false);
  pvalues = (int*) biomcmc_malloc (sizeof (int) * n_output); /* for openMP this should be a 2D matrix */
  for (i=0; i < n_output; i++) pvalues[i] = 0;
  for (i=0; i < sptree->ndistinct; i++) for (j=0; j < sptree->ndistinct; j++) for (k=0; k < 7; k++) 
    if ((i != j ) && (obs_distances[7*j + k] >= obs_distances[7*i + k])) pvalues[7*j + k]++;
  n_samples = sptree->ndistinct - 1; /* number of comparisons per tree, for calculating p-value */

  n_samples += update_randomized_distances (genetree->distinct[0], sptree->distinct[0], sptree->ndistinct, obs_distances, pvalues, n_reps);

  *output_distances = (double*) biomcmc_malloc (sizeof (double) * n_output); /* pointer used (and freed) by calling function */
  for (i=0; i < n_output; i++) (*output_distances)[i] = ((double)pvalues[i])/(double)(n_samples);

  del_topology_space (genetree);
  del_topology_space (sptree);
  if (pvalues) free (pvalues);
  if (obs_distances) free (obs_distances);
  biomcmc_random_number_finalize(); /* free the global variable */
  return n_output;
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

void
generate_output_distances (topology_space gtree, topology_space stree, double **distances, int *n, bool rescale)
{
  int i, j=0;
  splitset split;

  *n = stree->ndistinct * 7;
  (*distances) = (double*) biomcmc_malloc (sizeof (double) * *n);
  split = create_splitset_dSPR_genespecies (gtree->distinct[0], stree->distinct[0]);
  for (i=0; i < stree->ndistinct; i++) {
    init_tree_recon_from_species_topology (gtree->distinct[0], stree->distinct[i]);
    dSPR_gene_species (gtree->distinct[0], stree->distinct[i], split);
    (*distances)[7 * i + 0] = (double) gtree->distinct[0]->rec->ndups; 
    (*distances)[7 * i + 1] = (double) gtree->distinct[0]->rec->nloss; 
    (*distances)[7 * i + 2] = (double) gtree->distinct[0]->rec->ndcos;
    (*distances)[7 * i + 3] = (double) split->spr;
    (*distances)[7 * i + 4] = (double) split->spr_extra;
    (*distances)[7 * i + 5] = (double) split->rf;
    (*distances)[7 * i + 6] = (double) split->hdist;
    if (rescale) { // rec->sp_size is the effective number of species
      double maxvalue = 2. * (double)(gtree->distinct[0]->rec->sp_size + gtree->distinct[0]->nleaves) - 3.;
      for (j = 0; j < 7; j++) (*distances)[7 * i + j] /= maxvalue; 
    }
  }
  del_splitset (split);
  // TODO: If sptrees are repeated, we must transform distances[distinct] into distances_new[tree]
  return;
}

int
update_randomized_distances (topology gtree, topology stree, int n_obs, double *obs_dist, int *pvalues, int n_reps)
{
  int i, j, k;
  splitset split;
  double *tmp_dist;

  split = create_splitset_dSPR_genespecies (gtree, stree);
  init_tree_recon_from_species_topology (gtree, stree);
  tmp_dist = (double*) biomcmc_malloc (sizeof (double) * 7);

  for (i=0; i < n_reps; i++) {
    randomize_topology (stree); 
    dSPR_gene_species (gtree, stree, split);
    tmp_dist[0] = (double) gtree->rec->ndups; 
    tmp_dist[1] = (double) gtree->rec->nloss; 
    tmp_dist[2] = (double) gtree->rec->ndcos;
    tmp_dist[3] = (double) split->spr;
    tmp_dist[4] = (double) split->spr_extra;
    tmp_dist[5] = (double) split->rf;
    tmp_dist[6] = (double) split->hdist;
    /* pvalues[] is freq of distances as low as observed */
    for (j=0; j < n_obs; j++) for (k=0; k < 7; k++) if (obs_dist[7*j + k] >= tmp_dist[k]) pvalues[7*j + k]++;
  }

  if (tmp_dist) free (tmp_dist);
  del_splitset (split);
  return n_reps;
}

topology_space
topology_space_from_newick_string (const char *treelist_string)
{
  size_t str_size;
  char *str_start, *str_end;
  topology_space ts = NULL;

  str_start = strchr (treelist_string, '(');
  if (!str_start) biomcmc_error ("Couldn't read first tree (from potential list) in newick format");
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
