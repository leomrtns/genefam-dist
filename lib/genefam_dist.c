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

/* later I can create/delete sptree_ratchet w/ distinct parameters (for python) using same genefam_sptree struct */

genefam_sptree new_genefam_sptree_from_species_char_vector (const char_vector species_names, const int n_sptrees, const int n_ratchet);
void del_genefam_sptree (genefam_sptree gf);
genetree internal_new_genetree_from_gene_and_sp_topologies (topology g_tree, topology s_tree);
void internal_del_genetree (genetree gt);
double internal_calculate_genetree_distance (genetree gt, topology s_tree);

sptree_ratchet internal_new_sptree_ratchet (int n_genes, int n_ratchet, int n_proposal, topology_space sptree); 
void internal_del_sptree_ratchet (sptree_ratchet sr);

/*! \brief create a new topology_space from string with list of trees in newick format */
topology_space internal_topology_space_from_newick_string (const char *treelist_string, bool use_root_location);

genefam_sptree
new_genefam_sptree_from_species_char_vector (const char_vector species_names, const int n_sptrees, const int n_ratchet)
{
  genefam_sptree gf;
  int i;

  gf = (genefam_sptree) biomcmc_malloc (sizeof (struct genefam_sptree_struct));
  gf->n_sptrees = n_sptrees;
  if (gf->n_sptrees < 4) gf->n_sptrees = 4;  
  gf->n_genefams = 0;
  gf->genefam = NULL;
  gf->best_trees = NULL;
  gf->next_sptree = 0;
  gf->sptree = (topology*) biomcmc_malloc (gf->n_sptrees * sizeof (topology));

  char_vector_remove_duplicate_strings (species_names); 
  char_vector_longer_first_order (species_names); // order species names from longest to shortest (speed up gene/spnames comparison) 
  for (i = 0; i < gf->n_sptrees; i++) {
    gf->sptree[i] = new_topology (species_names->nstrings);
    gf->sptree[i]->taxlabel = species_names; gf->sptree[i]->taxlabel->ref_counter++; 
    randomize_topology (gf->sptree[i]);
  } 
  return gf;
}

void
del_genefam_sptree (genefam_sptree gf)
{
  int i;
  if (!gf) return;
  if (gf->sptree) {
    for (i = gf->n_sptrees - 1; i >= 0 ; i--) free_topology (gf->sptree[i]);
    free (gf->sptree);
  }
  if (gf->genefam) {
    for (i = gf->n_genefams - 1; i >= 0 ; i--) internal_del_genetree(gf->genefam[i]);
    free (gf->genefam);
  }
  internal_del_sptree_ratchet (gf->best_trees);
  free (gf);
}

int
add_genefam_from_string()
{ // FIXME (just copied from find_best_tree.c)
  // must check if possible to store idx into reconciliation 
  biomcmc_realloc ((int*) idx_gene_to_sp, genetre->distinct[0]->nleaves * sizeof (int));
  index_sptaxa_to_genetaxa (sptaxa, genetre->taxlabel, idx_gene_to_sp, ef);/* map species names to gene names and store into idx[] */
  valid_species_size = prepare_spdistmatrix_from_gene_species_map (pool->this_gene_spdist, idx_gene_to_sp, genetre->distinct[0]->nleaves);
  if (valid_species_size > 4) {
    char_vector_add_string (gfilename, arg_filename[i]);
    update_pooled_matrix_from_gene_tree (pool, genetre->distinct[0], idx_gene_to_sp);
  }
  del_topology_space (genetre);
  return 1; // success
}

int
find_initial_best_trees_from_genefam ()
{ // this function should use pool matrix, and create the sptree_ratchet best_trees
}

genetree 
internal_new_genetree_from_gene_and_sp_topologies (topology g_tree, topology s_tree)
{
  int i;
  genetree gt;
  gt = (genetree) biomcmc_malloc (sizeof (struct genetree_struct));
  gt->tree = g_tree; gt->tree->ref_counter++; // point to element of topology_space (usually)
  gt->split = create_splitset_dSPR_genespecies (gt->tree, s_tree);
  // init_tree_recon_from_species_topology() would order sptree names everytime (and already calculates distances...)
  if (!gt->tree->rec) gt->tree->rec = new_reconciliation (g_tree->nleaves, s_tree->nleaves); // if gene has rec, it better be from same sptree size!
  index_sptaxa_to_reconciliation (s_tree->taxlabel, g_tree->taxlabel, gt->tree->rec);
  for (i = 0; i < NDISTS; i++) {
    gt->minmax[i] = 1.e35;  // higher bound for min
    gt->minmax[i+NDISTS] = -1.e35; // lower bound for max
    gt->dist[i] = -1.;
  }
  del_char_vector (g_tree->taxlabel); // should just decrease ref_counter and then delete later, with topology_space
  for (i = 0; i < 32; i++) {
    randomize_topology (s_tree); 
    internal_calculate_genetree_distance (gt, s_tree);
  }
  return gt;
}

void
internal_del_genetree (genetree gt)
{
  if (!gt) return;
  del_splitset (gt->split);
  del_topology (gt->tree);
  free (gt);
  return;
}

double
internal_calculate_genetree_distance (genetree gt, topology s_tree) 
{ /* local_optima -> return zero if a new local optimum is found (new minimum for _some_ distance)  */
  int k;
  double g_score = 0.;
  gene_tree_reconcile (gt->tree, s_tree);
  dSPR_gene_species (gt->tree, s_tree, gt->split); 
  gt->dist[0] = (double) gt->tree->rec->ndups;
  gt->dist[1] = (double) gt->tree->rec->nloss;
  gt->dist[2] = (double) gt->tree->rec->ndcos;
  gt->dist[3] = (double) gt->split->rf;
  gt->dist[4] = (double) gt->split->hdist;
  gt->dist[5] = (double) gt->split->spr + gt->split->spr_extra;  // order same as guenomu, NOT same as genefam_dist or treesignal.py
  /* minimum values are always updated */
  for (k=0; k < NDISTS; k++) {
    if (gt->minmax[k] > gt->dist[k]) { gt->minmax[k] = gt->dist[k]; g_score = -1; } // min
    if (gt->minmax[k+NDISTS] < gt->dist[k])  gt->minmax[k+NDISTS] = gt->dist[k]; // max
  }
  if (g_score < 0) return 0.; // if new minimum is found; however if sptree has been seen before then regular score is calculated
  g_score = 0.;
  for (k=0; k < NDISTS; k++) g_score += biomcmc_log1p( (double)(gt->dist[k] - gt->minmax[k])/(double)(gt->minmax[k+NDISTS] - gt->minmax[k]) ); 
  return g_score;
}

sptree_ratchet
internal_new_sptree_ratchet (int n_genes, int n_ratchet, int n_proposal, topology_space sptree)
{ // n_minibatches will be given by optimiser function
  int i;
  sptree_ratchet sr; 
  sr = (sptree_ratchet) biomcmc_malloc (sizeof (struct sptree_ratchet_struct));
  sr->n_genesamples = n_genes;
  sr->n_ratchet = n_ratchet;
  sr->n_proposal = n_proposal;
  sr->best_score = 1.e35;
  if (sr->n_genesamples < 2) sr->n_genesamples = 2;
  if (sr->n_ratchet < 8) sr->n_ratchet = 8; // ratchet[0] has best score at first (due to initial_sorting() )
  if (sr->n_proposal < 2) sr->n_proposal = 2;
  sr->next_avail = sr->n_ratchet - 1; // idx of currently worse tree (which is just before best tree in a ratchet)

  sr->genesample = (genetree*) biomcmc_malloc (sr->n_genesamples * sizeof (genetree)); 
  sr->ratchet  = (topology*) biomcmc_malloc (sr->n_ratchet * sizeof (topology)); 
  sr->proposal = (topology*) biomcmc_malloc (sr->n_proposal * sizeof (topology));
  sr->ratchet_score  = (double*) biomcmc_malloc (sr->n_ratchet * sizeof (double)); 

  for (i=0; i < sr->n_ratchet; i++) {
    sr->ratchet[i] = init_topol_for_new_gene_sptrees (sptree->distinct[0]);
    randomize_topology (sr->ratchet[i]); // should be filled outside this function (to recycle optimal sptrees from previous iteration)
  }
  for (i=0; i < sr->n_proposal; i++) sr->proposal[i] = init_topol_for_new_gene_sptrees (sptree->distinct[0]);
  for (i=0; i < sr->n_genes; i++) sr->genesample[i] = NULL; // must be created later, when we receive the gene topol_spaces 

  initialise_ratchet_from_topol_space (sr, sptree);
  return sr;
}

void
internal_del_sptree_ratchet (sptree_ratchet sr)
{
  int i;
  if (!sr) return;
  if (sr->genesample) {
    for (i = sr->n_genesamples - 1; i >= 0 ; i--) internal_del_genetree (sr->genesample[i]);
    free (sr->genesample);
  }
  if (sr->ratchet) {
    for (i = sr->n_ratchet - 1; i >= 0 ; i--) del_topology (sr->ratchet[i]);
    free (sr->ratchet);
  }
  if (sr->proposal) {
    for (i = sr->n_proposal - 1; i >= 0 ; i--) del_topology (sr->proposal[i]);
    free (sr->proposal);
  }
  if (sr->ratchet_score) free (sr->ratchet_score);
  free (sr);
  return;
}

topology_space
internal_topology_space_from_newick_string (const char *treelist_string, bool use_root_location)
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
    add_string_with_size_to_topology_space (&ts, str_start, str_size, use_root_location);
    if (str_end) str_start = strchr (str_end, '('); 
    else str_start = NULL;
  }
  return ts;
}

/*** First implementation of treesignal relied on external species trees, and didn't have access to C structs 
 *   Legacy code below 
 ***/

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

  genetree = internal_topology_space_from_newick_string (gtree_str, false); /* bool -> use_root_location */
  sptree = internal_topology_space_from_newick_string (splist_str, true);
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


  sptree = internal_topology_space_from_newick_string (splist_str, true);
  genetree = internal_topology_space_from_newick_string (gtree_str, false);

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

  genetree = internal_topology_space_from_newick_string (gtree_str, false);
  sptree = internal_topology_space_from_newick_string (splist_str, true);

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
  char *s, *output_tree_string, *tmp_string;

  biomcmc_random_number_init(0ULL);
  strees = internal_topology_space_from_newick_string (splist_str, true);
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

