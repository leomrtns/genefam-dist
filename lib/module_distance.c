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

#include "module_distance.h"

/*! \brief calculate distances between gene tree and collection of sptrees, and storing into newly allocated distances[] */ 
void generate_output_distances (topology_space gtree, topology_space stree, double **distances, int *n);

int
genefam_module_treesignal_fromtrees (const char *gtree_str, const char *splist_str, double **output_distances)
{
  int n_output = 0;
  topology_space genetree = NULL, sptree = NULL;
  size_t str_size;
  char *str_start, *str_end;

  str_start = strchr (gtree_str, '(');
  if (!str_start) biomcmc_error ("Couldn't read gene tree in newick format");
  str_end   = strchr (str_start, ';');
  if (str_end == NULL) str_size = strlen (str_start);
  else str_size = str_end - str_start;
  add_string_with_size_to_topology_space (&genetree, str_start, str_size);

  str_start = strchr (splist_str, '(');
  if (!str_start) biomcmc_error ("Couldn't read first species tree in newick format");
  while (str_start) {
    str_end = strchr (str_start, ';');
    if (str_end == NULL) str_size = strlen (str_start); /* function strchrnul() does this, but may not be portable? */
    else str_size = str_end - str_start; 
    add_string_with_size_to_topology_space (&sptree, str_start, str_size);
    if (str_end) str_start = strchr (str_end, '('); 
    else str_start = NULL;
  }
//  char *s; s = topology_to_string_by_id (genetree->distinct[0], NULL);  printf ("DEBUG2 g: %s\n", s); free(s); 
//  s = topology_to_string_by_id (sptree->distinct[0], NULL);    printf ("DEBUG2 s: %s\n", s); free(s); 

  generate_output_distances (genetree, sptree, output_distances, &n_output);
  del_topology_space (genetree);
  del_topology_space (sptree);
  printf("OUTPUT %d\n", n_output);
  return n_output;
}

void
generate_output_distances (topology_space gtree, topology_space stree, double **distances, int *n)
{
  int i;
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
  }
  del_splitset (split);
  return;
}
