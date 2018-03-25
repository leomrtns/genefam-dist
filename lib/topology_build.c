/* 
 * This file is part of guenomu, a hierarchical Bayesian procedure to estimate the distribution of species trees based
 * on multi-gene families data.
 * Copyright (C) 2009  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 *
 * Guenomu is free software; you can redistribute it and/or modify it under the terms of the GNU General Public 
 * License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more 
 * details (file "COPYING" or http://www.gnu.org/copyleft/gpl.html).
 */

#include "topology_build.h"

/*! \brief create internal node with given children (children coalesce into parent node) */
void create_parent_node_from_children (topology tree, int parent, int lchild, int rchild);
/*! \brief recursive function that applies spr (almost always an nni) on a subtree */
bool topology_apply_shortspr_weighted_subtree (topology tree, topol_node lca, double *prob, double scale, bool update_done);

void
randomize_topology (topology tree)
{ 
  int i, lchild, rchild, parent = tree->nleaves, *idx = tree->index, n_idx = tree->nleaves;

  /* tree->index is also used by quasi_randomize_topology(), and here we tell it the info was destroyed */
  tree->quasirandom = false;

  for (i=0; i < n_idx; i++) idx[i] = i;

  while (n_idx > 2) { /* draw two nodes to be connected (like star-topology refinement) */
    i = biomcmc_rng_unif_int (n_idx);
    rchild = idx[i];
    idx[i] = idx[--n_idx]; /* avoid replacement */
    i = biomcmc_rng_unif_int (n_idx);
    lchild = idx[i];
    idx[i] = parent; /* new node we are about to create (available to sampling on next iteration) */
    create_parent_node_from_children (tree, parent, lchild, rchild);
    parent++; 
  }
  /* now idx[] has only two elements and "parent" is now (2*ntax - 2) */
  create_parent_node_from_children (tree, parent, idx[0], idx[1]);
  tree->root = tree->nodelist[parent];
  tree->root->sister = tree->root;
  tree->root->up = NULL;

  update_topology_sisters   (tree);
  topology_apply_spr_unrooted (tree, false); // for some reason the above function does not sample uniformly from unrooted
}

void
quasi_randomize_topology (topology tree, int sample_type)
{ 
  int i, lchild, rchild, parent = tree->nleaves, swap;
  int *taxa_idx = tree->index, 
      *r_idx    = tree->index + tree->nleaves, 
      *l_idx    = tree->index + 2 * tree->nleaves,
      *idx      = tree->index + 3 * tree->nleaves; /* vectors */

  if ((!sample_type) || (!tree->quasirandom)) { sample_type = 14; tree->quasirandom = true; } /* starts a new random tree and prepare vectors */

  /* despite the weird numbering (8,2,4,1), pls do not change the order: first initialize taxa_idx, only then shuffle, etc.
   * this numbering is related to the frequency with which we could do these moves, if user calls using 
   * sample_type=1,2,3,4,5... */
  if (sample_type & 8) for (i = 0; i < tree->nleaves; i++) taxa_idx[i] = i; /* initialization */
  if (sample_type & 2) { /* change order of leaves (does not change shape) */
    for (i = tree->nleaves - 1; i > 0; i--) { /* Knuth shuffle: notice that zero is excluded from loop, and unif(0,i) is with i included */
      lchild = biomcmc_rng_unif_int (i + 1); /* unif(0,i-1) -- that is, i excluded -- would be Sattolo's algorithm */
      swap = taxa_idx[lchild]; taxa_idx[lchild] = taxa_idx[i]; taxa_idx[i] = swap;
    }
  }
  if (sample_type & 4) { /* sample a new shape (does not change leaves) */ 
    for (i = 0; i < tree->nleaves - 2; i++) {
      r_idx[i] = biomcmc_rng_unif_int (tree->nleaves - i);
      l_idx[i] = biomcmc_rng_unif_int (tree->nleaves - i - 1);
    }
  }
  if (sample_type & 1) { /* deterministic cycle over choices (may change shape) */
    for (i = 0; i < tree->nleaves - 2; i++) {
      r_idx[i]--; /* should always be between 0 and nleaves - i - 1 */
      if (r_idx[i] < 0) r_idx[i] = tree->nleaves - i - 1;
      l_idx[i]--; /* should always be between 0 and nleaves - i - 2 */
      if (l_idx[i] < 0) l_idx[i] = tree->nleaves - i - 2;
    }
  }

  /* for all sample types we must copy the indexes to the "disposable" vector */
  for (i = 0; i < tree->nleaves; i++) idx[i] = taxa_idx[i]; 

  for (i = 0; i < tree->nleaves - 2; i++) {/* draw two nodes to be connected (like star-topology refinement) */
    rchild          = idx[ r_idx[i] ];
    idx[ r_idx[i] ] = idx[tree->nleaves - i - 1]; /* avoid replacement */
    lchild          = idx[ l_idx[i] ];
    idx[ l_idx[i] ] = parent; /* new node we are about to create (available to sampling on next iteration) */
    create_parent_node_from_children (tree, parent, lchild, rchild);
    parent++; 
  }
  /* now idx[] has only two elements and "parent" is now (2*ntax - 2) */
  create_parent_node_from_children (tree, parent, idx[0], idx[1]);
  tree->root = tree->nodelist[parent];
  tree->root->sister = tree->root;
  tree->root->up = NULL;

  update_topology_sisters   (tree);
  update_topology_traversal (tree);
}

void
upgma_from_distance_matrix (topology tree, distance_matrix dist, bool single_linkage) 
{ /* always upper diagonal (that is, only i < j in d[i][j]) */
  int i, j, idx_i, parent = tree->nleaves, idx_j, min_row, min_col, row, col, idx_col, n_idx = tree->nleaves,
      *idx = tree->index,                            /* indexes in UPGMA */
      *idxtree = tree->index + tree->nleaves,        /* indexes in tree (since have values > nleaves) */
      *min_by_row = tree->index + 2 * tree->nleaves; /* column having min value for each row */
  double *dst_by_row, dst_row, new_dist, *gsize, *height, gs1, gs2;

  /* tree->index is also used by quasi_randomize_topology(), and here we tell it the info was destroyed */
  tree->quasirandom = false;
  
  if (!tree->blength) topology_malloc_blength (tree);

  gsize      = (double *) biomcmc_malloc (n_idx * sizeof (double)); /* number of leaves below node */
  height     = (double *) biomcmc_malloc (n_idx * sizeof (double)); /* distance from node to tips (ultrametric) */ 
  dst_by_row = (double *) biomcmc_malloc (n_idx * sizeof (double)); /* min value itself for each row */

  for (i=0; i < n_idx; i++) { 
    idx[i]     = i; /* index to actual vector element for UPGMA distance matrix */
    idxtree[i] = i; /* index to actual vector element for tree nodes */
    gsize[i]   = 1.; /* group size (# elements below node in UPGMA */
    height[i]  = 0.; /* distance of node from present */
  }

  /* arbitrary initial values for vector of row-wise minimum distances and overall minimum */
  for (i=0; i < dist->size - 1; i++) dst_by_row[i] = 1.e35; 
  /* update vector with minimum distances per row */
  for (j=1; j < n_idx; j++) for (i=0; i < j; i++)
    if (dist->d[i][j] < dst_by_row[i]) { dst_by_row[i] = dist->d[i][j]; min_by_row[i] = j; }

  while (n_idx > 2) { /* draw two nodes to be connected */
    /* find vector elements with overall minimum distance */ 
    dst_row = 1.e35;
    for (i=0; i < n_idx; i++) if ((idx[i] < dist->size - 1) && (dst_by_row[idx[i]] < dst_row)) { 
      dst_row = dst_by_row[idx[i]]; min_row = i; min_col = min_by_row[idx[i]]; 
    }
    if (dst_row < 1.e-35) dst_row = 1.e-35;

    /* choose pair with smallest distance */
    i = min_row;  j = min_col; /* relative to idx[] */
    /* tree node creation */
    idx_i      = idxtree[i];
    idx_j      = idxtree[j];
    idxtree[i] = parent; /* new node we are about to create (available to sampling on next iteration) */
    idxtree[j] = idxtree[--n_idx]; /* avoid replacement (should come last since we may have i == n_idx-1 */
    create_parent_node_from_children (tree, parent, idx_i, idx_j);
    parent++; /* parent on next iteration */

    /* calculate branch lengths (UPGMA distances refer to tip, that's why we keep node UPGMA dists on  height[]) */
    gs1 = dst_row/2. - height[idx[i]]; 
    gs2 = dst_row/2. - height[idx[j]]; 
    if (gs1 < 1.e-35) gs1 = 1.e-35;
    if (gs2 < 1.e-35) gs2 = 1.e-35;
    tree->blength[idx_i] = gs1; /* UPGMA distance */ 
    tree->blength[idx_j] = gs2; /* UPGMA distance */
    height[idx[i]] = dst_row/2.;

    /* idx_j nd idx_j have indexes of actual matrix elements */
    idx_i     = idx[i]; /* unlike index of tree nodes, we don't change idx[i] which will hold dists for internal node */
    idx_j     = idx[j]; 
    idx[j]    = idx[n_idx]; /* avoid replacement */

    /* update distance matrix */
    dst_by_row[idx_i] = 1.e35; /* will need update */
    gs1 = (gsize[idx_i] + gsize[idx_j]);
    for (i=0; i < n_idx; i++) { 
      /* calculates distances to new node */
      if (single_linkage) {  /* a.k.a nearest neighbor clustering. Distance to new node is minimum between elements */
        if (idx[i] < idx_j) new_dist = dist->d[idx[i]][idx_j]; /* upper diagonal (d[row][col] --> row < col) */
        else                new_dist = dist->d[idx_j][idx[i]];
        if (idx[i] < idx_i) { row = idx[i]; col = idx_i; idx_col = min_row; }
        else                { col = idx[i]; row = idx_i; idx_col = i; }
        if ((row < col) && (new_dist < dist->d[row][col])) dist->d[row][col] = new_dist; /* skip row==col; if new > dist, then dist=dist (=MIN()) */
      }
      else { /* UPGMA: distance to new node is average between elements */
        if (idx[i] < idx_j) new_dist = gsize[idx_j] * dist->d[idx[i]][idx_j];
        else                new_dist = gsize[idx_j] * dist->d[idx_j][idx[i]];
        if (idx[i] < idx_i) { row = idx[i]; col = idx_i; idx_col = min_row; }
        else                { col = idx[i]; row = idx_i; idx_col = i; }
        if (row < col) dist->d[row][col] = (new_dist + (gsize[idx_i] * dist->d[row][col]))/gs1; /* UPGMA distance -- skip row==col */
      }

      /* check if any new minimum is found */  
      if ((dist->d[row][col] < dst_by_row[row])) { dst_by_row[row] = dist->d[row][col]; min_by_row[row] = idx_col; }

      /* rows whose minimum value was min_row need to be updated */
      if ((idx[i] < (dist->size - 1)) && ((min_by_row[idx[i]] == min_row) || (min_by_row[idx[i]] == min_col) || (min_by_row[idx[i]] >= n_idx))) {
        dst_by_row[idx[i]] = 1.e35;
        for (row = 0; row < n_idx; row++) /* "row" is just a recycled var (like i,j,k,l); no special meaning */
          if (((col=idx[row]) > idx[i]) && (dist->d[idx[i]][col] < dst_by_row[idx[i]])) { 
            dst_by_row[idx[i]] = dist->d[idx[i]][col]; min_by_row[idx[i]] = row; /* min_by_row has only row < col */
          }
      }
    } 

    gsize[idx_i] += gsize[idx_j];/* update group size of new internal node */

  } // while (n_idx > 2)
  /* now idx[] has only two elements and "parent" is now (2*ntax - 2) */
  create_parent_node_from_children (tree, parent, idxtree[0], idxtree[1]);
  tree->root = tree->nodelist[parent];

  if (idx[0] < idx[1]) dst_row = dist->d[idx[0]][idx[1]];
  else                 dst_row = dist->d[idx[1]][idx[0]];
  tree->blength[idxtree[0]] = dst_row/2. - height[idx[0]]; /* UPGMA distance */ 
  tree->blength[idxtree[1]] = dst_row/2. - height[idx[1]]; /* UPGMA distance */

  update_topology_sisters   (tree);
  update_topology_traversal (tree);

  if (gsize)   free (gsize);
  if (height)  free (height);
  if (dst_by_row) free (dst_by_row);
}

void
bionj_from_distance_matrix (topology tree, distance_matrix dist) 
{ /* always use upper diagonal of distance_matrix(that is, only i < j in d[i][j]) */
  int i, j, parent = tree->nleaves, n_idx = tree->nleaves, i1, i2, b1, b2, // b1, b2 are best, b1 < b2
      *idx = tree->index,                            /* indexes in UPGMA */
      *idxtree = tree->index + tree->nleaves;        /* indexes in tree (since have values > nleaves) */
  double **delta, Q_min, Q_ij, var_1_2, diff_1_2, blen_1, blen_2, lambda;

  /* tree->index is also used by quasi_randomize_topology(), and here we tell it the info was destroyed */
  tree->quasirandom = false;
  
  if (!tree->blength) topology_malloc_blength (tree);

  delta = (double **) biomcmc_malloc (n_idx * sizeof (double*)); /* delta matrix with dists, variances and sums of distances */
  for (i=0; i < n_idx; i++) delta[i] = (double *) biomcmc_malloc (n_idx * sizeof (double));
  /* delta has distances in upper diag, var in lower diag and sum in diagonal (opposite of original BIONJ C program!) */
  for (i=0; i < n_idx; i++) for (j=i+1; j<n_idx; j++) delta[i][j] = delta[j][i] = dist->d[i][j]; // only upper diagonal of dist is used 
  for (i=0; i < n_idx; i++) delta[i][i] = 0.; 

  for (i=0; i < n_idx; i++) { 
    idx[i]     = i; /* index to actual vector element for UPGMA distance matrix */
    idxtree[i] = i; /* index to actual vector element for tree nodes */
  }

  while (n_idx > 2) { /* choose two nodes to be connected */
    /* update sums: diagonal values of delta will be the sum of distances */
    for (i=0; i < n_idx; i++) {
      delta[idx[i]][idx[i]] = 0;
      for (j=0; j < n_idx; j++) if (j!=i) { // idx(i) < idx(j) for dissimilarities
        if (idx[i] < idx[j]) delta[idx[i]][idx[i]] += delta[idx[i]][idx[j]]; 
        else                 delta[idx[i]][idx[i]] += delta[idx[j]][idx[i]];
      }
    }
    /* find pair that minimises agglomerative criterion -- matrix Q_ij */ 
    Q_min = 1.e64;
    for (i=0; i < n_idx; i++) for (j=0; j < i; j++) {
      if (idx[i] < idx[j]) { i1 = i; i2 = j; } // idx[i1] < idx[i2] always
      else                 { i1 = j; i2 = i; }
      Q_ij = (double)(n_idx - 2) * delta[idx[i1]][idx[i2]] - delta[idx[i1]][idx[i1]] - delta[idx[i2]][idx[i2]];
      if (Q_ij < Q_min - 1.e-8) { Q_min = Q_ij; b1 = i1; b2 = i2; }
    }
    //for (i=0; i < n_idx; i++) {
    //  for (j=0; j < n_idx; j++) { printf ("%9.8lf ", delta[idx[i]][idx[j]]); } printf (" <-- \n");
    //}
    //printf ("chosen: %d %d\n", b1, b2);
    diff_1_2 = (delta[idx[b1]][idx[b1]] - delta[idx[b2]][idx[b2]])/(double)(n_idx-2);
    blen_1 = 0.5 * (delta[idx[b1]][idx[b2]] + diff_1_2);
    blen_2 = 0.5 * (delta[idx[b1]][idx[b2]] - diff_1_2);
    /* calculate lambda */
    var_1_2 = delta[idx[b2]][idx[b1]];  // variance between b1 and b2
    if(var_1_2 < 1.e-12) lambda=0.5; // delta[b2][b1] is var between b1 and b2
    else {
      lambda = 0.;
      for (i=0; i< n_idx; i++) if(b1 != i && b2 != i) {
        if (idx[i] < idx[b1]) lambda += delta[idx[b1]][idx[i]]; // lambda += (var(b1,i) - var(b2,i)
        else                  lambda += delta[idx[i]][idx[b1]];
        if (idx[i] < idx[b2]) lambda -= delta[idx[b2]][idx[i]];
        else                  lambda -= delta[idx[i]][idx[b2]];
      }
      lambda = 0.5 + lambda/(2.*(double)(n_idx-2)* delta[idx[b2]][idx[b1]]);
    }
    if(lambda > 1.0) lambda = 1.0;
    if(lambda < 0.0) lambda = 0.0;
    /* update distances and variances of b1, which will be new node (b2 will be replaced) */
    for (i=0; i< n_idx; i++) if(b1 != i && b2 != i) {
      if (idx[b1] < idx[i]) {i1 = b1; i2 = i;}
      else                  {i2 = b1; i1 = i;} // idx[i1] < idx[i2] always
      /* Distance update --> i<j in delta[i][j] */
      delta[idx[i1]][idx[i2]] = lambda * (delta[idx[i1]][idx[i2]] - blen_1);
      if (idx[b2] < idx[i]) delta[idx[i1]][idx[i2]] += (1. - lambda) * (delta[idx[b2]][idx[i]] - blen_2); // distance(b2,i)
      else                  delta[idx[i1]][idx[i2]] += (1. - lambda) * (delta[idx[i]][idx[b2]] - blen_2);
      /* Variance update --> i > j in delta[i][j] */
      delta[idx[i2]][idx[i1]] =  lambda * (delta[idx[i2]][idx[i1]] - (1.-lambda) * var_1_2);
      if (idx[b2] < idx[i]) delta[idx[i2]][idx[i1]] += (1. - lambda) * (delta[idx[i]][idx[b2]]); // variance(b2,i)
      else                  delta[idx[i2]][idx[i1]] += (1. - lambda) * (delta[idx[b2]][idx[i]]); // variance(b2,i)
    }
    // do we need to make sure that b1 < b2 (not only idx[b1]<idx[b2]) otherwise we may be updating last index n_idx-1 which will  be  neglected ??
    /* tree node creation */
    create_parent_node_from_children (tree, parent, idxtree[b1], idxtree[b2]);
    tree->blength[idxtree[b1]] = blen_1;
    tree->blength[idxtree[b2]] = blen_2;
    idxtree[b1] = parent; /* new node we just created (available to sampling on next iteration) */
    idxtree[b2] = idxtree[--n_idx]; /* avoid replacement (should come last since we may have i == n_idx-1 */
    parent++; /* parent on next iteration */
    idx[b2]    = idx[n_idx]; /* avoid replacement */

  } // while (n_idx > 2)

  /* now idx[] has only two elements and "parent" is now (2*ntax - 2) */
  create_parent_node_from_children (tree, parent, idxtree[0], idxtree[1]);
  tree->root = tree->nodelist[parent];

  if (idx[0] < idx[1]) tree->blength[idxtree[0]] = tree->blength[idxtree[1]] = delta[idx[0]][idx[1]];
  else                 tree->blength[idxtree[0]] = tree->blength[idxtree[1]] = delta[idx[1]][idx[0]];

  update_topology_sisters   (tree);
  update_topology_traversal (tree);
  if (delta) {
    for (i = n_idx - 1; i >= 0; i--) if (delta[i]) free (delta[i]);
    free (delta);
  }
}

void
fill_species_dists_from_gene_dists (distance_matrix spdist, distance_matrix gendist, int *sp_id, bool use_upper_gene)
{
  int i, j, k, row, col, *freq;

  freq = (int*) biomcmc_malloc (spdist->size * sizeof (int));
  for (i = 0; i < spdist->size; i++) freq[i] = 0; /* species frequency for this gene */
  for (i = 0; i < gendist->size; i++) freq[ sp_id[i] ]++; /* used to calculate mean */
  for (i = 0; i < spdist->size; i++) {
    for (j = 0; j <= i; j++)     spdist->d[i][j] = 0.; /* lower diag are mean values */
    for (;j < spdist->size; j++) spdist->d[i][j] = 1.e35; /* upper diag are minimum values */
  }
  
  for (j=1; j < gendist->size; j++) for (i=0; i < j; i++) if (sp_id[i] != sp_id[j]) {
    if (sp_id[i] < sp_id[j]) { row = sp_id[i]; col = sp_id[j]; } /* [row][col] of sptree is upper triangular for minimum */
    else                     { row = sp_id[j]; col = sp_id[i]; }
    if (!use_upper_gene) { k = i; i = j; j = k; } /* then i should be larger than j -- swap values */
    if (gendist->d[i][j] < spdist->d[row][col]) spdist->d[row][col] = gendist->d[i][j]; /* upper diag = minimum */
    spdist->d[col][row] += gendist->d[i][j]; /* lower diag = mean */
  }

  for (i = 0; i < spdist->size; i++) for (j = 0; j < i; j++) if (freq[i] && freq[j]) spdist->d[i][j] /= (double)(freq[i] * freq[j]); 

#ifdef BIOMCMC_PRINT_DEBUG
  for (i=0; i < gendist->size; i++) printf ("spdistfromgene %d\t -> %d\n", i, sp_id[i]);
  for (j=1; j < spdist->size; j++)  for (i=0; i < j; i++) 
    printf ("spdistfromgene (%d\t%d)\t%lf\n", i, j, spdist->d[i][j]); 
#endif
  free (freq);
}

void
update_species_dists_from_spdist (distance_matrix global, distance_matrix local, int *spexist)
{ 
  int i, j;
  if (global->size != local->size) biomcmc_error ("species distance matrices have different sizes within and across loci");

  for (i = 0; i < local->size; i++) for (j = 0; j < i; j++) if (spexist[i] && spexist[j]) { 
    if (global->d[j][i] > local->d[j][i]) global->d[j][i] = local->d[j][i]; /* upper triangular => minimum */
    global->d[i][j] += local->d[i][j]; /* just the sum; to have the mean we need to divide by representativity of each species across loci */
    // // guenomu receives another matrix // if (counter) { counter->d[i][j] += 1.; counter->d[j][i] += 1.; }
  }
}

int
prepare_spdistmatrix_from_gene_species_map (spdist_matrix spdist, int *sp_id, int n_sp_id)
{
  int i, number_of_species_present_in_gene = 0;
  for (i = 0; i < spdist->size; i++) spdist->species_present[i] = false; // update presence mask of species in gene 
  for (i = 0; i < n_sp_id; i++) spdist->species_present[ sp_id[i] ] = true;
  for (i = 0; i < spdist->size; i++) if (spdist->species_present[i]) number_of_species_present_in_gene++;
  return number_of_species_present_in_gene;
}

void
fill_spdistmatrix_from_gene_dists (spdist_matrix spdist, distance_matrix gendist, int *sp_id, bool use_upper_gene)
{ // more compact than functions above (which could be eliminated in future versions)
  int i, j, i2, j2, idx, row, col, n_pairs = spdist->size*(spdist->size-1)/2;

  for (i = 0; i < n_pairs; i++) {
    spdist->mean[i] = 0;
    spdist->min[i] = 1.e35;
    spdist->count[i] = 0;
  }
  
  for (j=1; j < gendist->size; j++) for (i=0; i < j; i++) if (sp_id[i] != sp_id[j]) {
    if (sp_id[i] < sp_id[j]) { row = sp_id[i]; col = sp_id[j]; } /* make sure that row < col */
    else                     { row = sp_id[j]; col = sp_id[i]; }
    i2 = i; j2 = j;  // i2 < j2 if upper and i2 > j2 if lower diag is used
    if (!use_upper_gene) { i2 = j; j2 = i; } /* then i2 should be larger than j2 -- swap values */
    idx = col * (col-1)/2 + row; /* index in spdist */
    if (gendist->d[i2][j2] < spdist->min[idx]) spdist->min[idx] = gendist->d[i2][j2];
    spdist->mean[idx] += gendist->d[i2][j2];
    spdist->count[idx]++;
  }

  for (i = 0; i < n_pairs; i++) if (spdist->count[i]) spdist->mean[i] /= spdist->count[i];
  return;
}

void
update_spdistmatrix_from_spdistmatrix (spdist_matrix global, spdist_matrix local)
{ 
  int i, j, idx;
  if (global->size != local->size) biomcmc_error ("species spdist matrices have different sizes within and across loci");

  for (j = 1; j < local->size; j++) for (i = 0; i < j; i++) if (local->species_present[i] && local->species_present[j]) { 
    idx = j * (j-1) /2 + i; // index in 1D vector without diagonals (for diagonals replace -1 for +1 BTW) 
    global->mean[idx] += local->mean[idx]; // global only stores average across locals (min => within locus)
    global->min[idx] += local->min[idx];
    global->count[idx]++;
  }
  for (i = 0; i < global->size; i++) global->species_present[i] |= local->species_present[i]; // overall presence of species 
}

void
create_parent_node_from_children (topology tree, int parent, int lchild, int rchild)
{
  tree->nodelist[parent]->left  = tree->nodelist[lchild];
  tree->nodelist[parent]->right = tree->nodelist[rchild];
  tree->nodelist[lchild]->up = tree->nodelist[parent];
  tree->nodelist[rchild]->up = tree->nodelist[parent];
  tree->nodelist[rchild]->sister = tree->nodelist[lchild];
  tree->nodelist[lchild]->sister = tree->nodelist[rchild];
}

void
topology_apply_rerooting (topology tree, bool update_done)
{
  int i, n1 = 0, n_valid = 0, n_invalid = 0,
      *valid = tree->index, *invalid = tree->index + tree->nnodes; /* index has size 4*nleaves = 2*nnodes + 2 */

  /* tree->index is also used by quasi_randomize_topology(), and here we tell it the info was destroyed */
  tree->quasirandom = false;

  invalid[n_invalid++] = tree->root->id; /* root and its children are the only forbidden regraft nodes */ 
  invalid[n_invalid++] = tree->root->left->id;
  invalid[n_invalid++] = tree->root->right->id;
  /* create valid[] vector of indexes by exclusion of invalid[] */ 
  qsort (invalid, n_invalid, sizeof (int), compare_int_increasing);
  for (i = 0; i < tree->nnodes; i++) { /* every node "i" is OR valid OR invalid */
    if ((n1 < n_invalid) && (i == invalid[n1])) n1++; /* node "i" is invalid: skip it (that's why we qsorted) */
    else valid[n_valid++] = i;                        /* node "i" is valid */
  }

  n1 = valid[ biomcmc_rng_unif_int (n_valid) ]; /* draw regraft node */

  apply_spr_at_nodes_LCAprune (tree, tree->root, tree->nodelist[n1], update_done);
  update_topology_traversal (tree);
}

void
topology_apply_shortspr (topology tree, bool update_done)
{
  bool success = false;
  double scale = 1./tree->nleaves;
  int i;
  /* tree->index is also used by quasi_randomize_topology(), and here we tell it the sequence was destroyed */
  tree->quasirandom = false; /* IOW, quasi_randomize() will have to start from scratch */
  if (!tree->traversal_updated) update_topology_traversal (tree);

  for (i = 0; (!success) && (i < 4); i++) {
    success = topology_apply_shortspr_weighted_subtree (tree, tree->root, NULL, scale, update_done);
    scale *= 2.;
  }
  if (!success) topology_apply_shortspr_weighted_subtree (tree, tree->root, NULL, 1., update_done);
  update_topology_traversal (tree);
}

void
topology_apply_shortspr_weighted (topology tree, double *prob, bool update_done)
{
  bool success = false;
  double scale = 1., *localprob = NULL;
  int i;

  /* tree->index is also used by quasi_randomize_topology(), and here we tell it the sequence was destroyed */
  tree->quasirandom = false; /* IOW, quasi_randomize() will have to start from scratch */
  if (!tree->traversal_updated) update_topology_traversal (tree);

  for (i = 0; (!success) && (i < 8); i++) {
    success = topology_apply_shortspr_weighted_subtree (tree, tree->root, prob, scale, update_done);
    scale *= 2.;
  }
  if (success) { update_topology_traversal (tree); return; }
  /* something is wrong, prob values too small maybe. */
  localprob = (double*) biomcmc_malloc (tree->nnodes * sizeof(double));
  for (i = 0; i < tree->nnodes; i++) localprob[i] = prob[i] + 1.;
  scale = 1./(double) (tree->nleaves);
  for (i = 0; (!success) && (i < 8); i++) {
    success = topology_apply_shortspr_weighted_subtree (tree, tree->root, localprob, scale, update_done);
    scale *= 2.;
  }
  update_topology_traversal (tree);
  if (localprob) free (localprob);
}

bool
topology_apply_shortspr_weighted_subtree (topology tree, topol_node lca, double *prob, double scale, bool update_done)
{  /* theoretically it works even if lca is root, but may leave the (unrooted version of the) tree unchanged */ 
  bool success = false;
  double prob_occurence, p_l, p_r;

  if (!lca->left->internal) return false; /* if left is a leaf then right is also a leaf (since no swaps were applied to lca yet */
  else                      success |= topology_apply_shortspr_weighted_subtree (tree, lca->left, prob, scale, update_done);

  if (lca->right->internal) success |= topology_apply_shortspr_weighted_subtree (tree, lca->right, prob, scale, update_done);
  else {  /* three leaves => only 3 possible topols => only 2 possible swaps: ((A,B),C) -> swap to (A,(B,C)) OR (B,(A,C)) */
    /* 1) the 2 possible swaps can be done by chosing one of the two leaves at subtree as regraft: A or B in ((A,B),C)
     * 2) although we called this function on children, we can still trust left subtree is larger than right */ 
    if (prob) p_l = scale * prob[lca->left->id]; 
    else      p_l = scale; /* if prob == NULL then we just use scale as common prob of swap */
    if (biomcmc_rng_unif_pos32 () < p_l) { /* branch swap will occur */
      if (biomcmc_rng_unif_int (2)) apply_spr_at_nodes_LCAprune (tree, lca, lca->left->left, update_done);
      else                          apply_spr_at_nodes_LCAprune (tree, lca, lca->left->right, update_done);
      return true;
    }
    return success; /* failed on this node, but may have succeeded before */
  }
  /* otherwise both children are internal */
  prob_occurence = biomcmc_rng_unif_pos32 ();
  if (prob) {
    p_l = scale * prob[lca->left->id];
    p_r = scale * prob[lca->right->id];
  }
  else p_l = p_r = scale;
  /* only left = L*(1-R) ; only right = R*(1-L) ; both = L*R ; none = (1-L)*(1-R) 
   * like in the prob line (stick break reprsentation):
   * |--------------|------*-----|-------------|--------| 
   * | only left    |    both    |  only right |  none  |        (1)
   * |---------------------*-------------------|         
   * |   left or both      |   right or both   |                 (2)
   * (note that "both" is split half between right and left) */
  if (prob_occurence < p_l + p_r - p_l * p_r) { /* at least one event happens [L(1-R)+R(1-L)+L*R] */
    topol_node np, nr;
    topol_node prunenode[2], regraftnode[3];
    if (prob_occurence < p_l - (p_l * p_r / 2.)) { /* only left subtree OR both with half prob [L(1-R)+L*R/2] -- marked (2) above */
      regraftnode[0] = lca->left; /* after SPR will be at right side */
      regraftnode[1] = lca->left->left;
      regraftnode[2] = lca->left->right;
      prunenode[0] = lca->right->left;
      prunenode[1] = lca->right->right;
      if (biomcmc_rng_unif_int (2)) apply_spr_at_nodes_LCAprune (tree, lca, lca->left->left, update_done);
      else                          apply_spr_at_nodes_LCAprune (tree, lca, lca->left->right, update_done);
    }
    else { /* only right subtree or both with other half prob [R(1-L) + L*R/2] -- marked (2) above  */
      regraftnode[0] = lca->right; /* after SPR will be at left side */
      regraftnode[1] = lca->right->left;
      regraftnode[2] = lca->right->right;
      prunenode[0] = lca->left->left;
      prunenode[1] = lca->left->right;
      if (biomcmc_rng_unif_int (2)) apply_spr_at_nodes_LCAprune (tree, lca, lca->right->left, update_done);
      else                          apply_spr_at_nodes_LCAprune (tree, lca, lca->right->right, update_done);
    }
    if ((prob_occurence > (p_l - (p_l * p_r))) && (prob_occurence < p_l)) { /* both events happen -- marked (1) above */
      np =   prunenode[ biomcmc_rng_unif_int (2)]; 
      nr = regraftnode[ biomcmc_rng_unif_int (3)]; 
      apply_spr_at_nodes_notLCAprune (tree, np, nr, update_done); 
    }
    return true;
  }
  return success; /* failed on this node, but may have succeeded before */
}

void
topology_apply_spr_on_subtree (topology tree, topol_node lca, bool update_done)
{ 
  int i, n1, n2, *valid = tree->index, *invalid, *regraft, n_valid = 0, n_invalid = 0, n_regraft = 0;
  topol_node first_child = lca;

  /* tree->index is also used by quasi_randomize_topology(), and here we tell it the info was destroyed */
  tree->quasirandom = false;

  if (lca == tree->root) biomcmc_error ("root node is not eligible for SPR move (maybe root->left or root->right?)");
  if (!tree->traversal_updated) update_topology_traversal (tree);
  if (lca->split->n_ones < 3) return; /* it should not complain if subtree is a "leaf" or has only two children */
 
  if (lca->split->n_ones == 3) { /* three leaves => only 3 possible topols => only 2 possible swaps */
    /* 1) the 2 possible swaps can be done by chosing one of the two leaves at subtree as regraft: A or B in ((A,B),C)
     * 2) lca->left is internal (by design of traversal, see bipartition_is_larger() function); */
    if (biomcmc_rng_unif_int (2)) apply_spr_at_nodes_LCAprune (tree, lca, lca->left->left, update_done);
    else                          apply_spr_at_nodes_LCAprune (tree, lca, lca->left->right, update_done);
    return;
  }

  /* find number of nodes below subtree (using postorder info) */
  while (first_child->internal) first_child = first_child->left; /* walk down the subtree (left is first, in postorder) */ 

  /* all nodes below lca node (inclusive) are eligible prune nodes */
  for (i = first_child->up->mid[0]; i <= lca->mid[0]; i++) { 
    valid[n_valid++] = tree->postorder[i]->id; /* postorder[] stores only internal nodes */
    if (!tree->postorder[i]->left->internal)  valid[n_valid++] = tree->postorder[i]->left->id;
    if (!tree->postorder[i]->right->internal) valid[n_valid++] = tree->postorder[i]->right->id;
  }
#ifdef BIOMCMC_DEBUG
  if (n_valid != (2*lca->split->n_ones - 1)) biomcmc_error ("%d nodes eligible in subtree with %d leaves (SPR)",n_valid, lca->split->n_ones);
#endif
  n1 = valid[ biomcmc_rng_unif_int (n_valid) ]; /* draw prune node */

  /* n_valid is always smaller than nnodes so we're safe (tree->index has only 4*nleaves = 2*nnodes + 2 elements */
  invalid = tree->index + n_valid; /* invalid[] vector points to after last element of valid[] */

  /* remove prune node's possible four immediate neighbors (left, right, up and sister), lca included */
  invalid[n_invalid++] = n1; /* trivial restriction (regraft != prune) */
  if (tree->nodelist[n1]->internal) { 
    invalid[n_invalid++] = tree->nodelist[n1]->left->id;
    invalid[n_invalid++] = tree->nodelist[n1]->right->id;
  }
  if (tree->nodelist[n1] != lca) {
    invalid[n_invalid++] = tree->nodelist[n1]->up->id;
    invalid[n_invalid++] = tree->nodelist[n1]->sister->id;
  }

  regraft = invalid + n_invalid; /* regraft[] starts after last element of invalid[] */

  /* create regraft[] vector of indexes by exclusion of invalid[] from valid[] */ 
  qsort (invalid, n_invalid, sizeof (int), compare_int_increasing);
  qsort (  valid,   n_valid, sizeof (int), compare_int_increasing);
  n2 = 0;
  for (i = 0; i < n_valid; i++) { /* every node is OR valid OR invalid */
    if ((n2 < n_invalid) && (valid[i] == invalid[n2])) n2++; /* skip invalid */
    else regraft[n_regraft++] = valid[i];
  }

  n2 = regraft[ biomcmc_rng_unif_int (n_regraft) ]; /* regraft node (here valid[] already has idx in postorder) */

  /* this wrapper function will decide if prune node is LCA or not of regraft */
  apply_spr_at_nodes (tree, tree->nodelist[n1], tree->nodelist[n2], update_done);
  update_topology_traversal (tree);
  return;
}

void
topology_apply_spr (topology tree, bool update_done)
{ 
  uint32_t n_left, n_right;

  if (cant_apply_swap (tree)) return; /* checks if any subtree has more than two leaves */
  n_left  = tree->root->left->split->n_ones;
  n_right = tree->root->right->split->n_ones;

  if     (n_right < 3) topology_apply_spr_on_subtree (tree, tree->root->left, update_done);
  else if (n_left < 3) topology_apply_spr_on_subtree (tree, tree->root->right, update_done);
  else if (biomcmc_rng_unif_int (n_right + n_left) < n_left) topology_apply_spr_on_subtree (tree, tree->root->left, update_done);
  else                                                       topology_apply_spr_on_subtree (tree, tree->root->right, update_done);
}

void
topology_apply_spr_unrooted (topology tree, bool update_done)
{ /* neglects root node and applies SPR that changes eq. unrooted topology */
  int i, n1, n2, *valid = tree->index, *invalid, *regraft, n_valid = 0, n_invalid = 0, n_regraft = 0;

  if (tree->nleaves < 4) return; /* There is only one unrooted triplet */
  if (!tree->traversal_updated) update_topology_traversal (tree);
  /* tree->index is also used by quasi_randomize_topology(), and here we tell it the info was destroyed */
  tree->quasirandom = false;
  if (tree->nleaves == 4) { /* special case, furthermore can be simplified */
    if (!tree->root->right->internal) { /* (((a,b),c),d) --> a or b regrafted to d (OR,equiv, regrafted to c) */
      valid[n_valid++] = tree->root->left->left->left->id; 
      valid[n_valid++] = tree->root->left->left->right->id;
      n1 = valid[ biomcmc_rng_unif_int (n_valid) ]; /* draw prune node */
      apply_spr_at_nodes_notLCAprune (tree, tree->nodelist[n1], tree->root->right, update_done);
      update_topology_traversal (tree);
      return;
    }
    else { /* ((a,b),(c,d)) --> a or b regrafted to c (again, since it's unrooted, it's equiv to regrafting to d) */
      valid[n_valid++] = tree->root->left->left->id; 
      valid[n_valid++] = tree->root->left->right->id;
      n1 = valid[ biomcmc_rng_unif_int (n_valid) ]; /* draw prune node */
      apply_spr_at_nodes_notLCAprune (tree, tree->nodelist[n1], tree->root->right->left, update_done);
      update_topology_traversal (tree);
      return;
    }
  }
  /* all nodes except root are eligible prune nodes */
  for (i = 0; i < tree->nnodes; i++) if (i != tree->root->id) valid[n_valid++] = i;
  n1 = valid[ biomcmc_rng_unif_int (n_valid) ]; /* draw prune node */

  /* tree->index has (4*nleaves) = (2*nnodes + 2) elements so we're safe since n_valid <= nnodes */
  invalid = tree->index + n_valid; /* invalid[] vector points to after last element of valid[] */

  /* remove prune node's possible four immediate neighbors (left, right, up and sister) */
  invalid[n_invalid++] = n1; /* trivial restriction (regraft != prune) */
  invalid[n_invalid++] = tree->nodelist[n1]->up->id; /* since prune cannot be root node */
  invalid[n_invalid++] = tree->nodelist[n1]->sister->id;
  if (tree->nodelist[n1]->internal) { 
    invalid[n_invalid++] = tree->nodelist[n1]->left->id;
    invalid[n_invalid++] = tree->nodelist[n1]->right->id;
  }
  if ((tree->nodelist[n1]->up == tree->root) && (tree->nodelist[n1]->sister->internal)) { /* equiv. to rerooting */ 
    invalid[n_invalid++] = tree->nodelist[n1]->sister->left->id;
    invalid[n_invalid++] = tree->nodelist[n1]->sister->right->id;
  }
  else if (tree->nodelist[n1]->up->up == tree->root) { /* equiv. to rerooting */ 
    invalid[n_invalid++] = tree->nodelist[n1]->up->sister->id;
  }

  regraft = invalid + n_invalid; /* regraft[] starts after last element of invalid[] */

  /* create regraft[] vector of indexes by exclusion of invalid[] from all nodes (including root now) */ 
  qsort (invalid, n_invalid, sizeof (int), compare_int_increasing);
  n2 = 0;
  for (i = 0; i < tree->nnodes; i++) { /* every node is OR valid OR invalid */
    if ((n2 < n_invalid) && (invalid[n2] == i)) n2++; /* skip invalid */
    else regraft[n_regraft++] = i;
  }

  n2 = regraft[ biomcmc_rng_unif_int (n_regraft) ]; /* regraft node (here valid[] already has idx in postorder) */

  /* this wrapper function will decide if prune node is LCA or not of regraft */
  apply_spr_at_nodes (tree, tree->nodelist[n1], tree->nodelist[n2], update_done);
  update_topology_traversal (tree);
  return;
}

void
topology_apply_nni (topology tree, bool update_done)
{
  int i, n1, n2, *valid = tree->index, *invalid = tree->index + tree->nnodes, n_valid = 0, n_invalid = 0, lca_idx = 0;
  /* please read comment on function topology_apply_spr() about the two vectors and the need to order them. I always
   * find this solution ugly, but at least I think is correct. The bugs didn't appear here with the tempting
   * no-replacement algorithm, but in apply_spr() so I suspect it was only with the regraft sampling (leomrtns). 
   *
   * UPDATE 2010.07.27: this function seems to be buggy (does not respect rooting) */

  /* tree->index is also used by quasi_randomize_topology(), and here we tell it the info was destroyed */
  tree->quasirandom = false;

  invalid[n_invalid++] = tree->root->id; /* root is forbidden prune node (otherwise we have rerooting) */
  /* if a child of root is pruned, the regraft must be on subtree below (at least two nodes apart from prune) */
  if ((!tree->root->left->internal) || (!tree->root->left->left->internal && !tree->root->left->right->internal)) 
    invalid[n_invalid++] = tree->root->left->id;
  if ((!tree->root->right->internal) || (!tree->root->right->left->internal && !tree->root->right->right->internal)) 
    invalid[n_invalid++] = tree->root->right->id;

  /* create valid[] vector of indexes by exclusion of invalid[] */ 
  if (n_invalid > 1) qsort (invalid, n_invalid, sizeof (int), compare_int_increasing);
  n1 = 0;
  for (i = 0; i < tree->nnodes; i++) { /* every node is OR valid OR invalid */
    if ((n1 < n_invalid) && (i == invalid[n1])) n1++; /* node "i" is invalid, skip (that's why we qsorted) */
    else valid[n_valid++] = i;                        /* node "i" is valid */
  }

  n1 = valid[ biomcmc_rng_unif_int (n_valid) ]; /* draw prune node */

  n_valid = 0;
  /* NNI: at most four possible regraft nodes of LCA type */
  if (tree->nodelist[n1]->internal) {
    if (tree->nodelist[n1]->left->internal) {
      valid[n_valid++] = tree->nodelist[n1]->left->left->id;
      valid[n_valid++] = tree->nodelist[n1]->left->right->id;
    }
    if (tree->nodelist[n1]->right->internal) {
      valid[n_valid++] = tree->nodelist[n1]->right->left->id;
      valid[n_valid++] = tree->nodelist[n1]->right->right->id;
    }
    lca_idx = n_valid; /* up to this index, spr is of LCA type */
  }

  if (tree->nodelist[n1]->up != tree->root) { /* root's children are only allowed to LCA-type moves */ 
    valid[n_valid++] = tree->nodelist[n1]->up->sister->id;
    if (tree->nodelist[n1]->up->up != tree->root) 
    valid[n_valid++] = tree->nodelist[n1]->up->up->id;
    if (tree->nodelist[n1]->sister->internal) {
    valid[n_valid++] = tree->nodelist[n1]->sister->left->id;
    valid[n_valid++] = tree->nodelist[n1]->sister->right->id;
    }
  }

  i = biomcmc_rng_unif_int (n_valid); /* draw regraft node */
  n2 = valid[i];

  if (i < lca_idx) apply_spr_at_nodes_LCAprune (tree, tree->nodelist[n1], tree->nodelist[n2], update_done);
  else          apply_spr_at_nodes_notLCAprune (tree, tree->nodelist[n1], tree->nodelist[n2], update_done);
  update_topology_traversal (tree);
}

bool
cant_apply_swap (topology tree)
{
  if (!tree->traversal_updated) update_topology_traversal (tree);
  if ((tree->root->left->split->n_ones < 3) && (tree->root->left->split->n_ones < 3)) return true;
  return false;
}

