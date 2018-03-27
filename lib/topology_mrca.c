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

#include "topology_mrca.h"

void gene_tree_reconcile_unrooted (topology gene, topology species);
void prepare_for_loss_calculation (topology gene, topology species);

/* * * MRCA functions * * */

void 
new_mrca_for_topology (topology t)
{
  int i, j;

  if (t->mrca) biomcmc_error ("trying to malloc() mrca structure already allocated for topology");
  /* distance between nodes i and j is given by lca[i-1][j] (i should be larger than j) */
  t->mrca = (topol_node**) biomcmc_malloc ( (t->nnodes - 1) * sizeof (topol_node*));

  for (i=0; i < t->nnodes - 1; i++) t->mrca[i] = (topol_node*) biomcmc_malloc ( (i+1) * sizeof (topol_node));
  for (i=0; i < t->nnodes - 1; i++) for (j=0; j <= i; j++) t->mrca[i][j] = NULL;
}

topol_node
mrca_between_nodes (topology topol, int i, int j)
{
  topol_node p;
  bipartition lowsplit;
  int i1, j1;

  if (i == j) return topol->nodelist[i];
  if (j > i) { int tmp = i; i = j; j = tmp; } /* i should always be larger than j */

  if ((!topol->traversal_updated)) update_topology_traversal (topol); 

  if (topol->mrca[i-1][j]) return topol->mrca[i-1][j];

  /* Minimize sequencial search by choosing node closer to root (topol_node::level = distance from root) */
  if (topol->nodelist[i]->level > topol->nodelist[j]->level) 
   { i1 = i; j1 = j; } /* node[i1] will be further from root */
  else 
   { i1 = j; j1 = i; } /* node[i1] will be further from root */
 
  p = topol->nodelist[j1]; /* here we start to actually search for the lca */

  /* if node[i1] - which is further from root - is internal we must compare fully both bipartitions - O(n) */
  if ((topol->nodelist[i1]->internal) && (lowsplit = topol->nodelist[i1]->split))
    while ((p) && (!bipartition_contains_bits (p->split, lowsplit))) p = p->up; /* climb up the tree */
  else /* if node[i1] is a leaf the comparison is faster - O(1) - for very large trees */
    while ((p) && (!bipartition_is_bit_set (p->split, i1))) p = p->up; /* climb up the tree */

  if (!p) biomcmc_error ("Couldn't find the MRCA. Possible bug related to root node."); /* this shouldn't happen(R) */

  return topol->mrca[i-1][j] = p;
}

/* * * gene tree/species tree functions * * */

void
init_tree_recon_from_species_topology (topology gene, topology species)
{  /* safe, since uses index_sptaxa_to_genetaxa() which does not assume sorted spnames */
  if (!gene->rec) gene->rec = new_reconciliation (gene->nleaves, species->nleaves);
  if (!species->mrca) new_mrca_for_topology (species);
  index_sptaxa_to_genetaxa (species->taxlabel, gene->taxlabel, gene->rec->sp_id, NULL);
  initialize_reconciliation_sp_count (gene->rec, species->taxlabel->nstrings, gene->nleaves);
  gene_tree_reconcile (gene, species); 
}

void
init_tree_recon_from_species_names (topology gene, char_vector sptaxlabel)
{  /* safe, since uses index_sptaxa_to_genetaxa() which does not assume sorted spnames */
  /* used by dSPR (topology_splitset) only (for now) */ 
  if (!gene->rec) gene->rec = new_reconciliation (gene->nleaves, sptaxlabel->nstrings);
  index_sptaxa_to_genetaxa (sptaxlabel, gene->taxlabel, gene->rec->sp_id, NULL);
  initialize_reconciliation_sp_count (gene->rec, sptaxlabel->nstrings, gene->nleaves);
}

void
index_sptaxa_to_genetaxa (char_vector species, char_vector gene, int *sp_idx_in_gene, empfreq ef_external)
{ /* We could have used gene->rec (reconciliation) instead of gene->rec->sp_id (int*) as input parameter, but in 
     some cases we don't care about the gene tree (e.g. MaxTree algorithm, where we build the species tree directly 
     from the gene distances) */
  int i, j, n_index = gene->nstrings, *index;
  empfreq ef;

  /* Search first largest species names (so that for example "ecoli" will match only if "ecoliII" doesn't) */
  if (ef_external) ef = ef_external; /* empfreq initialized outside (by calling function, e.g) */
  else ef = new_empfreq_sort_decreasing (species->nchars, species->nstrings, 1); /* local ("1" means size_t) */

  index = (int*) biomcmc_malloc (n_index * sizeof (int));
  for (i=0; i < n_index; i++) { 
    sp_idx_in_gene[i] = -1; /* initialize mapping */ 
    index[i] = i;  /* scan gene leaves _without_ replacement */
  }

  for (i=0; i < species->nstrings; i++) for (j=0; j < n_index; j++) /* search sp name in all unmapped gene names */
    if ((gene->nchars[index[j]] >= species->nchars[ ef->i[i].idx ]) && /* ordered species names (by nchars[]) */
        (strcasestr (gene->string[index[j]], species->string[ ef->i[i].idx ]))) { 
      /* found species name within gene name; we have a mapping */
      sp_idx_in_gene[ index[j] ] = ef->i[i].idx;
      index[j] = index[--n_index]; // with this the whole search takes O(N ln N),
      j--; // index[j] is now a distinct element (the last)
    }

  if (n_index) {
    fprintf (stderr, "Couldn't find species for genes:\n");
    for (i=0; i < n_index; i++) fprintf (stderr, " \"%s\"\n", gene->string[index[i]]);
    biomcmc_error ("gene names should contain the name of species");
  }

  if (!ef_external) del_empfreq (ef); 
  if (index) free (index);
}

void
index_sptaxa_to_reconciliation (char_vector species, char_vector gene, reconciliation rec)
{/* largest species names already sorted (so that for example "ecoli" will match only after "ecoliII" didn't) */
  int i, j, n_index = gene->nstrings, *index;

  index = (int*) biomcmc_malloc (n_index * sizeof (int));
  for (i=0; i < n_index; i++) { 
    rec->sp_id[i] = -1; /* initialize mapping */ 
    index[i] = i;  /* scan gene leaves _without_ replacement */
  }

  for (i=0; i < species->nstrings; i++) for (j=0; j < n_index; j++) /* search sp name in all unmapped gene names */
    if ((gene->nchars[index[j]] >= species->nchars[i]) && /* ordered species names (by nchars[]) */
        (strcasestr (gene->string[index[j]], species->string[i]))) { 
      /* found species name within gene name; we have a mapping */
      rec->sp_id[ index[j] ] = i;
      index[j] = index[--n_index]; // with this the whole search takes O(N ln N),
      j--; // index[j] is now a distinct element (the last)
    }

  if (n_index) {
    fprintf (stderr, "Couldn't find species for genes:\n");
    for (i=0; i < n_index; i++) fprintf (stderr, " \"%s\"\n", gene->string[index[i]]);
    biomcmc_error ("gene names should contain the name of ordered species");
  }
  if (index) free (index);

  initialize_reconciliation_sp_count (rec, species->nstrings, gene->nstrings);
}

void
initialize_reconciliation_sp_count (reconciliation rec, int n_sp, int n_idx)
{
  int i; /* rec->sp_id[i] is the species index for gene i */
  for (i = 0; i < n_sp; i++)  rec->sp_count[i] = 0; /* representativity of each species in gene family */
  for (i = 0; i < n_idx; i++) rec->sp_count[ rec->sp_id[i] ]++; /* update species frequencies */

  rec->sp_size = 0;
  for (i = 0; i < n_sp; i++) if (rec->sp_count[i]) rec->sp_size++;
  rec->size_diff = 2 * (n_idx - rec->sp_size); /* term to calculate deepcoals */ 
}

void
initialize_reconciliation_from_species_tree (topology gene, topology species)
{
  int i;
  if (!species->mrca) new_mrca_for_topology (species);
  for (i=0; i < gene->nleaves; i++) gene->rec->map_d[i] = species->nodelist[ gene->rec->sp_id[i] ];
  gene->rec->current_sptree = species;
}

void
gene_tree_reconcile (topology gene, topology species)
{
  int i, g_id;
  topol_node map_lchild, map_rchild;

  if (!   gene->traversal_updated) update_topology_traversal (gene);
  if (!species->traversal_updated) update_topology_traversal (species);
  if (species != gene->rec->current_sptree) initialize_reconciliation_from_species_tree (gene, species);

  prepare_for_loss_calculation (gene, species);

  for (i=0; i < gene->nleaves-1; i++) {
    g_id = gene->postorder[i]->id; /* node ID on gene tree */
    map_lchild = gene->rec->map_d[ gene->postorder[i]->left->id ]; /* gene->map[] are nodes on species tree */
    map_rchild = gene->rec->map_d[ gene->postorder[i]->right->id ];
    gene->rec->map_d[g_id] = mrca_between_nodes (species, map_lchild->id, map_rchild->id);

    /* cummulative number of duplications below node, following e.g. Bioinformatics.2001.821 */
    gene->rec->ndup_d[g_id] = gene->rec->ndup_d[gene->postorder[i]->left->id] + gene->rec->ndup_d[gene->postorder[i]->right->id]; 
    /* cummul. number of losses below node, following SIAM.2000.729 (not over edges like ACMTransComputBiolBioinfo.2010.14) */
    gene->rec->nlos_d[g_id] = gene->rec->nlos_d[gene->postorder[i]->left->id] + gene->rec->nlos_d[gene->postorder[i]->right->id];

    if ((gene->rec->map_d[g_id] == map_lchild) || (gene->rec->map_d[g_id] == map_rchild)) {
      gene->rec->ndup_d[g_id]++;

      if (map_lchild != map_rchild) { /* if all three are the same there are no losses */
        if (map_lchild == gene->rec->map_d[g_id]) /* number of intermediate nodes in sptree + 1 */
          gene->rec->nlos_d[g_id] += (map_rchild->mid[4] - gene->rec->map_d[g_id]->mid[4]);
        else
          gene->rec->nlos_d[g_id] += (map_lchild->mid[4] - gene->rec->map_d[g_id]->mid[4]);
      }
    } // if (node is a duplication)
    else /* number of nodes, in sptree, between map and map's children ("-2" since difference in levels=1 means NO interm) */
      gene->rec->nlos_d[g_id] += (map_lchild->mid[4] + map_rchild->mid[4] - 2 * gene->rec->map_d[g_id]->mid[4] - 2);
  }

  gene_tree_reconcile_unrooted (gene, species);
}

void
gene_tree_reconcile_unrooted (topology gene, topology species)
{ /* total n_dups rooted at node should be ndup_u + ndup_d + Indicator{ mrca(map_d, map_u) } */
  int i, g_id, r_left = gene->root->left->id, r_right = gene->root->right->id, thisloss, thisdups, thiscoal, 
      min_coal = 0xffffff, min_dups = 0xffffff, min_loss = 0xffffff; /* large number */
  topol_node map_root, map_up, map_sister;

  gene->rec->map_u[r_left]  = gene->rec->map_d[r_right]; /* neglect root node; r_left and r_right have same info */
  gene->rec->map_u[r_right] = gene->rec->map_d[r_left];
  gene->rec->ndup_u[r_left]  = gene->rec->ndup_d[r_right];
  gene->rec->ndup_u[r_right] = gene->rec->ndup_d[r_left];
  gene->rec->nlos_u[r_left]  = gene->rec->nlos_d[r_right];
  gene->rec->nlos_u[r_right] = gene->rec->nlos_d[r_left];


  /* as "rooted" version, but replacing left and right by up and sister (obs: postorder[gene->nleaves-2] => root) */
  for (i = gene->nleaves-3; i >= 0; i--) if ((gene->postorder[i]->id != r_left) && (gene->postorder[i]->id != r_right)) {
    g_id = gene->postorder[i]->id; /* node ID on gene tree */
    map_up     = gene->rec->map_u[ gene->postorder[i]->up->id ]; /* gene->map[] are nodes on species tree */
    map_sister = gene->rec->map_d[ gene->postorder[i]->sister->id ];
    gene->rec->map_u[g_id] = mrca_between_nodes (species, map_up->id, map_sister->id);

    gene->rec->ndup_u[g_id] = gene->rec->ndup_u[gene->postorder[i]->up->id] + gene->rec->ndup_d[gene->postorder[i]->sister->id]; 
    gene->rec->nlos_u[g_id] = gene->rec->nlos_u[gene->postorder[i]->up->id] + gene->rec->nlos_d[gene->postorder[i]->sister->id]; 
    
    if ((gene->rec->map_u[g_id] == map_up) || (gene->rec->map_u[g_id] == map_sister)) {
      gene->rec->ndup_u[g_id]++;

      if (map_up != map_sister) { /* if all three are the same there are no losses */
        if (map_sister == gene->rec->map_u[g_id]) /* then map_up is distinct from mrca(up,sister) */ 
          gene->rec->nlos_u[g_id] += (map_up->mid[4] - gene->rec->map_u[g_id]->mid[4]);
        else
          gene->rec->nlos_u[g_id] += (map_sister->mid[4] - gene->rec->map_u[g_id]->mid[4]);
      }
    } // if (node is a duplication)
    else /* number of nodes, in sptree, between map and map's "children" ("-2" since difference in levels=1 means NO interm) */
      gene->rec->nlos_u[g_id] += (map_sister->mid[4] + map_up->mid[4] - 2 * gene->rec->map_u[g_id]->mid[4] - 2);
  }

  /* we went in preorder over internal nodes; at last, we go over leaves */
  for (i = 0; i < gene->nleaves; i++) if ((gene->nodelist[i]->id != r_left) && (gene->nodelist[i]->id != r_right)) {
    map_up     = gene->rec->map_u[ gene->nodelist[i]->up->id ]; /* gene->map[] are nodes on species tree */
    map_sister = gene->rec->map_d[ gene->nodelist[i]->sister->id ];
    gene->rec->map_u[i] = mrca_between_nodes (species, map_up->id, map_sister->id);

    gene->rec->ndup_u[i] = gene->rec->ndup_u[gene->nodelist[i]->up->id] + gene->rec->ndup_d[gene->nodelist[i]->sister->id]; 
    gene->rec->nlos_u[i] = gene->rec->nlos_u[gene->nodelist[i]->up->id] + gene->rec->nlos_d[gene->nodelist[i]->sister->id];

    if ((gene->rec->map_u[i] == map_up) || (gene->rec->map_u[i] == map_sister)) {
      gene->rec->ndup_u[i]++;

      if (map_up != map_sister) { /* if all three are the same there are no losses */
        if (map_sister == gene->rec->map_u[i]) /* then map_up is distinct from mrca(up,sister) */ 
          gene->rec->nlos_u[i] += (map_up->mid[4] - gene->rec->map_u[i]->mid[4]);
        else
          gene->rec->nlos_u[i] += (map_sister->mid[4] - gene->rec->map_u[i]->mid[4]);
      }
    } // if (node is a duplication)
    else /* number of nodes, in sptree, between map and map's "children" ("-2" since difference in levels=1 means NO interm) */
      gene->rec->nlos_u[i] += (map_sister->mid[4] + map_up->mid[4] - 2 * gene->rec->map_u[i]->mid[4] - 2);
  }

  /* create virtual roots at every edge to calc the total number of dups (above + below) */
  for (i = 0; i < gene->nnodes; i++) if ((i != r_right) && (i != gene->root->id)) { /* only r_left (or r_right) is needed */
    map_root = mrca_between_nodes (species, gene->rec->map_u[i]->id, gene->rec->map_d[i]->id);
    /* duplications */ 
    thisdups = gene->rec->ndup_u[i] + gene->rec->ndup_d[i]; /* g_id now is the number of dups at virtual root */
    if ((map_root == gene->rec->map_u[i]) || (map_root == gene->rec->map_d[i])) thisdups++;
    /* losses */
    thisloss = gene->rec->nlos_u[i] + gene->rec->nlos_d[i]; /* number of losses at virtual root */
    if      ((map_root == gene->rec->map_u[i]) && (map_root != gene->rec->map_d[i]))
      thisloss += (gene->rec->map_d[i]->mid[4] - map_root->mid[4]); // "d(a(g),g) + 1" in SIAM.2000.729 
    else if ((map_root != gene->rec->map_u[i]) && (map_root == gene->rec->map_d[i]))
      thisloss += (gene->rec->map_u[i]->mid[4] - map_root->mid[4]); // "d(a(g),g) + 1" in SIAM.2000.729  
    else if ((map_root != gene->rec->map_u[i]) && (map_root != gene->rec->map_d[i]))
      thisloss += (gene->rec->map_u[i]->mid[4] + gene->rec->map_d[i]->mid[4] - 2 * map_root->mid[4] - 2);
    // else "loss(g) = 0" following SIAM.2000.729
    /* deep coalescences = loss - 2 x dups + 2 x |leaf difference between gene and species trees| */
    thiscoal = thisloss - 2 * thisdups + gene->rec->size_diff;

    if (thisdups < min_dups) min_dups = thisdups;
    if (thisloss < min_loss) min_loss = thisloss; 
    if (thiscoal < min_coal) min_coal = thiscoal;
  }
  gene->rec->ndups = min_dups;
  gene->rec->nloss = min_loss;
  gene->rec->ndcos = min_coal;
}

void
prepare_for_loss_calculation (topology gene, topology species)
{ 
  int i, c_l, c_r;

  for (i = 0; i < species->nleaves; i++) species->nodelist[i]->mid[2] = gene->rec->sp_count[i];
  for (i = 0; i < species->nleaves - 1; i++) { 
    /* mid[2] is the "effective" subtree cardinality: corrects for duplicates and absent species */
    c_l = species->postorder[i]->left->mid[2]; 
    c_r = species->postorder[i]->right->mid[2];
    species->postorder[i]->mid[2] = c_l + c_r; 

    /* mid[3] indicates if node is active or not (0 = pruned; 1 = normal; 0xffff = dummy node (only one child) */
    if ((!c_l) && (!c_r)) species->postorder[i]->mid[3] = 0; /* both children absent from gene: exclude this node */ 
    else if ((c_l) && (c_r)) species->postorder[i]->mid[3] = 1; /* both children present: this node is eligible */
    /* only one active child: this node will simply duplicate values from single valid child */
    //else species->postorder[i]->mid[3] = (c_l? 2: 3); /* id=2 if active child is left, id=3 if active is right node */
    else species->postorder[i]->mid[3] = 0xffff; 
  }

  /* in preorder: mid[4] has level (distance from root) taking into account only active species */
  if (species->root->mid[3] == 1) species->root->mid[4] = 0;
  else                            species->root->mid[4] = -1;

  for (i = species->nleaves-3; i >= 0; i--) {
    if (species->postorder[i]->mid[3] == 1)
      species->postorder[i]->mid[4] = species->postorder[i]->up->mid[4] + 1;
    else /* works if mid[3] is larger than one; if mid[3] is zero this node won't be mapped into anyway */
      species->postorder[i]->mid[4] = species->postorder[i]->up->mid[4];
  }
  for (i = 0; i < species->nleaves; i++) if (species->nodelist[i]->mid[2]) /* only leaves with one or more homologs */
    species->nodelist[i]->mid[4] = species->nodelist[i]->up->mid[4] + 1;
}

 /* OLD NOTE about DEEPCOAL (now we use nloss - 2xndups):
  * following ThanNahleh.PLoSComputBiol.2009.preprint, for "regular nodes" extra lineages above a node are equal to
  * |subtree rooted at node| - # coalescences below subtree (including root node of subtree) - 1 
  * No details are given for leaves, but on MolPhylEvol.1997.349 they mention the need for creating artificial nodes
  * for duplicated species and removal of subtrees with species unrepresented in gene family. For deep coal, we look
  * only at the species tree [1], where the mapping represents the coalescences. 
  * [1] the exception are maybe the species with more than one copy, but I didn't look at that. */ 
