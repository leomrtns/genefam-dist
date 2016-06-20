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

/*FIXME:: function declarations have wrong names; (sptree_split X RF_split or whatever -- now I'm concerned with splitspr) */

#include "topology_splitset.h"
#include "topology_mrca.h"

/*! \brief function used by qsort for a vector of bipartitions (from smaller to larger) */
int compare_splitset_bipartition_increasing (const void *a1, const void *a2);
/*! \brief function used by qsort for a vector of bipartitions (from larger to smaller) */
int compare_splitset_bipartition_decreasing (const void *a1, const void *a2);
/*! \brief copies bipartion information of internal nodes to splitset, "unroots" the bipartitions and order the resulting vectors */
void prepare_split_from_topologies (topology t1, topology t2, splitset split, int recycle_t1);
void prepare_genetree_sptree_split (topology gene, topology species, splitset split);
void split_add_gene_subtree (splitset split, int taxa);
int dSPR_topology_lowlevel (splitset split, int rf_or_hdist);
/*! \brief vector with pointers to bipartitions that are identical on both trees */
void split_create_agreement_list (splitset split);
void split_remove_agree_edges (splitset split, bipartition *b, int *nb);
void split_remove_duplicates (bipartition *b, int *nb);
void split_compress_agreement (splitset split);
void split_create_disagreement_list (splitset split);
void split_disagreement_assign_match (splitset split);
void split_find_small_disagreement (splitset split);
void split_remove_small_disagreement (splitset split);
void split_minimize_subtrees (splitset split);
void split_remove_redundant_bit (splitset split, int id);
void split_replace_bit (splitset split, int to, int from);
void split_new_size (splitset split, int size, bool update_bipartitions);
void bipvector_swap_position (bipartition *b, int i1, int i2);

/* For leafset-bijective trees AND gene/sptree RF/hdist calculation, we use splitset */

splitset 
new_splitset (int nleaves)
{
  splitset split;
  int i;

  split = (splitset) biomcmc_malloc (sizeof (struct splitset_struct));
  split->size = split->spsize = nleaves - 1; /* actually can be nleaves-3, but we use disagree[] as temporary for gene/species comparison */
  split->disagree = split->agree = NULL; /* used only when calculating dSPR/RF */
  split->n_agree = split->n_disagree = 0;
  split->prune = NULL; 
  split->h = NULL; /* hungarian method for bipartite matching; unused unless for tree distance */
  split->match = false;
  split->spr = split->spr_extra = split->rf = split->hdist = 0;

  split->g_split = (bipartition*) biomcmc_malloc (split->size * sizeof (bipartition));
  split->s_split = (bipartition*) biomcmc_malloc (split->size * sizeof (bipartition));
  split->g_split[0] = new_bipartition (nleaves);
  split->s_split[0] = new_bipartition (nleaves);
  for (i = 1; i < split->size; i++) {
    split->g_split[i] = new_bipartition_from_bipsize (split->g_split[0]->n); /* use same bipsize */
    split->s_split[i] = new_bipartition_from_bipsize (split->s_split[0]->n);
  }
  split->n_g = split->n_s = 0; /* number of elements actively in use */

  /* the pointer below is needed to store leaf info of gene tree on sptree; this pointer is used by del_splitset() */
  split->sp0   = split->s_split; /* used for gene/species RF calculation  (sp0 is longer than s_split) */
  return split;
}

splitset
new_splitset_RF (int nleaves)
{ /* nleaves should be relative to gene family */ 
  int i;
  splitset split = new_splitset (nleaves);

  split->agree = (bipartition*) biomcmc_malloc (split->size * sizeof (bipartition));
  split->agree[0] = new_bipartition (nleaves); // this bipsize will be recycled below 
  for (i = 1; i < split->size; i++) split->agree[i] = new_bipartition_from_bipsize (split->agree[0]->n);
  split->prune = new_bipartition_from_bipsize (split->agree[0]->n);
  split->h = new_hungarian (split->size);
  split->match = false; /* do we want to calculate the assignment matching cost (using hungarian() )? */

  return split;
}

void
del_splitset (splitset split)
{
  int i;
  if (!split) return;

  del_bipartition (split->prune);
  if (split->disagree) {
    for (i = split->size * split->size - 1; i >= 0; i--) del_bipartition (split->disagree[i]);
    free (split->disagree);
  }
  if (split->agree) {
    for (i = split->size - 1; i >= 0; i--) del_bipartition (split->agree[i]);
    free (split->agree);
  }
  if (split->g_split) {
    for (i = split->size - 1; i >= 0; i--) del_bipartition (split->g_split[i]);
    free (split->g_split);
  }
  if (split->sp0) {
    for (i = split->spsize - 1; i >= 0; i--) del_bipartition (split->sp0[i]);
    free (split->sp0);
  }
  del_hungarian (split->h);
  free (split);
}


int
compare_splitset_bipartition_increasing (const void *a1, const void *a2)
{ /* similar to bipartition_is_larger() */
  bipartition *b1 = (bipartition *) a1;
  bipartition *b2 = (bipartition *) a2;
  int i;

  if ((*b1)->n_ones > (*b2)->n_ones) return 1;
  if ((*b1)->n_ones < (*b2)->n_ones) return -1;

  for (i = (*b1)->n->ints - 1; (i >= 0) && ((*b1)->bs[i] == (*b2)->bs[i]); i--); /* find position of distinct bipartition elem*/
  if (i < 0) return 0; /* identical bipartitions */
  if ((*b1)->bs[i] > (*b2)->bs[i]) return 1;
  else return -1;
}

int
compare_splitset_bipartition_decreasing (const void *a1, const void *a2)
{ 
  return compare_splitset_bipartition_increasing (a2, a1);
}

void
prepare_split_from_topologies (topology t1, topology t2, splitset split, int recycle_t1)
{
  int i;
  if ((!recycle_t1) && (!t1->traversal_updated)) update_topology_traversal (t1); /* recycle_t1 > 0 --> we just used t1 in prev iteration */
  if (!t2->traversal_updated) update_topology_traversal (t2);
  /* the vector elements share a single bitstring size, that is modified by e.g. dSPR calculation */
  bipsize_resize (split->g_split[0]->n, split->g_split[0]->n->original_size); 
  bipsize_resize (split->s_split[0]->n, split->s_split[0]->n->original_size); 

  /* heavy child (more leaves or leaf with larger ID in case of tie) at left, after traversal */
  for (i=0; i < t1->nleaves-3; i++) {
    /* Q: why i from 0 to nleaves-4 ?       A: [nleaves-2] is root; [nleaves-3] is OR
     * left child of root  - in which case right child is leaf and left bipartition will be singleton; OR
     * right child of root - which will have same info as (previously visited) left child and therefore redundant */
    bipartition_copy (split->s_split[i], t2->postorder[i]->split);
    bipartition_flip_to_smaller_set (split->s_split[i]);
  }
  if (!recycle_t1) for (i=0; i < t1->nleaves-3; i++) {  /* Q: why i from 0 to nleaves-4 ? --> same as above */ 
    bipartition_copy (split->g_split[i], t1->postorder[i]->split);
    bipartition_flip_to_smaller_set (split->g_split[i]);
  }

  split->n_g = split->n_s = i;
  //for (i = 0; i < split->n_s; i++) bipartition_print_to_stdout (split->s_split[i]); printf ("::DEBUG::  :: ::S\n");
  if (!recycle_t1) qsort (split->g_split, split->n_g, sizeof (bipartition), compare_splitset_bipartition_increasing);
  qsort (split->s_split, split->n_s, sizeof (bipartition), compare_splitset_bipartition_increasing);
}

bool
topology_is_equal_unrooted (topology t1, topology t2, splitset split, int recycle_t1)
{
  int i;
  if ((split->size != t1->nleaves - 1) || (split->size != t2->nleaves - 1)) return false; /* can be an error, t1 or t2 has wrong size */

  prepare_split_from_topologies (t1, t2, split, recycle_t1);
  for (i=0; i < split->n_g; i++) if (!bipartition_is_equal (split->g_split[i], split->s_split[i])) return false;
  return true;
}

splitset
create_splitset_RF_genespecies (topology gene, topology species)
{
  int i, j; 
  splitset split = new_splitset_RF (gene->nleaves);

  /* gene->rec->sp_id[i] has species ID of each gene leaf i, and rec->sp_count[j] has multiplicity of each species j */
  if (!gene->rec) init_tree_recon_from_species_names (gene, species->taxlabel);

  /* species tree before removal of missing species, has at most 2*sp_nleaves-1 bipartitions, and after their removal _but_
   * considering subtrees spanned by gene tree (the mul-tree) has at most 2*gene_nleaves-1 bipartitions. We don't need
   * the leaves (trivial bipartitions) -- therefore the 2 x nleaves is overkill --, but sp0[] needs sp_nleaves. */
  if (gene->nleaves > species->nleaves) split->spsize = species->nleaves + gene->nleaves;
  else                                  split->spsize = 2 * species->nleaves; /* +sp leaves for the first elems */

  split->s_split = (bipartition*) biomcmc_realloc ((bipartition*) split->s_split, split->spsize * sizeof (bipartition));
  for (i = split->size; i < split->spsize; i++) split->s_split[i] = new_bipartition_from_bipsize (split->s_split[0]->n);
  split->sp0 = split->s_split; /* E.g.  sp0 points to [0|1|2|3|4|5|6|7] while s_split points to [4|5|6|7] */
  split->s_split = split->s_split + species->nleaves; /* pointer to nleaves-th element */

  for (i=0; i < species->nleaves; i++) { /* initialize bipartitions at leaves with gene leave information */
    bipartition_zero (split->sp0[i]);
    for (j = 0; j < gene->nleaves; j++) if (gene->rec->sp_id[j] == i) bipartition_set (split->sp0[i], j);
  }
  /* we could exclude species absent from gene tree (more memory efficient) [1], but it's easier to leave them since indexes
   * are preserved (when calling postorder[]->left->id for instance), otherwise we should need a index vector to map tree
   * IDs to positions here. Everything works fine here since bitstrings always refer to gene labels (that is, absent species are 
   * all zero or repeated bitstrings, when copying sptree structure).
   * [1] means swap_bipartition, s_split -= removed_leaves, realloc (sp0, smaller_size)  
   * */
  return split;
}

void
prepare_genetree_RF_split (topology gene, topology species, splitset split)
{
  int i;
  if (!gene->traversal_updated)    update_topology_traversal (gene);
  if (!species->traversal_updated) update_topology_traversal (species);

  bipsize_resize (split->g_split[0]->n, split->g_split[0]->n->original_size); 
  bipsize_resize (split->s_split[0]->n, split->s_split[0]->n->original_size); 

  /* use sp0[] that has info about the species leaves (coded for gene ids as a multifurcation) */ 
  for (i=0; i < species->nleaves - 1; i++) bipartition_zero (split->s_split[i]); /* some might be skipped below */
  for (i=0; i < species->nleaves - 1; i++) {
    bipartition_OR (split->sp0[species->postorder[i]->id], /* postorder[i]->id is always larger than nleaves */
                    split->sp0[species->postorder[i]->left->id], /* sp0 has elems not accessible by s_split (w/ leaves info) */
                    split->sp0[species->postorder[i]->right->id], false);
  }
  split->n_s = species->nleaves - 1;
  /* some spleaves may be zero (absent on gene) leading to internal nodes with < 2 elems, thus we need to "minimize" the species subtrees 
   * (not to mention some other skipped node in postorder) */ 
  for (i = 0; i < split->n_s; i++) { /* here we use s_split[], that points to internal bipartitions only */
    bipartition_flip_to_smaller_set (split->s_split[i]);
    if (split->s_split[i]->n_ones < 2) { split->n_s--; bipvector_swap_position (split->s_split, i, split->n_s); i--; }
  }
  split_remove_duplicates (split->s_split, &(split->n_s));
  for (i=0; i < species->nleaves; i++) if (split->sp0[i]->n_ones > 1) { /* add multifurcations representing species on gene tree */
    bipartition_copy (split->s_split[split->n_s], split->sp0[i]); /* idea from RF distance (arxiv.1210.2665) */
    bipartition_flip_to_smaller_set (split->s_split[split->n_s++]); /* this species has many many copies, that is most of leaves are same species... */
  }
  /* gene tree bipartitions are simpler */
  split->n_g = gene->nleaves - 3;
  for (i=0; i < split->n_g; i++) {
    bipartition_copy (split->g_split[i], gene->postorder[i]->split);
    bipartition_flip_to_smaller_set (split->g_split[i]);
  }
}

/* empirical observation: split->spr below tends to overestimate (very rarely it overestimates), while 
 * spr + spr_extra may overestimate, have a higher variability but has the impressive property of recognizing even a
 * lot of SPRs. A good compromise seems to weight the contribution of spr_extra. BTW this variable counts the number
 * of "swapped" prune edges -- that is, a prune bipartition representing subtrees that are in opposite sides of the
 * edge, that usually mean 2 SPRs like (A,B)-+-(C,D)  <-> (A,C)-+-(B,D) */  
/* split->spr += (split->spr_extra/2); // I'll comment it out and leave it to the calling function */

int
dSPR_topology (topology t1, topology t2, splitset split)
{
  split->hdist = split->spr = split->spr_extra = split->rf = 0;
  if (topology_is_equal_unrooted (t1, t2, split, false)) return 0;
  split->match = true;  return dSPR_topology_lowlevel (split, 0);
}

int
dSPR_topology_rf (topology t1, topology t2, splitset split)
{
  split->hdist = split->spr = split->spr_extra = split->rf = 0;
  if (topology_is_equal_unrooted (t1, t2, split, false)) return 0;
  split->match = true;  return dSPR_topology_lowlevel (split, 1);
}

int
dSPR_topology_hdist (topology t1, topology t2, splitset split)
{
  split->hdist = split->spr = split->spr_extra = split->rf = 0;
  if (topology_is_equal_unrooted (t1, t2, split, false)) return 0;
  split->match = true;  return dSPR_topology_lowlevel (split, 2);
}

int
RF_gene_species (topology gene, topology species, splitset split)
{
  split->hdist = split->spr = split->spr_extra = split->rf = 0;
  prepare_genetree_RF_split (gene, species, split);
  split->match = true;  return dSPR_topology_lowlevel (split, 1);
}

int
RF_gene_species_hdist (topology gene, topology species, splitset split)
{
  split->hdist = split->spr = split->spr_extra = split->rf = 0;
  prepare_genetree_RF_split (gene, species, split);
  split->match = true;  return dSPR_topology_lowlevel (split, 2);
}

int
RF_gene_species_hdist (topology gene, topology species, splitset split)
{ /* Not a good approximation if no bijection gene/sptree. In such a case, try splitspr version */
  split->hdist = split->spr = split->spr_extra = split->rf = 0;
  prepare_genetree_RF_split (gene, species, split);
  split->match = true;  return dSPR_topology_lowlevel (split, 2);
}

int
dSPR_topology_lowlevel (splitset split, int rf_or_hdist)
{
  int i = 0, mismatch = -1;

  split->hdist = split->spr = split->spr_extra = split->rf = 0;
  split->n_agree = split->n_disagree = 0;
  bipsize_resize (split->disagree[0]->n, split->g_split[0]->n->bits); 
  bipsize_resize (split->agree[0]->n,    split->g_split[0]->n->bits); 
  //for (i = 0; i < split->n_g; i++) bipartition_print_to_stdout (split->g_split[i]); printf ("G   ::DEBUG 0 ::\n");
  //for (i = 0; i < split->n_s; i++) bipartition_print_to_stdout (split->s_split[i]); printf ("S\n");

  i++; /* to trick -Werror, since we don't use it unless for debug */
  while (mismatch) {
    split_create_agreement_list (split);  // vector of identical bipartitions
    split_compress_agreement (split);     // iterative replacement of cherry by new leaf

//    for (i = 0; i < split->n_g; i++) bipartition_print_to_stdout (split->g_split[i]); printf ("G   ::DEBUG::\n");
//    for (i = 0; i < split->n_s; i++) bipartition_print_to_stdout (split->s_split[i]); printf ("S\n");
//    for (i = 0; i < split->n_agree; i++) bipartition_print_to_stdout (split->agree[i]); printf ("A\n");

    if (mismatch == -1) split->rf = split->n_g + split->n_s;
    if (rf_or_hdist == 1) return split->rf; /* do not calculate SPR, exit now  (1=RF, 2=hdist) */
    mismatch = (split->n_g > 0) && (split->n_s > 0); // all edges were in agreement
    if (!mismatch) return split->spr;
    
    split_create_disagreement_list (split); // vector of smallest disagreements
    split_disagreement_assign_match (split); /* assignment matching between edges using hungarian method (split->hdist after first time) */
    if (rf_or_hdist == 2) return split->hdist; /* do not calculate SPR, exit now */
    
    split_remove_duplicates (split->disagree, &(split->n_disagree)); // some elements are equal; this function also qsorts 
    split_find_small_disagreement (split);  // could also be one leaf only 
    
    //for (i = 0; i < split->n_disagree; i++) { bipartition_print_to_stdout (split->disagree[i]); printf ("\n"); }
    //printf ("{%d} prune: ", split->n_disagree); bipartition_print_to_stdout (split->prune); printf ("\n");

    split->spr++;
    split_remove_small_disagreement (split);

    split_minimize_subtrees (split);
    mismatch = (split->n_g > 0) && (split->n_s > 0); // all edges were in agreement
  }
  return split->spr;
}

void
split_create_agreement_list (splitset split)
{
  int s, g;
  for (g = 0; g < split->n_g; g++) for (s = 0; s < split->n_s; s++)
    if (bipartition_is_equal (split->g_split[g], split->s_split[s])) {
      bipartition_copy (split->agree[split->n_agree++], split->g_split[g]); 
      split->n_g--; bipvector_swap_position (split->g_split, g, split->n_g); /* if we don't swap them, we lose ref to "old" value on g_split[] */
      split->n_s--; bipvector_swap_position (split->s_split, s, split->n_s); 
      g--; s = split->n_s; /* pretend loop finished, examine again with new values */
    }
  split_remove_agree_edges (split, split->g_split, &(split->n_g));
  split_remove_agree_edges (split, split->s_split, &(split->n_s));
}

void
split_remove_agree_edges (splitset split, bipartition *b, int *nb)
{
  int i, a;
  for (i = 0; i < (*nb); i++) for (a = 0; a < split->n_agree; a++) 
    if (bipartition_is_equal (b[i], split->agree[a])) {
      (*nb)--; 
      bipvector_swap_position (b, i, (*nb));
      i--;
      a = split->n_agree; /* loop again over new value */
    }
}

void
split_remove_duplicates (bipartition *b, int *nb)
{
  int i, j;
  bipartition pivot;
  if ((*nb) < 2) return; /* only if we have a vector with > 1 element */
  qsort (b, (*nb), sizeof (bipartition), compare_splitset_bipartition_increasing);

  for (i = (*nb) - 1; i >= 1; i--)
    if (bipartition_is_equal (b[i], b[i-1])) {
      pivot = b[i]; /* do not lose a pointer to this element */
      for (j = i; j < (*nb)-1; j++) b[j] = b[j+1];
      b[j] = pivot; /* j = (*nb) - 1, which will become obsolete through next line -->  (*nb)-- */
      (*nb)--;
    }
}

void
split_compress_agreement (splitset split)
{
  int i, j, pair[2];

  for (i = 0; i < split->n_agree; i++) if (split->agree[i]->n_ones == 2) { /* cherry in common, can be represented by just one leaf */
    bipartition_to_int_vector (split->agree[i], pair, 2);
    split_remove_redundant_bit (split, pair[1]);
    split_new_size (split,split->agree[0]->n->bits - 1, false); /* false = do not recalculate every bipartition's last elem */
    bipartition_resize_vector (split->agree, split->n_agree);
    for (j = 0; j < split->n_agree; j++) { /* minimize subtree size and remove single leaves */
      bipartition_flip_to_smaller_set (split->agree[j]); /* agree only */
      if (split->agree[j]->n_ones < 2) bipvector_swap_position (split->agree, j--, --split->n_agree);
    }
    i = -1; /* redo all iterations, with new info (agree[] will be smaller) */
  }
  bipartition_resize_vector (split->g_split, split->n_g);
  bipartition_resize_vector (split->s_split, split->n_s);
}

void
split_create_disagreement_list (splitset split)
{
  int g, s;

  for (g = 0; g < split->n_g; g++) for (s = 0; s < split->n_s; s++) { 
    //bipartition_XOR (split->disagree[split->n_disagree], split->g_split[g], split->s_split[s], true); /* true means to calculate n_ones */
    //bipartition_flip_to_smaller_set (split->disagree[split->n_disagree++]);
    bipartition_XOR (split->disagree[g * split->n_s + s], split->g_split[g], split->s_split[s], true); /* true means to calculate n_ones */
    bipartition_flip_to_smaller_set (split->disagree[g * split->n_s + s]);
  }
  split->n_disagree = split->n_g * split->n_s;
}

void
split_disagreement_assign_match (splitset split)
{ /* also calculates split->hdist */ 
  int g, s, max_n, sum = 0;
  
  if (split->n_g > split->n_s) max_n = split->n_g;
  else                         max_n = split->n_s;
  if (max_n < 2) return;
  
  hungarian_reset (split->h);
  for (g = 0; g < split->n_g; g++) for (s = 0; s < split->n_s; s++)  
    hungarian_update_cost (split->h, g, s, split->disagree[g * split->n_s + s]->n_ones);
  hungarian_solve (split->h, max_n);
  /* now split->h->col_mate will have the pairs */
  /* if we do the matching below it becomes much faster, but we may miss the best prune subtrees in a few cases (do not compromise the algo) */
  split->n_disagree = 0;
  for (g = 0; g < max_n; g++) if ((g < split->n_g) && ( split->h->col_mate[g] < split->n_s)) { /* some matchings might be to dummy edges */
    bipartition_XOR (split->disagree[split->n_disagree], split->g_split[g], split->s_split[split->h->col_mate[g]], true); /* true means to calculate n_ones */
    bipartition_flip_to_smaller_set (split->disagree[split->n_disagree++]);
    sum += split->disagree[split->n_disagree-1]->n_ones;
  }
//  if (split->match) { split->hdist = (100 * split->h->final_cost)/split->h->initial_cost; split->match = false; }
  if (split->match) { split->hdist = split->h->final_cost+split->h->initial_cost; split->match = false; }
}

void
split_find_small_disagreement (splitset split)
{
  bipartition dis;
  int a, d;

  bipartition_copy (split->prune, split->disagree[0]); /* smallest, in case we don't find a better one in loop below */ 
  if (split->prune->n_ones < 2) return;

  dis = new_bipartition_from_bipsize (split->disagree[0]->n);
  for (d = 0; d < split->n_disagree; d++) for (a = 0; a < split->n_agree; a++) {
    if ((split->disagree[d]->n_ones == split->agree[a]->n_ones) || 
        (split->disagree[d]->n_ones == (split->agree[a]->n->bits - split->agree[a]->n_ones))) {
      bipartition_XOR (dis, split->disagree[d], split->agree[a], true);
      if      (!dis->n_ones)               { bipartition_copy (split->prune, split->disagree[d]); d = split->n_disagree; a = split->n_agree; }
      else if (dis->n_ones == dis->n->bits) { bipartition_NOT (split->prune, split->disagree[d]); d = split->n_disagree; a = split->n_agree; }
      //      bipartition_XORNOT (dis, split->disagree[d], split->agree[a], true); /* equiv to above XOR + NOT */
      //      if (!dis->n_ones) { bipartition_NOT (split->prune, split->disagree[d]); del_bipartition (dis); return; }
    }
  }
  /* check if prune nodes are all on same side of a tree or if they are actually two SPRs (one from each tree) */
  for (d = 0; d < split->n_g; d++) {
    if (!bipartition_contains_bits (split->g_split[d], split->prune)) {
      bipartition_NOT (dis, split->g_split[d]);
      if (!bipartition_contains_bits (dis, split->prune)) { split->spr_extra++; d = split->n_g; }
    }
  } 
  del_bipartition (dis); 
}

void
split_remove_small_disagreement (splitset split)
{
  int *index, i, j = split->prune->n_ones - 1, k = 0, size = split->agree[0]->n->bits; 

  index = (int*) biomcmc_malloc (split->prune->n_ones * sizeof (int));
  bipartition_to_int_vector (split->prune, index, split->prune->n_ones);

  for (i = size - 1; i >= (size - split->prune->n_ones); i--) {
    if (index[k] >= (size - split->prune->n_ones)) i = -1;
    else {
      if (i == index[j]) j--;
      else split_replace_bit (split, index[k++], i); 
    }
  }

  split_new_size (split,size - split->prune->n_ones, true); 
  if (index) free (index);
}

void
split_minimize_subtrees (splitset split)
{
  int i;

  for (i = 0; i < split->n_s; i++) {      
    bipartition_flip_to_smaller_set (split->s_split[i]);
    if (split->s_split[i]->n_ones < 2) { split->n_s--; bipvector_swap_position (split->s_split, i, split->n_s); i--; }
  }
  for (i = 0; i < split->n_g; i++) {
    bipartition_flip_to_smaller_set (split->g_split[i]);
    if (split->g_split[i]->n_ones < 2) { split->n_g--; bipvector_swap_position (split->g_split, i, split->n_g); i--; }
  }
  for (i = 0; i < split->n_agree; i++) {
    bipartition_flip_to_smaller_set (split->agree[i]);
    if (split->agree[i]->n_ones < 2) { split->n_agree--; bipvector_swap_position (split->agree, i, split->n_agree); i--; }
  }
}

void
split_remove_redundant_bit (splitset split, int id)
{
  int  last = split->agree[0]->n->bits-1;
  if (id < last) split_replace_bit (split, id, last);
}

void
split_replace_bit (splitset split, int to, int from)
{
  if (from <= to) return;
  /*not needed for disagree[] */
  bipartition_replace_bit_in_vector (split->agree,   split->n_agree, to, from, true);
  bipartition_replace_bit_in_vector (split->g_split, split->n_g,     to, from, true);
  bipartition_replace_bit_in_vector (split->s_split, split->n_s,     to, from, true);
}

void
split_new_size (splitset split, int size, bool update_bipartitions)
{
  bipsize_resize (split->g_split[0]->n, size);
  bipsize_resize (split->s_split[0]->n, size);
  bipsize_resize (split->agree[0]->n, size);
  bipsize_resize (split->disagree[0]->n, size);
  if (update_bipartitions) {
    bipartition_resize_vector (split->g_split, split->n_g);
    bipartition_resize_vector (split->s_split, split->n_s);
    bipartition_resize_vector (split->agree, split->n_agree);
  }
}

void
bipvector_swap_position (bipartition *b, int i1, int i2)
{
  bipartition pivot = b[i1];
  b[i1] = b[i2];
  b[i2] = pivot;
}

/* For gene/species trees dSPR calculation, we use splitspr instead of splitset */
#ifdef _THIS_IS_COMMENTED_OUT_SINCE_UNFINISHED
sprsplit
new_sprsplit (topology gene, topology species)
{
  sprsplit spl;
  int i, min_size, max_size;

  spl = (sprsplit) biomcmc_malloc (sizeof (struct sprsplit_struct));
  spl->gsize  = gene->nleaves - 3;
  spl->spsize = species->nleaves - 3;
  spl->spr = split->spr_extra = split->rf = split->hdist = 0;

  spl->gene = (bipartition*) biomcmc_malloc (spl->gsize  * sizeof (bipartition));
  spl->g_sp = (bipartition*) biomcmc_malloc (2 * spl->gsize  * sizeof (bipartition)); /* gene, but with sp info */
  spl->sp   = (bipartition*) biomcmc_malloc (spl->spsize * sizeof (bipartition));
  spl->gene[0] = new_bipartition (gene->nleaves); /* info from species on each side */
  spl->g_sp[0] = new_bipartition (species->nleaves);
  spl->sp[0]   = new_bipartition (species->nleaves);
  for (i = 1; i < spl->gsize;  i++)     spl->gene[i] = new_bipartition_from_bipsize (spl->gene[0]->n); /* use same bipsize */
  for (i = 1; i < 2 * spl->gsize;  i++) spl->g_sp[i] = new_bipartition_from_bipsize (spl->g_sp[0]->n); /* use same bipsize */
  for (i = 1; i < spl->spsize; i++) spl->sp[i]   = new_bipartition_from_bipsize (spl->sp[0]->n); /* use same bipsize */
  spl->n_g = spl->n_s = 0; /* number of elements actively in use */

  if (gene->nleaves > species->nleaves) { min_size = species->nleaves - 1; max_size = gene->nleaves - 1; } 
  else                                  { max_size = species->nleaves - 1; min_size = gene->nleaves - 1; } 
  spl->asize = min_size; /* agree[] vector */
  spl->dsize = (gene->nleaves - 1) * (species->npeaves - 1); /* disagree[] vector */
  spl->h = new_hungarian (max_size);
  spl->match = false; /* do we want to calculate the assignment matching cost (using hungarian() )? */

  spl->agree    = (bipartition*) biomcmc_malloc (spl->asize * sizeof (bipartition));
  spl->disagree = (bipartition*) biomcmc_malloc (spl->dsize * sizeof (bipartition));
  spl->agree[0]    = new_bipartition (species->nleaves); // this bipsize will be recycled below 
  spl->disagree[0] = new_bipartition (species->nleaves); 
  for (i = 1; i < spl->asize; i++) spl->agree[i]    = new_bipartition_from_bipsize (spl->agree[0]->n);
  for (i = 1; i < spl->dsize; i++) spl->disagree[i] = new_bipartition_from_bipsize (spl->disagree[0]->n);
  spl->prune     = new_bipartition_from_bipsize (split->gene[0]->n);
  spl->spcounter = new_bipartition_from_bipsize (split->sp[0]->n);

  spl->map = (bipartition*) biomcmc_malloc (species->nleaves * sizeof (bipartition));
  spl->map[0] = new_bipartition (gene->nleaves); 
  for (i = 1; i < species->nleaves; i++) spl->map[i] = new_bipartition_from_bipsize (spl->map[0]->n);

  return spl;
}

void
del_sprsplit (sprsplit spl)
{
  int i;
  if (!spl) return;

  if (spl->map) {
    for (i = spl->spcounter->n->original_size - 1; i >= 0; i--) del_bipartition (spl->map[i]);
    free (spl->map);
  }
  if (spl->disagree) {
    for (i = spl->dsize - 1; i >= 0; i--) del_bipartition (spl->disagree[i]);
    free (sp->disagree);
  }
  if (spl->agree) {
    for (i = spl->asize - 1; i >= 0; i--) del_bipartition (spl->agree[i]);
    free (spl->agree);
  }
  if (spl->gene) {
    for (i = spl->gsize - 1; i >= 0; i--) del_bipartition (spl->gene[i]);
    free (spl->gene);
  }
  if (spl->g_sp) {
    for (i = 2 * spl->gsize - 1; i >= 0; i--) del_bipartition (spl->g_sp[i]);
    free (spl->g_sp);
  }
  if (spl->sp) {
    for (i = spl->spsize - 1; i >= 0; i--) del_bipartition (spl->sp[i]);
    free (spl->sp);
  }
  del_bipartition (spl->prune);
  del_bipartition (spl->spcounter);
  del_hungarian (spl->h);
  free (spl);
}

sprsplit
create_splitspr_genespecies (topology gene, topology species)
{
  int i, j;  //STOPHERE
  sprsplit split = new_splitspr (gene, species);
  /* gene->rec->sp_id[i] has species ID of each gene leaf i, and rec->sp_count[j] has multiplicity of each species j */
  if (!gene->rec) init_tree_recon_from_species_names (gene, species->taxlabel);
 /* initialize bipartitions at leaves with species information */
  for (i=0; i < gene->nleaves; i++) bipartition_set (split->gene[i], gene->rec->sp_id[i]);
  return split;
}

void
prepare_genetree_dSPR_split (topology gene, topology species, splitspr s)
{
  int i, leap = gene->nleaves - 1;
  if (!gene->traversal_updated)    update_topology_traversal (gene);
  if (!species->traversal_updated) update_topology_traversal (species);
/* TODO: the translation is done through one mask for each species, e.g. mask 11000 if gene leaves 0 and 1 belong to this species */
  bipsize_resize (s->gene[0]->n, s->gene[0]->n->original_size); 
  bipsize_resize (s->sp[0]->n,   s->sp[0]->n->original_size); 
  for (i = gene->nleaves; i < s->gsize; i++) bipartition_zero (s->gene[i]);
  /* we need both bipartitions below and above for each gene node -- we cannot use ~bipartition since same species on both sides 
   * below(node) = below(left) + below(right); above(node) = above(up) + below (sister) */
  for (i=0; i < gene->nleaves - 2; i++) { /* postorder[nleaves - 2] => root */
    bipartition_OR (s->gene[gene->postorder[i]->id], /* postorder[i]->id is always larger than nleaves */
                    s->gene[gene->postorder[i]->left->id], /* gene has elems not accessible by g1 and g2 (w/ leaves info) */
                    s->gene[gene->postorder[i]->right->id], false);
  }
  /* root node will have all species below and none above (due to bipartition_zero() ) */
  for (i=gene->nleaves - 4; i >= 0; i--) { /* now we go in preorder, using bipartition 'below' to update 'above' */
    bipartition_OR (s->gene[leap + gene->postorder[i]->id], /* "+leap" to store info 'above' node */
                    s->gene[leap + gene->postorder[i]->up->id], /* leap+root_id will have zero */
                    s->gene[gene->postorder[i]->sister->id], false); /* this is "below" node, no leap */
  }
  s->g1 = s->gene + gene->nleaves; /* points to first internal node */
  s->g2 = s->g1 + leap; /* points to first internal node */
  s->n_g = gene->nleaves - 1;
  for (i=0; i < s->n_g; i++) {
    if ((!s->g2[i]->n_ones) && ((!s->g1[i]->n_ones) || (s->g1[i]->n_ones == s->g1[i]->n->bits))) {
      s->n_g--; /* remove root and unused nodes */
      bipvector_swap_position (s->s1, i, s->n_g);
      bipvector_swap_position (s->s2, i, s->n_g); 
    }
  }
  /* sptree bipartitions are simpler (for why nleaves-3, see comment on prepare_split_from_topologies(): one is root and other is redundant/singleton) */
  s->n_s = species->nleaves - 3;
  for (i=0; i < s->n_s; i++) {
    bipartition_copy (s->sp[i], species->postorder[i]->split);
    bipartition_flip_to_smaller_set (s->sp[i]);
  }
  s->hdist = s->spr = s->spr_extra = s->rf = 0;
}

int
spr_gene_species (topology gene, topology species, splitspr s)
{
  prepare_genetree_dSPR_split (gene, species, s);
  s->match = true;  return spr_topology_genespecies_lowlevel (s, 0);
}

int
spr_gene_species_hdist (topology gene, topology species, splitspr s)
{
  prepare_genetree_dSPR_split (gene, species, s);
  split->match = true;  return spr_topology_genespecies_lowlevel (s, 2);
}

int
spr_topology_genespecies_lowlevel (splitset split, int rf_or_hdist)
{
  int i = 0, mismatch = -1;

  split->hdist = split->spr = split->spr_extra = split->rf = 0;
  split->n_agree = split->n_disagree = 0;
  bipsize_resize (split->disagree[0]->n, split->g_split[0]->n->bits); 
  bipsize_resize (split->agree[0]->n,    split->g_split[0]->n->bits); 
  //for (i = 0; i < split->n_g; i++) bipartition_print_to_stdout (split->g_split[i]); printf ("G   ::DEBUG 0 ::\n");
  //for (i = 0; i < split->n_s; i++) bipartition_print_to_stdout (split->s_split[i]); printf ("S\n");

  i++; /* to trick -Werror, since we don't use it unless for debug */
  while (mismatch) {
    split_create_agreement_list (split);  // vector of identical bipartitions
    split_compress_agreement (split);     // iterative replacement of cherry by new leaf

//    for (i = 0; i < split->n_g; i++) bipartition_print_to_stdout (split->g_split[i]); printf ("G   ::DEBUG::\n");
//    for (i = 0; i < split->n_s; i++) bipartition_print_to_stdout (split->s_split[i]); printf ("S\n");
//    for (i = 0; i < split->n_agree; i++) bipartition_print_to_stdout (split->agree[i]); printf ("A\n");

    if (mismatch == -1) split->rf = split->n_g + split->n_s;
    if (rf_or_hdist == 1) return split->rf; /* do not calculate SPR, exit now  (1=RF, 2=hdist) */
    mismatch = (split->n_g > 0) && (split->n_s > 0); // all edges were in agreement
    if (!mismatch) return split->spr;
    
    split_create_disagreement_list (split); // vector of smallest disagreements
    split_disagreement_assign_match (split); /* assignment matching between edges using hungarian method (split->hdist after first time) */
    if (rf_or_hdist == 2) return split->hdist; /* do not calculate SPR, exit now */
    
    split_remove_duplicates (split->disagree, &(split->n_disagree)); // some elements are equal; this function also qsorts 
    split_find_small_disagreement (split);  // could also be one leaf only 
    
    //for (i = 0; i < split->n_disagree; i++) { bipartition_print_to_stdout (split->disagree[i]); printf ("\n"); }
    //printf ("{%d} prune: ", split->n_disagree); bipartition_print_to_stdout (split->prune); printf ("\n");

    split->spr++;
    split_remove_small_disagreement (split);

    split_minimize_subtrees (split);
    mismatch = (split->n_g > 0) && (split->n_s > 0); // all edges were in agreement
  }
  return split->spr;
}

void
split_create_agreement_list (splitset split)
{
  int s, g;
  for (g = 0; g < split->n_g; g++) for (s = 0; s < split->n_s; s++)
    if (bipartition_is_equal (split->g_split[g], split->s_split[s])) {
      bipartition_copy (split->agree[split->n_agree++], split->g_split[g]); 
      split->n_g--; bipvector_swap_position (split->g_split, g, split->n_g); /* if we don't swap them, we lose ref to "old" value on g_split[] */
      split->n_s--; bipvector_swap_position (split->s_split, s, split->n_s); 
      g--; s = split->n_s; /* pretend loop finished, examine again with new values */
    }
  split_remove_agree_edges (split, split->g_split, &(split->n_g));
  split_remove_agree_edges (split, split->s_split, &(split->n_s));
}

void
split_remove_agree_edges (splitset split, bipartition *b, int *nb)
{
  int i, a;
  for (i = 0; i < (*nb); i++) for (a = 0; a < split->n_agree; a++) 
    if (bipartition_is_equal (b[i], split->agree[a])) {
      (*nb)--; 
      bipvector_swap_position (b, i, (*nb));
      i--;
      a = split->n_agree; /* loop again over new value */
    }
}

void
split_remove_duplicates (bipartition *b, int *nb)
{
  int i, j;
  bipartition pivot;
  if ((*nb) < 2) return; /* only if we have a vector with > 1 element */
  qsort (b, (*nb), sizeof (bipartition), compare_splitset_bipartition_increasing);

  for (i = (*nb) - 1; i >= 1; i--)
    if (bipartition_is_equal (b[i], b[i-1])) {
      pivot = b[i]; /* do not lose a pointer to this element */
      for (j = i; j < (*nb)-1; j++) b[j] = b[j+1];
      b[j] = pivot; /* j = (*nb) - 1, which will become obsolete through next line -->  (*nb)-- */
      (*nb)--;
    }
}

void
split_compress_agreement (splitset split)
{
  int i, j, pair[2];

  for (i = 0; i < split->n_agree; i++) if (split->agree[i]->n_ones == 2) { /* cherry in common, can be represented by just one leaf */
    bipartition_to_int_vector (split->agree[i], pair, 2);
    split_remove_redundant_bit (split, pair[1]);
    split_new_size (split,split->agree[0]->n->bits - 1, false); /* false = do not recalculate every bipartition's last elem */
    bipartition_resize_vector (split->agree, split->n_agree);
    for (j = 0; j < split->n_agree; j++) { /* minimize subtree size and remove single leaves */
      bipartition_flip_to_smaller_set (split->agree[j]); /* agree only */
      if (split->agree[j]->n_ones < 2) bipvector_swap_position (split->agree, j--, --split->n_agree);
    }
    i = -1; /* redo all iterations, with new info (agree[] will be smaller) */
  }
  bipartition_resize_vector (split->g_split, split->n_g);
  bipartition_resize_vector (split->s_split, split->n_s);
}

void
split_create_disagreement_list (splitset split)
{
  int g, s;

  for (g = 0; g < split->n_g; g++) for (s = 0; s < split->n_s; s++) { 
    //bipartition_XOR (split->disagree[split->n_disagree], split->g_split[g], split->s_split[s], true); /* true means to calculate n_ones */
    //bipartition_flip_to_smaller_set (split->disagree[split->n_disagree++]);
    bipartition_XOR (split->disagree[g * split->n_s + s], split->g_split[g], split->s_split[s], true); /* true means to calculate n_ones */
    bipartition_flip_to_smaller_set (split->disagree[g * split->n_s + s]);
  }
  split->n_disagree = split->n_g * split->n_s;
}

void
split_disagreement_assign_match (splitset split)
{ /* also calculates split->hdist */ 
  int g, s, max_n, sum = 0;
  
  if (split->n_g > split->n_s) max_n = split->n_g;
  else                         max_n = split->n_s;
  if (max_n < 2) return;
  
  hungarian_reset (split->h);
  for (g = 0; g < split->n_g; g++) for (s = 0; s < split->n_s; s++)  
    hungarian_update_cost (split->h, g, s, split->disagree[g * split->n_s + s]->n_ones);
  hungarian_solve (split->h, max_n);
  /* now split->h->col_mate will have the pairs */
  /* if we do the matching below it becomes much faster, but we may miss the best prune subtrees in a few cases (do not compromise the algo) */
  split->n_disagree = 0;
  for (g = 0; g < max_n; g++) if ((g < split->n_g) && ( split->h->col_mate[g] < split->n_s)) { /* some matchings might be to dummy edges */
    bipartition_XOR (split->disagree[split->n_disagree], split->g_split[g], split->s_split[split->h->col_mate[g]], true); /* true means to calculate n_ones */
    bipartition_flip_to_smaller_set (split->disagree[split->n_disagree++]);
    sum += split->disagree[split->n_disagree-1]->n_ones;
  }
//  if (split->match) { split->hdist = (100 * split->h->final_cost)/split->h->initial_cost; split->match = false; }
  if (split->match) { split->hdist = split->h->final_cost+split->h->initial_cost; split->match = false; }
}

void
split_find_small_disagreement (splitset split)
{
  bipartition dis;
  int a, d;

  bipartition_copy (split->prune, split->disagree[0]); /* smallest, in case we don't find a better one in loop below */ 
  if (split->prune->n_ones < 2) return;

  dis = new_bipartition_from_bipsize (split->disagree[0]->n);
  for (d = 0; d < split->n_disagree; d++) for (a = 0; a < split->n_agree; a++) {
    if ((split->disagree[d]->n_ones == split->agree[a]->n_ones) || 
        (split->disagree[d]->n_ones == (split->agree[a]->n->bits - split->agree[a]->n_ones))) {
      bipartition_XOR (dis, split->disagree[d], split->agree[a], true);
      if      (!dis->n_ones)               { bipartition_copy (split->prune, split->disagree[d]); d = split->n_disagree; a = split->n_agree; }
      else if (dis->n_ones == dis->n->bits) { bipartition_NOT (split->prune, split->disagree[d]); d = split->n_disagree; a = split->n_agree; }
      //      bipartition_XORNOT (dis, split->disagree[d], split->agree[a], true); /* equiv to above XOR + NOT */
      //      if (!dis->n_ones) { bipartition_NOT (split->prune, split->disagree[d]); del_bipartition (dis); return; }
    }
  }
  /* check if prune nodes are all on same side of a tree or if they are actually two SPRs (one from each tree) */
  for (d = 0; d < split->n_g; d++) {
    if (!bipartition_contains_bits (split->g_split[d], split->prune)) {
      bipartition_NOT (dis, split->g_split[d]);
      if (!bipartition_contains_bits (dis, split->prune)) { split->spr_extra++; d = split->n_g; }
    }
  } 
  del_bipartition (dis); 
}

void
split_remove_small_disagreement (splitset split)
{
  int *index, i, j = split->prune->n_ones - 1, k = 0, size = split->agree[0]->n->bits; 

  index = (int*) biomcmc_malloc (split->prune->n_ones * sizeof (int));
  bipartition_to_int_vector (split->prune, index, split->prune->n_ones);

  for (i = size - 1; i >= (size - split->prune->n_ones); i--) {
    if (index[k] >= (size - split->prune->n_ones)) i = -1;
    else {
      if (i == index[j]) j--;
      else split_replace_bit (split, index[k++], i); 
    }
  }

  split_new_size (split,size - split->prune->n_ones, true); 
  if (index) free (index);
}

void
split_minimize_subtrees (splitset split)
{
  int i;

  for (i = 0; i < split->n_s; i++) {      
    bipartition_flip_to_smaller_set (split->s_split[i]);
    if (split->s_split[i]->n_ones < 2) { split->n_s--; bipvector_swap_position (split->s_split, i, split->n_s); i--; }
  }
  for (i = 0; i < split->n_g; i++) {
    bipartition_flip_to_smaller_set (split->g_split[i]);
    if (split->g_split[i]->n_ones < 2) { split->n_g--; bipvector_swap_position (split->g_split, i, split->n_g); i--; }
  }
  for (i = 0; i < split->n_agree; i++) {
    bipartition_flip_to_smaller_set (split->agree[i]);
    if (split->agree[i]->n_ones < 2) { split->n_agree--; bipvector_swap_position (split->agree, i, split->n_agree); i--; }
  }
}

void
split_remove_redundant_bit (splitset split, int id)
{
  int  last = split->agree[0]->n->bits-1;
  if (id < last) split_replace_bit (split, id, last);
}

void
split_replace_bit (splitset split, int to, int from)
{
  if (from <= to) return;
  /*not needed for disagree[] */
  bipartition_replace_bit_in_vector (split->agree,   split->n_agree, to, from, true);
  bipartition_replace_bit_in_vector (split->g_split, split->n_g,     to, from, true);
  bipartition_replace_bit_in_vector (split->s_split, split->n_s,     to, from, true);
}

void
split_new_size (splitset split, int size, bool update_bipartitions)
{
  bipsize_resize (split->g_split[0]->n, size);
  bipsize_resize (split->s_split[0]->n, size);
  bipsize_resize (split->agree[0]->n, size);
  bipsize_resize (split->disagree[0]->n, size);
  if (update_bipartitions) {
    bipartition_resize_vector (split->g_split, split->n_g);
    bipartition_resize_vector (split->s_split, split->n_s);
    bipartition_resize_vector (split->agree, split->n_agree);
  }
}

#endif /* of commented out functions */

