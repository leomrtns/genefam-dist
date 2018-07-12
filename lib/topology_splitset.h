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

/*! \file topology_splitset.h
 *  \brief Functions that use only the split bipartitions of topologies -- treating them as unrooted usually. 
 *
 *  Use a "splitset" structure that copies the bipartition information of all nodes (so that original trees are
 *  untouched) and then modifies this splitset -- which is in fact composed of two vectors of bipartitions like in the
 *  recombination program biomc2.
 */

#ifndef _biomcmc_topology_splitset_h_
#define _biomcmc_topology_splitset_h_

#include "topology_common.h"

typedef struct splitset_struct* splitset;

struct splitset_struct
{
  int size, spsize, spr, spr_extra, rf, hdist, hdist_reduced; /*! \brief spr, extra prunes for spr, rf distances and hdist=assignment cost */
  int n_g, n_s, n_agree, n_disagree;
  bipartition *g_split, *s_split, *agree, *disagree, *sp0; /* sp0 points to vec[0], s_split points to vec[x] */
  bipartition prune;
  hungarian h; /* hungarian method for solving the assignment between edges */
  bool match;  /*! \brief do we want to calculate the minimum cost assignment */
};

/*! \brief Allocate space for splitset structure (two vectors of bipartitions), for simple comparisons */
splitset new_splitset (int nleaves); 
/*! \brief Allocate space for splitset structure for dSPR calculation (also allocates aux vectors) */
splitset new_splitset_dSPR (int nleaves);
/*! \brief free memory allocated for splitset structure */
void del_splitset (splitset split);
/*! \brief Splitset structure for dSPR calculation for gene/species tree mapping */
splitset create_splitset_dSPR_genespecies (topology gene, topology species);
/*! \brief unrooted version of topologies comparison (on same leaf set), based on splitset previously allocated; recycle_t1 > 0 if we just compared it */
bool topology_is_equal_unrooted (topology t1, topology t2, splitset split, int recycle_t1);
/*! \brief approximate dSPR between unrooted topologies on same leaf set */
int dSPR_topology (topology t1, topology t2, splitset split);
/*! \brief RF distance (symmetric difference) between unrooted topologies on same leaf set */
int dSPR_topology_rf (topology t1, topology t2, splitset split);
/*! \brief h distance (optimal cost by the Hungarian method) between unrooted topologies on same leaf set */
int dSPR_topology_hdist (topology t1, topology t2, splitset split);
/*! \brief approximate dSPR between unrooted gene and species trees (leafset mapping) */
int dSPR_gene_species (topology gene, topology species, splitset split);
/*! \brief RF distance between unrooted gene and species trees (leafset mapping) */
int dSPR_gene_species_rf (topology gene, topology species, splitset split);
/*! \brief h distance (edge disagreement assignment cost) between unrooted gene and species trees (leafset mapping) */
int dSPR_gene_species_hdist (topology gene, topology species, splitset split);

#endif
