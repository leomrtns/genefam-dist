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

/*! \file genefam_dist.h 
 *  \brief biomcmc library interface to external programs, specific to genefam_dist repo.
 *
 *  The idea is for genefam_dist is to be a general library for the analysis of multi-labelled gene family trees, 
 *  including treesignal but also for species tree inference and distance calculation.
 *  This library started as a branching from the biomcmc library (from the guenomu software)
 */

#ifndef _genefam_dist_h_
#define _genefam_dist_h_

#include "topology_space.h"
#include "topology_mrca.h"
#include "topology_build.h"
#include "prob_distribution.h" 
#include "nexus_common.h" // opaque library called by alignment.c, but should also be visible to other progs 

#ifdef THESE_ARE_COMMENTS
#include "lowlevel.h" // called by bipartition.h, hashtable.h, random_number.h, topology_common.h etc
#include "hashtable.h"    // called by alignment.h topology_space.h common.h 
#include "bipartition.h"  // called by topology_common.h
#include "topology_common.h" // called by topology_mrca.h, topology_build.h and topology_splitset.h
#include "topology_splitset.h" // called by topology_space.h
#include "random_number_gen.h" // called by random_number.h
#include "random_number.h"   // called by topology_build.h
#include "empirical_frequency.h" // called by topology_mrca.h, nexus_common.h  and alignment.h 
#endif // of THESE_ARE_COMMENTS

#define NDISTS 6 

typedef struct genetree_struct* genetree;
typedef struct sptree_ratchet_struct* sptree_ratchet; // gene_sptrees in find_best_tree.c 
typedef struct genefam_sptree_struct* genefam_sptree; // accessible to other libraries

struct genetree_struct
{
  topology tree;
  double minmax[2 * NDISTS], dist[NDISTS]; // values are integers; might have a variable "scale" to allow for quick or no scaling 
  splitset split;
};

struct sptree_ratchet_struct  
{  // later this will have weights[NDISTS] to allow for 'random' score weighting assignment
  int n_genesamples, n_ratchet, n_proposal;
  int next_avail; // idx to ratchet elements 
  genetree *genesample;
  topology *ratchet, *proposal; 
  double *ratchet_score, best_score; 
};

struct genefam_sptree_struct
{
  int n_genefams, n_sptrees, next_sptree;
  genetree *genefam;
  topology *sptree;
  sptree_ratchet best_trees; // later this will be a vector since can be multi-threaded
};

/*** First implementation of treesignal relied on external species trees, and didn't have access to C structs 
 *   Legacy code below 
 ***/

/*! \brief given a gene tree and a group of species trees, both in newick format, return the spectrum of unnormalized distances */
int genefam_module_treesignal_fromtrees (const char *gtree_str, const char *splist_str, double **output_distances);
/*! \brief given a gene tree and a group of species trees, both in newick format, return the spectrum of rescaled distances */
int genefam_module_treesignal_fromtrees_rescale (const char *gtree_str, const char *splist_str, double **output_distances);
/*! \brief given a gene tree and a group of species trees, both in newick format, return frequencies of random sptrees with distances smaller than group */
int genefam_module_treesignal_fromtrees_pvalue (const char *gtree_str, const char *splist_str, int n_reps, double **output_distances);
/*! \brief given set of trees, return trees with SPR neighbours (as newick string). Useful for generating reference sptrees */
char* genefam_module_randomise_trees_with_spr (const char *splist_str, int n_copies, int n_spr);
/*! \brief generate a set of trees (as newick string), all separated from previous one by a number of SPR moves */
char* genefam_module_generate_spr_trees (int n_leaves, int n_iter, int n_spr);

#endif
