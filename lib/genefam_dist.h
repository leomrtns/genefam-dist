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
 *  The idea is for genefam_dist is to be general for several sofware, including treesignal but also parallel projects by
 *  leomrtns. This library started as a branching from the biomcmc library (from the guenomu software)
 */

#ifndef _genefam_dist_h_
#define _genefam_dist_h_

#include "topology_space.h"
#include "topology_splitset.h"
#include "topology_mrca.h"
#include "topology_build.h"
#include "prob_distribution.h" 
#include "nexus_common.h" // opaque library called by alignment.c, but should also be visible to other progs 

#ifdef THESE_ARE_COMMENTS
#include "lowlevel.h" // called by bipartition.h, hashtable.h, random_number.h, topology_common.h etc
#include "hashtable.h"    // called by alignment.h topology_space.h common.h 
#include "bipartition.h"  // called by topology_common.h
#include "topology_common.h" // called by topology_space.h, topology_mrca.h, topology_build.h and topology_splitset.h
#include "random_number_gen.h" // called by random_number.h
#include "random_number.h"   // called by topology_build.h
#include "empirical_frequency.h" // called by topology_mrca.h and alignment.h
#endif // of THESE_ARE_COMMENTS

/*! \brief given a gene tree and a group of species trees, both in newick format, return the spectrum of unnormalized distances */
int genefam_module_treesignal_fromtrees (const char *gtree_str, const char *splist_str, double **output_distances);
/*! \brief given a gene tree and a group of species trees, both in newick format, return frequencies of random sptrees with distances smaller than group */
int genefam_module_treesignal_fromtrees_pvalue (const char *gtree_str, const char *splist_str, int n_reps, double **output_distances);
/*! \brief generate a set of trees (as newick string), all separated from previous one by a number of SPR moves */
char* genefam_module_generate_spr_trees (int n_leaves, int n_iter, int n_spr);

#endif
