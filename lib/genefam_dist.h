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

/*! \file genefam_dist.h 
 *  \brief biomcmc library interface to external programs, specific to genefam_dist repo.
 *
 *  The idea is for biomcmc to be general for several sofware , including guenomu but also parallel projects by
 *  leomrtns. This particular version of the library is visible to the user as "genefam_dist", but the internal functions 
 *  treat it as "biomcmc" such that we can more easily export successful snippets to other programs/repositories.
 */

#ifndef _biomcmc_h_
#define _biomcmc_h_

#include "module_distance.h"
#include "topology_build.h"
#include "prob_distribution.h" 
#include "nexus_common.h" // opaque library called by alignment.c, but should also be visible to other progs 

#ifdef THESE_ARE_COMMENTS
#include "topology_space.h" // called by module_distance.h
#include "topology_splitset.h" // called by module_distance.h
#include "topology_mrca.h" // called by module_distance.h
#include "lowlevel.h" // called by bipartition.h, hashtable.h, random_number.h, topology_common.h etc
#include "hashtable.h"    // called by alignment.h topology_space.h common.h 
#include "bipartition.h"  // called by topology_common.h
#include "topology_common.h" // called by topology_space.h, topology_mrca.h, topology_build.h and topology_splitset.h
#include "random_number_gen.h" // called by random_number.h
#include "random_number.h"   // called by topology_build.h
#include "empirical_frequency.h" // called by topology_mrca.h and alignment.h
#endif // of THESE_ARE_COMMENTS

#endif
