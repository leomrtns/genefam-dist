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

/*! \file likelihood.h 
 *  \brief likelihood computation of a phylogeny (partial likelihood vectors) given the topology 
 *
 */

#ifndef _biomcmc_likelihood_h_
#define _biomcmc_likelihood_h_

#include "phylogeny.h"
#include "topology_common.h"

/*! \brief ln(likelihood) of topology, updating all internal nodes */ 
extern void (*ln_likelihood) (phylogeny phy, topology tre);
void accept_likelihood (phylogeny phy, topology tre);

/*! \brief ln(likelihood) of topology, based on changed nodes by dynamically updating lk_vector */
extern void (*ln_likelihood_moved_branches) (phylogeny phy, topology tre);
void accept_likelihood_moved_branches (phylogeny phy, topology tre);

/*! \brief ln(likelihood) of topology, based on changed nodes by statically updating lk_vector (under calling 
 * function control) */
extern void (*ln_likelihood_moved_branches_at_lk_vector) (phylogeny phy, topology tre, int idx);
void accept_likelihood_moved_branches_at_lk_vector (phylogeny phy, topology tre, int idx, double likelihood);

/*! \brief set likelihood functions to neglect alignment data, constant at one (ln = 0) [Bayesian prior] */
void set_likelihood_to_prior (void);
/*! \brief explicitly tell program that we must calculate likelihoods (simulating posterior distribution); set by default */
void set_likelihood_to_posterior (void);

#endif
