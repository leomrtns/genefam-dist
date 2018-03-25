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

/*! \file topology_build.h 
 *  \brief Creation and modification of topologies  
 *
 */

#ifndef _biomcmc_topology_build_h_
#define _biomcmc_topology_build_h_

#include "topology_common.h"
#include "random_number.h"

/*! \brief low level function that generates a random tree (equiv. to random refinement of a star topology) */
void randomize_topology (topology tree);
/*! \brief generates a random topology if sample_type==0, but can reuse some info later to create a "correlated" tree */
void quasi_randomize_topology (topology tree, int sample_type);
/*! \brief lowlevel UPGMA (or single-linkage) function that depends on a topology and a matrix_distance */
void upgma_from_distance_matrix (topology tree, distance_matrix dist, bool single_linkage);
/*! \brief lowlevel bioNJ function (Gascuel and Cuong implementation) that depends on a topology and a matrix_distance */
void bionj_from_distance_matrix (topology tree, distance_matrix dist) ;
/*! \brief updates distances between species based on genes and gene-to-species mapping, with min on upper and mean on lower diagonal  */
void fill_species_dists_from_gene_dists (distance_matrix spdist, distance_matrix gendist, int *sp_id, bool use_upper_gene);
/*! \brief update global (over loci) species distances besed on local (within locus) species distances */
void update_species_dists_from_spdist (distance_matrix global, distance_matrix local, int *spexist);

int prepare_spdistmatrix_from_gene_species_map (spdist_matrix spdist, int *sp_id, int n_sp_id);
void fill_spdistmatrix_from_gene_dists (spdist_matrix spdist, distance_matrix gendist, int *sp_id, bool use_upper_gene);
void update_spdistmatrix_from_spdistmatrix (spdist_matrix global, spdist_matrix local);

/*! \brief random rerooting */
void topology_apply_rerooting (topology tree, bool update_done);
/*! \brief recursive SPR over all internal nodes, assuming common prob of swap  per node */
void topology_apply_shortspr (topology tree, bool update_done);
/*! \brief recursive SPR over all internal nodes, using prob[] vector as rough guide of error rate for node */
void topology_apply_shortspr_weighted (topology tree, double *prob, bool update_done);
/*! \brief random Subtree Prune-and-Regraft branch swapping for subtree below lca node */
void topology_apply_spr_on_subtree (topology tree, topol_node lca, bool update_done);
/*! \brief random Subtree Prune-and-Regraft branch swapping */
void topology_apply_spr (topology tree, bool update_done);
/*! \brief random Subtree Prune-and-Regraft branch swapping generalized (neglecting root) */
void topology_apply_spr_unrooted (topology tree, bool update_done);
/*! \brief random Nearest Neighbor Interchange branch swapping (SPR where regraft node is close to prune node) */
void topology_apply_nni (topology tree, bool update_done);
/*! \brief check if it is possible to apply SPR/NNI without rerooting (used by topology_apply_spr() and MCMC functions) */
bool cant_apply_swap (topology tree);

#endif
