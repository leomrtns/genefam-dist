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

/*! \file topology_mrca.h
 *  \brief Associates a table of most recent common ancestors to all node pairs of a topology.
 *  Gene family and species tree functions (but the MRCA info is only on the species tree).
 *
 *  A gene family is a group of sequences with possibly more than one copy for same species and at the same time no
 *  representative for some other species. 
 *  Based on the mapping between leaves of the gene family tree (which have the topology::rec) and an external species 
 *  topology (represented by the topology::mrca) we map every internal node of the gene tree to an internal node of 
 *  the species tree (reconciliation). The topology::mrca vector (a triangular matrix) must be reset whenever the
 *  topology changes (it is properly called by update_tree_traversal()), and then updated only when needed (by
 *  gene_tree_reconcile()). An alternative would be to reset and then update the whole matrix at once, but might be a
 *  waste in the case where not all species taxa are represented in the gene tree. This alternative version however has
 *  worst time O(Nlog N) while our current approach has worst case O(n^2) (N for species tree and n for gene tree)...
 *
 *  If our only information about the species are the taxon names, we can call the function char_vector_longer_first_order()
 *  from file alignment.h to order them and then we can call index_ordered_sptaxa_to_genetaxa() directly. In some cases,
 *  however, we also have topologies associated with the species names (for example if the species trees are read from a
 *  topology_space_struct). And theoretically the order of species names in a topology_space_struct could be given by an
 *  external hashtable_struct (for instance, if we have the alignment <b>and</b> the respective topologies then the
 *  indexing should be the same between leaves and sequences).
 *  In such scenarios we cannot change the order of char_vector elements and thus we must use an empfreq_struct to hold
 *  the mapping between original order and decreasing size order. The first strategy is used by the guenomu MCMC
 *  sampler and by the guenomu posterior probability analyser; the higher level general functions like 
 *  init_tree_recon_from_species_mrca() use the second approach. */

#ifndef _biomcmc_topology_mrca_h_
#define _biomcmc_topology_mrca_h_

#include "topology_common.h"

/*! \brief Allocate space for mrca structure (a triangular matrix) */
void new_mrca_for_topology (topology t); 
/*! \brief Finds the MRCA between the node pair if it isn't stored already and returns a pointer to node */
topol_node mrca_between_nodes (topology topol, int i, int j);

/*! \brief initialize reconciliation_struct based on existing topology - pointed by mrca */
void init_tree_recon_from_species_topology (topology gene, topology species);
/*! \brief initialize reconciliation_struct when we only know the species taxon names */
void init_tree_recon_from_species_names (topology gene, char_vector sptaxlabel);

/*! \brief find occurences of species->string[] inside gene->string[] filling indexes in sp_idx_in_gene.
 *
 *  The species are taxon names which may be associated with topologies or alignments, such that we can not reorder its
 *  elements from longer to shorter (which is essential for pattern finding). So we must use a local or external vector
 *  with the mapping between original and sorted species names. If this function (and not the
 *  index_sptaxa_to_reconciliation() below) is used, the initialization of rec->sp_count[] must be done by hand (by
 *  calling initialize_reconciliation_sp_count()). */
void index_sptaxa_to_genetaxa (char_vector species, char_vector gene, int *sp_idx_in_gene, empfreq ef_external);
/*! \brief find occurences of ordered species->string[] inside gene->string[] filling indexes in rec->sp_id[] and updating
 * rec->sp_count[]. 
 *
 *  This function can only be used whe wen know that the species names are ordered from longer to shorter. This happens
 *  only when there are no topologies (where the leaves are the indexes of taxlabel) associated with the species'
 *  taxlabels. So they cannot be used within a topology_space_struct but can be used if our primary data is only these
 *  taxlabels (like in the guenomu MCMC algrithm). Calls initialize_reconciliation_sp_count() automatically. */
void index_sptaxa_to_reconciliation (char_vector species, char_vector gene, reconciliation rec);

/*! \brief Fill rec->sp_count[] with the number of representatives of each species (idexed by rec->sp_id[]) */
void initialize_reconciliation_sp_count (reconciliation rec, int n_sp, int n_idx);

/*! \brief transform indexes found in index_sptaxa_to_genetaxa() to pointers to species nodes */
void initialize_reconciliation_from_species_tree (topology gene, topology species);

/*! \brief Find reconciliation map between gene and species trees */
void gene_tree_reconcile (topology gene, topology species);

#endif
