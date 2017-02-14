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

/*! \file topology_common.h 
 *  \brief General-purpose topology structures created from nexus_tree_struct (and low-level functions)
 *
 *  The topology structure should actually be called "tree" since it has information about branch lengths, but these
 *  functions neglect branch lenght information. Here we have functions that create split bipartitions for edges (stored
 *  by nodes below the edge) and compare distinct topologies based on these bipartitions. 
 *  We also have here the lowest-level function that apply an SPR on a topology (again, without caring about the branch
 *  length). 
 */

#ifndef _biomcmc_topology_common_h_
#define _biomcmc_topology_common_h_

#include "bipartition.h"
#include "hashtable.h"
#include "empirical_frequency.h"

typedef struct topol_node_struct* topol_node;
typedef struct topology_struct* topology;
typedef struct reconciliation_struct* reconciliation;

/*! \brief Information of a node (binary tree). */
struct topol_node_struct
{
  topol_node up, right, left, sister; /*! \brief Parent, children and sister nodes. */
  int id, level; /*! \brief Node ID (values smaller than nleaves indicate leaves) and distance from root. */
  int mid[5];    /*! \brief Mapping between nodes and postorder vectors [0,1] (postorder, undone); idx for deep coal [2,3] and losses [4] */
  bool internal, /*! \brief If internal node, TRUE; if leaf, FALSE. */ 
       u_done,   /*! \brief Has the topology up this edge (eq. to node) changed? (needed in likelihood calc) */
       d_done;   /*! \brief Has the topology down this edge (eq. to node) changed? (needed in likelihood calc) */
  bipartition split;    /*! \brief bipartition with information about leaves below node */ 
};

/*! \brief Binary unrooted topology (rooted at leaf with ID zero) */
struct topology_struct
{
  topol_node *nodelist; /*! \brief vector of nodes (the first \f$L\f$ are the leaves). */
  double *blength;      /*! \brief Branch lengths, with mean, min, max vectors for topology_space */
  int id;               /*! \brief topology ID (should be updated by hand, e.g. by functions in topology_space.c) */
  topol_node root;      /*! \brief Pointer to root node. */
  int nleaves;          /*! \brief Number of leaves \f$L\f$. */
  int nnodes;           /*! \brief Number of nodes, including leaves (\f$ 2L-1\f$ for a binary rooted tree). */
  topol_node undo_prune;   /*! \brief How to revert most recent SPR move (prune node). */
  topol_node undo_regraft; /*! \brief How to revert most recent SPR move (regraft node). */
  bool undo_lca;           /*! \brief revert SPR move is lca type or not */
  topol_node *postorder;   /*! \brief pointers to all internal nodes in postorder (from last to first is preorder) */
  topol_node *undone;      /*! \brief pointers to outdated nodes in postorder (from last to first is preorder) */
  int n_undone;  /*! \brief number of outdated nodes (which need likelihood calc etc) in topology_struct::undone. */
  unsigned int hashID1, hashID2; /*! \brief hash values of tree, ideally a unique value for each tree (collisions happen...) */
  bool traversal_updated;  /*! \brief zero if postorder[] vector needs update, one if we can use postdorder[] to traverse tree  */ 
  int ref_counter;         /*! \brief number of references of topology (how many places are pointing to it) */
  char_vector taxlabel;    /*! \brief Taxon names (just a pointer; actual values are setup by ::topology_space_struct or ::alignment_struct) */
  reconciliation rec;      /*! \brief in case of a gene tree, map nodes to foreign species tree */
  int *index;             /*! \brief sandbox vector used in spr moves / quasirandom tree shuffle just to avoid recurrent allocation */
  bool quasirandom;        /*! \brief tells if quasi-random structure was initialized (and topology_struct::idx is properly set) */
  /* Most recent common ancestor table for topology (assuming this tree is species tree) */ 
  topol_node **mrca; /*! \brief lower triangular matrix of topol_nodes (LCA between topol_node::id (i-1) and j) of size (nnodes-1)*/
};

/*! \brief mapping between gene tree nodes (this) and (external) species tree nodes */
struct reconciliation_struct
{
  topology current_sptree; /*! \brief pointer to species tree, to guarantee map is always consistent */
  topol_node *map_d; /*! \brief Mapping of all nodes from gene to species (the first gene::nnodes are fixed) */
  topol_node *map_u; /*! \brief Mapping of all nodes from gene to species, assuming gene tree is upside down (unrooted TESTING version) */
  int *sp_id,    /*! \brief mapping of gene (this tolopogy) leaf to ID of taxon in species tree */
      *sp_count, /*! \brief how many copies of each species are present in this gene (used by deepcoal) */
      sp_size,   /*! \brief effective number of species present in gene family */
      size_diff, /*! \brief twice the difference in number of leaves between gene tree and (reduced/effective) species tree */
      *dup,      /*! \brief indexes of duplication nodes on gene tree, and number of such nodes (unused for now) */
      *ndup_d,   /*! \brief number of duplications below node (edge above node, since struct assume node == edge above it) */
      *ndup_u,   /*! \brief number of duplications above node (edge upside down, thus "children" are 'up' and 'sister') */
      *nlos_d,   /*! \brief number of losses below node and edge above node */  
      *nlos_u,   /*! \brief number of losses above node, including edge above it */
      ndups,     /*! \brief minimum number of duplications over all possible rootings, acc. to reconciliation_struct::dup */
      nloss,     /*! \brief number of losses corresponding to rooting (edge) that minimizes duplications */
      ndcos;     /*! \brief total number of deep coalescences (from nloss - 2 X ndups + size_diff) */ 
};


/*! \brief Allocate space for new topology_struct */
topology new_topology (int nleaves);
/*! \brief Allocate vector for branch lengths (3 vectors: mean, min and max values observed in topol_space collection) */
void topology_malloc_blength (topology tree);
/*! \brief Free space allocated by topology_struct */
void del_topology (topology topol);

/* DEBUG function */
void debug_topol (topology tree);

/*! \brief Allocate space for new reconciliation_struct (other functions defined in topology_mrca.c) */
reconciliation new_reconciliation (int gene_nleaves, int sp_nleaves);
/*! \brief Create new reconciliation struct and copy values (except species tree info) from another struct */
reconciliation new_reconciliation_from_reconciliation (int gene_nleaves, int sp_nleaves, reconciliation from);
/*! \brief release allocated memory for reconciliation_struct */
void del_reconciliation (reconciliation r);

/*! \brief Copy information from topology_struct. 
 *
 * Since IDs do not change, this function only needs to update topol_node_struct::up, topol_node_struct::right, 
 * and topol_node_struct::left pointers and topol_node_struct::map_id from internal nodes; update of 
 * topol_node::sister is handled by function update_topology_sisters(). 
 * \param[in]  from_tree original topology_struct 
 * \param[out] to_tree (previously allocated) copied topology_struct
 */
void copy_topology_from_topology (topology to_tree, topology from_tree);

/*! \brief Update pointers to topol_node_struct::sister */
void update_topology_sisters (topology tree);

/*! \brief Update topol_node::preorder, topol_node::postorder, topol_node::bipartition 
 * and order siblings by number of descendants. */ 
void update_topology_traversal (topology tree);

/*! \brief Compare two topologies based on bipartitions (not on branch lengths) */
bool topology_is_equal (topology t1, topology t2);

/*! \brief Boolean if node2 is on the path of node1 to the root. */
bool node1_is_child_of_node2 (topol_node node1, topol_node node2);

/*! \brief allocate memory for a new distance_matrix that will be used on topologies */
distance_matrix new_distance_matrix_for_topology (int nleaves);

/*! \brief fill in distance_matrix with the patristic distances from topology (can be used with distinct branch length vectors to fill upper and lower diagonals */
void fill_distance_matrix_from_topology (distance_matrix dist, topology tree, double *blen, bool use_upper);

/*! \brief Print subtree in newick format to string using leaf IDs.
 *
 * Stores in string the tree in newick format, using leaf ID numbers (in practical applications needs a TRANSLATION 
 * nexus block). Memory allocation is handled by this function, but needs to be freed by the calling function. 
 * \param[in] tree tree to be printed
 * \param[in] blen vector with branch lengths (usually tree->blength)
 * \return a pointer to newly allocated string 
 */
char * topology_to_string_by_id (const topology tree, double *blen);

/*! \brief Print subtree in newick format to string creating names (based on leaf IDs.)
 *
 * Stores in string the tree in newick format, using newly-created names based on leaf ID numbers (useful for generating
 * random trees that must be read by other programs.) Memory allocation is handled by this function, but needs to be freed by the calling function. 
 * \param[in] tree tree to be printed
 * \param[in] blen vector with branch lengths (usually tree->blength)
 * \return a pointer to newly allocated string 
 */
char * topology_to_string_create_name (const topology tree, double *blen);

/*! \brief Print subtree in newick format to string using leaf names.
 *
 * Stores in string the tree in newick format, preserving sequence names if available.
 * Memory allocation is handled by this function, but needs to be freed by the calling 
 * function. 
 * \param[in] tree tree to be printed
 * \param[in] blen vector with branch lengths (usually tree->blength)
 * \return a pointer to newly allocated string */
char * topology_to_string_by_name (const topology tree, double *blen);

/*! \brief Prints subtree in dot format to file.
 *
 * Prints to file the tree in dot format (undirected graph). The dot format can be used with the 
 * <a href="http://www.graphviz.org/">graphviz</a>  suite of programs, and is not restricted to trees but can also 
 * handle arbitrary graph structures. Notice that we do not make use of the graphviz library, we simply create the 
 * text file graphviz programs take as input. Unfortunately, it is not helpful to print the nexus_tree_struct since the
 * program works basically with the topology_struct. On the other hand it is easy to change this function to make it 
 * work with topology_struct.
 * \param[in,out] fout pointer to file handler where tree is to be printed;
 * \param[in] label graph name or any other label;
 * \param[in] tree topology_struct to be printed; */
void graphviz_file_topology (FILE * fout, char *label, const topology tree);

/*! \brief Apply one subtree prune-and-regraft (SPR branch swapping) operation at specified nodes. 
 *
 * Each node is associated to one edge (the branch immediately above it), thus the location of the regraft node will 
 * impose the direction of pruning - the prune edge will always detach away from subtree containing regraft.
 * The actual SPR move needs to handle two cases: <b>prune node is in the path from regraft node to the root</b> 
 * (prune node is least common ancestor between prune and regraft) and <b>prune node is not in the path from 
 * regraft node to root</b> (prune and regraft nodes share a distinct common ancestor). When prune node is the root, 
 * the first case implies in rerooting. Checking against illegal moves (prune==regraft, prune==regraft->up, etc) should
 * be done previous to this function call. This function will call the corresponding lower-level one based on position
 * of prune node. If you know the direction of pruning (rerooting, e.g.) you can call the other two functions directly. 
 *
 *  \param[in,out] p topology over which to apply move
 *  \param[in] prune node to be pruned (detached). Direction determined by regraft
 *  \param[in] regraft node above which prune node will be reattached*/
void apply_spr_at_nodes (topology p, topol_node prune, topol_node regraft, bool update_done);
/*! \brief Apply one SPR branch swapping at specified nodes when prune subtree is above prune node.*/ 
void apply_spr_at_nodes_LCAprune (topology tree, topol_node prune, topol_node regraft, bool update_done);
/*! \brief Apply one SPR branch swapping at specified nodes when subtree to be pruned is below prune node.*/ 
void apply_spr_at_nodes_notLCAprune (topology tree, topol_node prune, topol_node regraft, bool update_done);
/*! \brief revert last SPR branch swapping */
void topology_undo_random_move (topology tree, bool update_done);
/*! \brief reset all d_done and u_done booleans to "true" (when rejecting a new state in MCMC) */
void clear_topology_flags (topology tree);
/*! \brief reset all d_done and u_done booleans to "false" (when updating a model parameter with MTM) */
void raise_topology_flags (topology tree);
/*! \brief revert last SPR branch swapping and clear flags (reject last proposal, in MCMC)  */
void topology_reset_random_move (topology tree);

/*! \brief store ID of each node's parent (in postorder) into vector, returning number of stored nodes */
int copy_topology_to_intvector_by_postorder (topology tree, int *ivec);
/*! \brief restore topological structure based on postordered ID vector, returning number of restored nodes */
int copy_intvector_to_topology_by_postorder (topology tree, int *ivec);
/*! \brief store ID of each node's parent into vector */
void copy_topology_to_intvector_by_id (topology tree, int *ivec);
/*! \brief restore topological structure based on ID vector */
void copy_intvector_to_topology_by_id (topology tree, int *ivec);

#endif
