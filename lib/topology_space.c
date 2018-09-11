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

// OBS: topology_space object should always have <distinct> set, and <freq> scaled to one 

#include "topology_space.h"
#include "nexus_common.h"

#define DEFAULTBLENGTH 1. /*!< \brief Default branch length. */

typedef struct nexus_node_struct* nexus_node;
typedef struct nexus_tree_struct* nexus_tree;

/*! \brief Node information for each tree in tree file. */
struct nexus_node_struct
{
  nexus_node up, right, left; /*! \brief Parent and children nodes. */
  int id;               /*! \brief Initial pre-order numbering of node. */
  double branch_length; /*! \brief Branch length from node to node->up. */
  char *taxlabel;       /*! \brief Leaf sequence name (pointer to actual string in nexus_treespace_struct). */
};

/*! \brief Data from each tree in tree file. */
struct nexus_tree_struct
{
  nexus_node *nodelist;  /*! \brief Vector with pointers to every internal node. */
  nexus_node *leaflist;  /*! \brief Vector with pointers to tree leaves. */
  nexus_node root;       /*! \brief Pointer to root node. */
  bool has_branches;     /*! \brief Boolean saying if tree has branch lengths or not. */
  int nleaves;           /*! \brief Number of leaves (number of seqs in nexus_alignment_struct). */
  /*! \brief Number of nodes, including leaves. Since the tree is binary and rooted, the number of nodes equals 
   * \f$ 2L-1\f$ where \f$L\f$ is the number of leaves. */
  int nnodes;
};


/*! \brief Allocates memory for nexus_tree_struct. */
nexus_tree new_nexus_tree (int nleaves);

/*! \brief Allocates memory for topology_space_struct (set of trees present in nexus file).  */
topology_space new_topology_space (void);

/*! \brief Frees memory used by tree. */
void del_nexus_tree (nexus_tree T);

/*! \brief Reads tree from file and store in treespace. 
 *
 * Comment to myself about the external hashtable:
 * Given a hashtable with names, we want to number tree taxlabels according to these hash values. 
 * After checking if all names from taxlabel have a corresponding hash key, this function will create a vector with 
 * the position, in hash, of elements of taxlabel. This is what I call a mapping (!!). For instance, if we have the 
 * hashtable and taxlabels, the mapping is given as in the following:
 * \verbatim
 *     hash["A"] = 0       taxlabel[0] = "C"             order[0] = 2 
 *     hash["B"] = 1       taxlabel[1] = "B"   lead to   order[1] = 1
 *     hash["C"] = 2  and  taxlabel[2] = "D"   mapping   order[2] = 3
 *     hash["D"] = 3       taxlabel[3] = "E"             order[3] = 4
 *     hash["E"] = 4       taxlabel[4] = "A"             order[4] = 0
 * \endverbatim
 * Using this ordering, all trees belonging to a topology_space_struct will be relabeled by the mapping. 
 * Notice than even in the same topology_space_struct distinct trees may have distinct leaf label orders, 
 * despite the taxlabel[] vector with leaf names is shared.
 * This is necessary since despite nexus_tree_struct has freedom about the order of nodes (including leaves), in 
 * topology_struct the order is defined. */
void add_tree_to_topology_space (topology_space tsp, const char *string, bool translate, hashtable external_hash, int **order, bool use_root_location);

/*! \brief Reads translation table (one line) of the form "number = taxa name" in tree file. */
void translate_taxa_topology_space (topology_space tsp, char *string, hashtable external_hash);

/*! \brief Copy information from nexus_tree struct to topology struct 
 *
 * Since topology nodes are related to likelihood vectors their IDs follow strict rules: 
 * - their IDs should not change, only the relations between them represented by pointers up, left, right, sister;
 * - their IDs are the position in the topology_struct::nodelist vector (with actual nodes);
 * - IDs smaller than number of leaves represent the leaves (indexed by same IDs);
 * \param[in] nxs_tree nexus_tree_struct with node IDs respecting topology rules
 * \param[out] tree (previously allocated) copied topology_struct
 */
void copy_topology_from_nexus_tree (topology tree, nexus_tree nxs_tree);

/*! \brief Creates nexus_tree structure. 
 *  
 * Given a string with a tree in newick format, initialize the tree.
 * \param[in] string string with tree information in newick format (with or without branch lengths)
 * \return newly created nexus_tree_struct */
nexus_tree new_nexus_tree_from_string (char **string);

/*! \brief Recursive function that creates a node based on parenthetic structure. */
nexus_node subtree_nexus_tree (nexus_tree tree, char *lsptr, char *rsptr, int *node_id, nexus_node up);

/*! \brief Reads leaf name (or number, if translation table is present). */
char* read_taxlabel ( const char *name_start, const char *name_end);

/*! \brief Preorder initialization of leaves. */
void create_leaflist_nexus_tree (nexus_tree tree, nexus_node this, int *id);

/*! \brief Counts the number of leaves and resolves (one) trifurcation of tree string. */
int number_of_leaves_in_newick (char **string, bool resolve_trifurcation);

/*! \brief Preorder initialization of _internal_ nodes; 'id' should be >= nleaves. */
void create_node_id_nexus_tree (nexus_node this, int *id);

/*! \brief Searches for (last reported) branch length on string or return default value.  */
double read_branch_length (char *right_string_ptr);

/*! \brief Returns position of innermost comma (divides string into two subtrees). */
int find_branch_split_newick (char *left_string_ptr, char *right_string_ptr);

/*! \brief string with original file name, with extension stripped -- be caredul not to overwrite it on program */
void store_filename_in_topology_space (topology_space tre, char *filename);

/*! \brief Auxiliary function for the python module */
void
add_string_with_size_to_topology_space (topology_space *tsp, char *long_string, size_t string_size, bool use_root_location)
{
  char *local_string;
  int i, index;
  nexus_tree tree;
  topology topol;
  /* read string into a nexus tree */
  local_string = (char*) biomcmc_malloc (sizeof (char) * (string_size + 1));
  strncpy (local_string, long_string, string_size + 1); /* adds '\0' only when long_string is smaller!! */
  local_string[string_size] = '\0'; /* not null-terminated by default */
  tree  = new_nexus_tree_from_string (&local_string);
  if (local_string) free (local_string);
  topol = new_topology (tree->nleaves);
  /* create taxlabels (shared across trees) or prepare nexus tree leaves to be reordered */
  if (!(*tsp)) {
    (*tsp) = new_topology_space();
    /* tsp->taxlabel will point to names of first tree. This info will be also available at the hashtable */ 
    (*tsp)->taxlabel = new_char_vector (tree->nleaves);
    (*tsp)->taxlabel_hash = new_hashtable (tree->nleaves);
    for (i=0; i< tree->nleaves; i++) { 
      char_vector_link_string_at_position ((*tsp)->taxlabel, tree->leaflist[i]->taxlabel, i); 
      insert_hashtable ((*tsp)->taxlabel_hash, (*tsp)->taxlabel->string[i], i);
      tree->leaflist[i]->taxlabel = NULL; // we don't need this copy anymore
      tree->leaflist[i]->id = i;
    }
  }
  else {
    if ((*tsp)->taxlabel->nstrings != tree->nleaves) {
      del_nexus_tree (tree); del_topology_space (*tsp); del_topology (topol);
      biomcmc_error ( "All trees must have same number of leaves\n");
    }
    for (i=0; i< tree->nleaves; i++) { /* nexus tree leaflist IDs must follow hashtable */
      index = lookup_hashtable ((*tsp)->taxlabel_hash, tree->leaflist[i]->taxlabel);
      if (index < 0) {
        del_nexus_tree (tree); del_topology (topol); del_topology_space (*tsp);
        biomcmc_error ( "Leaf names are not the same across all trees\n");
      }
      if (tree->leaflist[i]->taxlabel) free (tree->leaflist[i]->taxlabel); 
      tree->leaflist[i]->id = index;
    }
  }
  /* now topology is ready to receive information from nexus tree */
  create_node_id_nexus_tree (tree->root, &i);
  if (tree->has_branches) topology_malloc_blength (topol); /* then it will copy length values from nexus_tree */
  copy_topology_from_nexus_tree (topol, tree);
  topol->taxlabel = (*tsp)->taxlabel; /* taxlabel is shared among all topologies */
  (*tsp)->taxlabel->ref_counter++;    /* since it is shared, it cannot be deleted if still in use */

  add_topology_to_topology_space_if_distinct (topol, (*tsp), use_root_location);

  del_nexus_tree (tree);
  return;
}

void
add_topology_to_topology_space_if_distinct (topology topol, topology_space tsp, bool use_root_location)
{ // Assumes freqs are counts, not normalized to one (currently this function is used by other functions that dont care for freqs)
 // Also, (1) doesn't update branch lengths; and (2) doesn't care about tsp->taxlabel
  int i, found_id = -1;
  tsp->tree = (topology*) biomcmc_realloc ((topology*) tsp->tree, sizeof (topology) * (tsp->ntrees+1));

  /* comparison includes root location (faster than unrooted since uses hash) */ 
  for (i=0; (i < tsp->ndistinct) && (found_id < 0); i++) if (topology_is_equal (topol, tsp->distinct[i])) found_id = i;
  if ((!use_root_location) && (tsp->ndistinct) && (found_id < 0)) { /* if they look distinct (different root), then do slower unrooted calculation */ 
    splitset split = new_splitset (topol->nleaves);
    for (i=0; (i < tsp->ndistinct) && (found_id < 0); i++) if (topology_is_equal_unrooted (topol, tsp->distinct[i], split, i)) found_id = i;
    del_splitset (split);
  }

  if (found_id >= 0) {
    tsp->tree[tsp->ntrees] = tsp->distinct[found_id];
    tsp->freq[found_id] += 1.;
    del_topology (topol);
  } // if tree was found 
  else { // then topol is unique 
    topol->id = tsp->ndistinct++;
    tsp->freq =     (double*)   biomcmc_realloc ((double*)   tsp->freq,     sizeof (double) * (tsp->ndistinct));
    tsp->distinct = (topology*) biomcmc_realloc ((topology*) tsp->distinct, sizeof (topology) * (tsp->ndistinct));
    tsp->freq[topol->id] = 1.;
    tsp->distinct[topol->id] = topol;
    if (topol->id > 0) for (i=0; i < tsp->distinct[topol->id]->nleaves; i++) {
      /* the leaf bipartitions never change, so can be shared among all topologies */
      del_bipartition (tsp->distinct[topol->id]->nodelist[i]->split);
      tsp->distinct[topol->id]->nodelist[i]->split = tsp->distinct[0]->nodelist[i]->split;
      tsp->distinct[0]->nodelist[i]->split->ref_counter++;
    }
  }
  tsp->ntrees++;
}

topology_space
read_topology_space_from_file (char *seqfilename, hashtable external_taxhash, bool use_root_location)
{
  return read_topology_space_from_file_with_burnin_thin (seqfilename, external_taxhash, 0, 1, use_root_location);
}

topology_space
read_topology_space_from_file_with_burnin_thin (char *seqfilename, hashtable external_taxhash, int burnin, int thin, bool use_root_location)
{
  topology_space treespace=NULL;
  FILE *seqfile;
  char *line=NULL, *line_read=NULL, *needle_tip=NULL;
  bool option_begin_trees    = false,
       option_translate_perm = false,
       option_translate_temp = false,
       option_include_tree = false;
  size_t linelength = 0;
  double freq_sum = 0., this_freq, *freq = NULL; /* posterior frequency per tree (if present) */
  int iteration = 1, i = 0, n_freq = 0, *order_external = NULL; // leaves will follow external_taxhash if exists (malloc'ed by add_tree_to_topology_space())
  if (burnin < 0) burnin = 0;
  if (thin < 1) thin = 1;

  seqfile = biomcmc_fopen (seqfilename, "r");
  
  /* the variable *line_read should point always to the same value (no line++ or alike) */
  biomcmc_getline (&line_read, &linelength, seqfile);
  line = remove_nexus_comments (&line_read, &linelength, seqfile);
  while (!nonempty_string (line)) {
    /* skip (possibly not following NEXUS format) initial comments and blank lines */
    if (biomcmc_getline (&line_read, &linelength, seqfile) < 0) 
      biomcmc_error ("Premature end of NEXUS tree file %s\n", seqfilename);
    line = remove_nexus_comments (&line_read, &linelength, seqfile);
  }
  if (!strcasestr (line, "NEXUS")) 
    biomcmc_error ( "%s is not ot a Nexus tree file (first line should be \"#NEXUS\")\n", seqfilename);

  while (biomcmc_getline (&line_read, &linelength, seqfile) != -1) {
    /* read frequency ('posterior distribution) information, in mrbayes' .trprobs format -> before remove_comments */
    needle_tip = line_read;
    if ((iteration > burnin) && !(iteration%thin)) option_include_tree = true;
    else option_include_tree = false;

    if ( option_include_tree && (needle_tip = strcasestr (needle_tip, "TREE")) && 
         (needle_tip = strrchr (needle_tip, '=')) &&
         (sscanf (needle_tip, "= [ &W %lf ]", &this_freq) == 1) ) {
      freq = (double*) biomcmc_realloc ((double*) freq, (n_freq + 1) * sizeof (double));
      freq[n_freq++] = this_freq;
      freq_sum += this_freq;
    }

    line = remove_nexus_comments (&line_read, &linelength, seqfile);
    if (nonempty_string (line)) { 
      /* don't do anything untill 'BEGIN TREES' block */
      if ((!option_begin_trees) && (strcasestr (line, "BEGIN TREES"))) {
        option_begin_trees = true;
        treespace = new_topology_space ();
      } 

      else if (!option_translate_temp) {
        /* check if we need to translate taxon names; in any case see if we have trees to read */
        if (strcasestr (line, "TRANSLATE")) {
          option_translate_perm = true;
          option_translate_temp = true;
        }
        else if (strcasestr (line, "TREE") && (needle_tip = strcasestr (line, "="))) {
          needle_tip++; /* remove "=" from string */
          iteration++;
          if (option_include_tree)
            add_tree_to_topology_space (treespace, needle_tip, option_translate_perm, external_taxhash, &order_external, use_root_location);
        }
      }
    
      if (option_translate_temp) {
        /* we are reading translation table token <-> taxlabel */  
        translate_taxa_topology_space (treespace, line, external_taxhash);
        if (strchr (line, ';')) option_translate_temp = false;
      }
      
    } // if (line)
  } //while (biomcmc_getline)

  if (external_taxhash) {
    treespace->taxlabel_hash = external_taxhash;
    external_taxhash->ref_counter++; /* since we are sharing the hashfunction */

    // reorder taxlabels to conform to hashtable 
    char_vector_reorder_strings (treespace->taxlabel, order_external);
  }

  fclose (seqfile);
  if (order_external) free (order_external);
  if (line_read) free (line_read);

  /* char_vector_remove_empty_strings() also updates string lengths into taxlabel->nchars so we call it anyway */
  if (char_vector_remove_empty_strings (treespace->taxlabel)) /* this problem should have been detected before... */
    biomcmc_error ("empty taxon names in nexus tree file (reading problem or wrong/duplicate numbers in translate)");

  /* vector of frequencies */
  treespace->freq = (double*) biomcmc_malloc (treespace->ndistinct * sizeof (double));
  for (i=0; i < treespace->ndistinct; i++) treespace->freq[i] = 0.;

  /* temporarily use treespace->freq[] to calculate mean branch lengths and mean tree length for each topology */
  if (treespace->has_branch_lengths) {
    int j;
    /* count how many times tree was seen (neglect frequency over distinct trees; we are interested in mean lenghts for this one tree */
    for (i=0; i < treespace->ntrees; i++) treespace->freq[ treespace->tree[i]->id ] += 1.;
    /* each branch value is sum over all trees with same topol; must be scaled to mean value */
    for (i=0; i < treespace->ndistinct; i++) for (j=0; j < treespace->distinct[i]->nnodes; j++)
      treespace->distinct[i]->blength[j] /= treespace->freq[i];
    /* branch values are scaled to one; to recover original values we must multiply by tlen. Now we calculate avge tree length */
    for (i=0; i < treespace->ndistinct; i++) treespace->tlen[3*i] /= treespace->freq[i]; /* [mean, min, max] sum of branch lengths */
    for (i=0; i < treespace->ndistinct; i++) treespace->freq[i] = 0.; /* clean up temp values */
  }

  /* Now we store the proper frequency values (over all trees) in treespace->freq[], for weighted or unweighted files */
  if ((freq) && (n_freq == treespace->ntrees)) { /* we have posterior frequency information */
    for (i=0; i < treespace->ntrees; i++)    treespace->freq[ treespace->tree[i]->id ] += freq[i];
    for (i=0; i < treespace->ndistinct; i++) treespace->freq[i] /= freq_sum; /* normalize to one */
  }
  else {
    for (i=0; i < treespace->ntrees; i++)    treespace->freq[ treespace->tree[i]->id ] += 1.;
    for (i=0; i < treespace->ndistinct; i++) treespace->freq[i] /= (double) treespace->ntrees; /* BUG found in 2013.02.28 (was dividing by ndistinct) */
  }
  if (freq) free (freq);
  
  store_filename_in_topology_space (treespace, seqfilename);
  return treespace;
}

void
merge_topology_spaces (topology_space ts1, topology_space ts2, double weight_ts1, bool use_root_location)
{ /* ts1->tree is not correct anymore, should not be used after calling this function */
  int i, j, *idx, n_idx = ts1->ndistinct; 
  double total_freq = 0.;

  // TODO: check if taxlabel_hash is the same; branch lengths; mark where convergence could go.
  
  if (weight_ts1 <= 0.) weight_ts1 = 1.; // usually ts1->ntrees/ts2->ntrees

  idx = (int*) biomcmc_malloc (ts1->ndistinct * sizeof (int)); 
  for (i=0; i < ts1->ndistinct; i++) {
    ts1->freq[i] *= weight_ts1; // weight (number of trees) of ts1 relative to ts2
    idx[i] = i; /* index of trees from ts1 not compared to ts2 yet */
  }

  for (j=0; j < ts2->ndistinct; j++) {
    int found_id = -1;
    /* comparison includes root location (faster than unrooted since uses hash) */ 
    for (i=0; (i < n_idx) && (found_id < 0); i++) if (topology_is_equal (ts2->distinct[j], ts1->distinct[idx[i]] )) found_id = i;
    if ((!use_root_location) && (found_id < 0)) { /* if they look distinct (different root), then do slower unrooted calculation */ 
      splitset split = new_splitset (ts2->distinct[j]->nleaves);
      for (i=0; (i < n_idx) && (found_id < 0); i++) if (topology_is_equal_unrooted (ts2->distinct[j], ts1->distinct[idx[i]], split, i)) found_id = i;
      del_splitset (split);
    }
    if (found_id >= 0) {
      ts1->freq[ idx[found_id] ] += ts2->freq[j];
      idx [found_id] = idx[--n_idx]; // assuming each tree from ts1 can be found at most once in ts2 
    } // if tree not found 
    else { // tree ts2->distinct[j] is unique to ts2 
      int new_id = ts1->ndistinct++;
      ts2->distinct[j]->id = new_id;
      ts1->freq =     (double*)   biomcmc_realloc ((double*)   ts1->freq,     sizeof (double) * (ts1->ndistinct));
      ts1->distinct = (topology*) biomcmc_realloc ((topology*) ts1->distinct, sizeof (topology) * (ts1->ndistinct));
      ts1->freq[new_id] = ts2->freq[j];
      ts1->distinct[new_id] = ts2->distinct[j];
      ts2->distinct[j] = NULL;
      for (i=0; i < ts1->distinct[new_id]->nleaves; i++) {
        /* the leaf bipartitions never change, so can be shared among all topologies */
        ts1->distinct[new_id]->nodelist[i]->split = ts1->distinct[0]->nodelist[i]->split;
        ts1->distinct[0]->nodelist[i]->split->ref_counter++;
      }
    }
  } // for j in ts2->ndistinct

  for (i=0; i < ts1->ndistinct; i++) total_freq += ts1->freq[i];
  for (i=0; i < ts1->ndistinct; i++) ts1->freq[i] /= total_freq;
  if (idx) free (idx);
}

/*
void 
sort_topology_space_by_frequency(topology_space tsp, double *external_freqs) 
{ 
  double part_sum = 0., freq, *local_freqs = tsp->freq, *pivot_d;
  topology pivot_t;
  char *stree;
  empfreq_double efd;
  if (external_freqs) local_freqs = external_freqs; 
  efd = new_empfreq_double_sort_decreasing (local_freqs, tsp->ndistinct);
  // FIXME: stopped here: must change tlen[3 x ndistinct] and freq[] if external is NULL. tree[i] is ponter to distinct, wont change
} */

void
save_topology_space_to_trprobs_file (topology_space tsp, char *filename, double credible)
{ // works even if tsp->freq is unscaled (i.e. total counts and not frequencies summing to one)
  int i, idx;
  FILE *stream;
  double part_sum = 0., freq, *scaled_freqs;
  char *stree;
  empfreq_double efd;

  if (credible > 1.) credible = 1.;

  stream = biomcmc_fopen (filename, "w");
  fprintf (stream, "#NEXUS\n[While frequency 'p' is unscaled, 'P' and 'W' are scaled by credible=%.4lf]\n", credible);
  fprintf (stream, "\n\nBegin trees;\n Translate\n");
  fprintf (stream, "\t1  %s", tsp->taxlabel->string[0]);
  for (i=1; i < tsp->distinct[0]->nleaves; i++) fprintf (stream, ",\n\t%d  %s", i+1, tsp->taxlabel->string[i]);
  fprintf (stream, "\n;\n");

  scaled_freqs = (double*) biomcmc_malloc (sizeof (double) * tsp->ndistinct);
  for (i = 0; i < tsp->ndistinct; i++) part_sum += tsp->freq[i];
  for (i = 0; i < tsp->ndistinct; i++) scaled_freqs[i] = tsp->freq[i]/part_sum;

  efd = new_empfreq_double_sort_decreasing (scaled_freqs, tsp->ndistinct);

  part_sum = 0.;
  for (i = 0; (i < tsp->ndistinct) && (part_sum < 1.); i++) {
    idx = efd->d[i].idx; /* element ordered by frequency */
    freq = scaled_freqs[idx] / credible; /* rescaling s.t. new frequencies sum to one (and not to "credible") */
    part_sum += freq;
    stree = topology_to_string_by_id (tsp->distinct[idx], false); /* from topol to newick */
    fprintf (stream, "tree tree_%d \t[p= %.5lf, P= %.5lf] = [&W %.8lf] %s;\n", i, tsp->freq[idx], part_sum, freq, stree);
    free (stree);
  }
  fprintf (stream, "\nEnd;\n");

  fclose (stream);
  del_empfreq_double (efd);
  if (scaled_freqs) free (scaled_freqs);
}

int
estimate_treesize_from_file (char *seqfilename)
{
  FILE *seqfile;
  char *line=NULL, *line_read=NULL, *needle_tip=NULL;
  size_t linelength = 0;
  int this_size, size = 0, ntrees = 0; 

  seqfile = biomcmc_fopen (seqfilename, "r");
  /* the variable *line should point always to the same value (no line++ or alike) */
  while ((biomcmc_getline (&line_read, &linelength, seqfile) != -1) && (ntrees < 10)) {
    line = remove_nexus_comments (&line_read, &linelength, seqfile);
    if (strcasestr (line, "TREE") && (needle_tip = strcasestr (line, "="))) {
      needle_tip++; /* remove "=" from string */
      this_size  = number_of_leaves_in_newick (&needle_tip, false); /* false = do not attempt to change tree string */
      if (this_size) { size += this_size; ntrees++; }
    }
  }
  fclose (seqfile);
  if (line_read) free (line_read);

  if (!ntrees) return -1;
  return size/ntrees;
}

topology_space
new_topology_space (void) 
{
  topology_space tsp;
  
  tsp = (topology_space) biomcmc_malloc (sizeof (struct topology_space_struct));

  tsp->ndistinct = tsp->ntrees = 0;
  tsp->tree = tsp->distinct = NULL;
  tsp->taxlabel = NULL;
  tsp->taxlabel_hash = NULL;
  tsp->filename = NULL;
  tsp->has_branch_lengths = false;
  /* tsp->dinstinct vector is increased by add_tree_to_topology_space() while tsp->tree are pointers to dsp->distinct
   * tsp->taxlabel vector is setup by translate_taxa_topology_space() or add_tree_to_topology_space(), whichever 
   * is used
   * tsp->taxlabel_hash is setup by one of the above in absence of global external_hash -- will be replaced by a pointer to an external hastable (from
   * alignment or another tree file, for instance) */
  tsp->freq = NULL; /* will be created online (as we read more trees) or when reading "[&W 0.01]" posterior probability data */
  tsp->tlen = NULL; /* created only if trees have branch lengths */
  return tsp;
}

void
del_topology_space (topology_space tsp) 
{
  if (tsp) {
    int i;
    if (tsp->distinct) { 
      for (i=tsp->ndistinct-1; i>=0; i--) del_topology (tsp->distinct[i]);
      free (tsp->distinct);
    }
    if (tsp->tree)     free (tsp->tree);
    if (tsp->freq)     free (tsp->freq);
    if (tsp->tlen)     free (tsp->tlen);
    if (tsp->filename) free (tsp->filename);
    del_hashtable (tsp->taxlabel_hash);
    del_char_vector (tsp->taxlabel);
    free (tsp);
  }
}

nexus_tree
new_nexus_tree (int nleaves) 
{
  nexus_tree tree;
  int i;
  size_t sizeof_node = sizeof (struct nexus_node_struct);

  tree = (nexus_tree) biomcmc_malloc (sizeof (struct nexus_tree_struct));
  tree->nleaves = nleaves;
  tree->nnodes  = 2*nleaves - 1;
  tree->has_branches = false;
  
  tree->nodelist = (nexus_node*) biomcmc_malloc (tree->nnodes * sizeof (nexus_node));
  tree->leaflist = (nexus_node*) biomcmc_malloc (tree->nleaves * sizeof (nexus_node));
  
  /* tree->nodelist will store the actual nodes */
  for (i=0; i<tree->nnodes; i++) { 
    tree->nodelist[i] = (nexus_node) biomcmc_malloc (sizeof_node);
    tree->nodelist[i]->up = tree->nodelist[i]->right = tree->nodelist[i]->left = NULL;
    tree->nodelist[i]->taxlabel = NULL;
  }

  return tree;
}

void 
del_nexus_tree (nexus_tree T) 
{
  if (T) {
    int i;
    for (i=T->nnodes-1; i >=0; i--) if (T->nodelist[i]) free (T->nodelist[i]);
    if (T->nodelist) free (T->nodelist);
    if (T->leaflist) free (T->leaflist);
    free (T);
  }
}

void
add_tree_to_topology_space (topology_space tsp, const char *string, bool translate, hashtable external_hash, int **order, bool use_root_location)
{
  int i, index, *original_order;
  char *local_string;
  double treelength = 0.;
  splitset split;
  nexus_tree tree;
  topology topol;

  /* use local copy of string to avoid problems with biomcmc_getline() */
  local_string = (char*) biomcmc_malloc (sizeof (char) * (strlen (string) + 1));
  strcpy (local_string, string);
  tree  = new_nexus_tree_from_string (&local_string);
  if (local_string) free (local_string);

  topol = new_topology (tree->nleaves);

  if ((tsp->ntrees == 0) && (!translate)) { /* CASE 1: first tree read and no TRANSLATE command in nexus file */
    tsp->has_branch_lengths = tree->has_branches; /* first tree decides if branch lengths are taken into account */

    /* tsp->taxlabel will point to names of first tree. This info will be also available at the hashtable */ 
    tsp->taxlabel = new_char_vector (tree->nleaves); 

    for (i=0; i< tree->nleaves; i++) { 
      // leaf names of first nexus_tree will be shared by topol_space
      char_vector_link_string_at_position (tsp->taxlabel, tree->leaflist[i]->taxlabel, i); 
      tree->leaflist[i]->taxlabel = NULL; // we don't need this copy anymore
      tree->leaflist[i]->id = i;
    }

    if (external_hash) {
      /* leaves should be renumbered according to external taxlabel_hash */
      *order = (int*) biomcmc_malloc (2 * tsp->taxlabel->nstrings * sizeof (int)); /* two vecs, second is temporary */
      for (i=0; i < tsp->taxlabel->nstrings; i++) {
        /* map order in which taxlabels appear originally - where hashtable came from, e.g. the alignment file */
        (*order)[i] = lookup_hashtable (external_hash, tsp->taxlabel->string[i]);
        if ((*order)[i] < 0) {
          del_nexus_tree (tree); del_topology_space (tsp);
          biomcmc_error ( "tree label %s not found in sequence data\n", tsp->taxlabel->string[i]); 
        }
      }
    }
    else { /* no global (external) hash, must create local one */
      tsp->taxlabel_hash = new_hashtable (tree->nleaves);
      for (i=0; i< tree->nleaves; i++) insert_hashtable (tsp->taxlabel_hash, tsp->taxlabel->string[i], i);
    }
  }

  else if ((tsp->ntrees == 0) && (translate)) { /* CASE 2: first tree read and TRANSLATE command in nexus file */
    tsp->has_branch_lengths = tree->has_branches; /* first tree decides if branch lengths are taken into account */
    for (i=0; i< tree->nleaves; i++) {
      sscanf (tree->leaflist[i]->taxlabel, " %d ", &index);
      index--;
      if (index < 0 || index >= tree->nleaves) {
        del_nexus_tree (tree); del_topology_space (tsp);
        biomcmc_error ( "leaf number \'%d\' out of range (1...NTAX) in nexus tree TRANSLATE \n", index);
      }
      if (tree->leaflist[i]->taxlabel) free (tree->leaflist[i]->taxlabel); // we only care about the node ID 
      tree->leaflist[i]->id = index;
    }

    if (external_hash) { //same as before, remembering that taxlabel_hash doesn't exist 
      /* leaves should be renumbered according to external_hash */
      *order = (int*) biomcmc_malloc (2 * tsp->taxlabel->nstrings * sizeof (int));
      for (i=0; i < tsp->taxlabel->nstrings; i++) {
        /* map order in which taxlabels appear originally - where hashtable came from, e.g. the alignment file */
        (*order)[i] = lookup_hashtable (external_hash, tsp->taxlabel->string[i]);
        if ((*order)[i] < 0) {
          del_nexus_tree (tree); del_topology_space (tsp);
          biomcmc_error ( "tree label %s not found in external hash table with mapped names (from alignment, generally)\n", tsp->taxlabel[i]); 
        }
      }
    }
  }

  else if ((tsp->ntrees > 0) && (!translate)) { /* CASE 3: not the first tree read and no TRANSLATE command in nexus file */
    if (tsp->taxlabel->nstrings != tree->nleaves) {
      del_nexus_tree (tree); del_topology_space (tsp);
      biomcmc_error ( "number of leaves disagrees between trees of same file\n");
    }
    /* use hashtable to check if names are consistent and point all leaves to taxlabel vector */
    for (i=0; i< tree->nleaves; i++) {
      if (external_hash) index = lookup_hashtable (external_hash, tree->leaflist[i]->taxlabel);
      else          index = lookup_hashtable (tsp->taxlabel_hash, tree->leaflist[i]->taxlabel);
      if (index < 0) {
        del_nexus_tree (tree); del_topology_space (tsp);
        biomcmc_error ( "leaf names disagree between trees of same file\n");
      }
      if (tree->leaflist[i]->taxlabel) free (tree->leaflist[i]->taxlabel); 
      tree->leaflist[i]->id = index;
    }
  }

  else {             /*  CASE 4: not the first tree read and TRANSLATE command in nexus file */
    if (tsp->taxlabel->nstrings != tree->nleaves) {
      del_nexus_tree (tree); del_topology_space (tsp);
      biomcmc_error ( "number of leaves disagrees between tree and TRANSLATE command\n");
    }
    for (i=0; i< tree->nleaves; i++) {
      sscanf (tree->leaflist[i]->taxlabel, " %d ", &index);
      index--;
      if (index < 0 || index >= tree->nleaves) {
        int nleaveslocal = tree->nleaves;
        del_nexus_tree (tree); del_topology_space (tsp);
        biomcmc_error ( "leaf number \'%d\' out of range (1...NTAX) in nexus tree TRANSLATE -- NTAX = %d\n", index + 1, nleaveslocal);
      }
      if (tree->leaflist[i]->taxlabel) free (tree->leaflist[i]->taxlabel); 
      tree->leaflist[i]->id = index;
    }
  }

  create_node_id_nexus_tree (tree->root, &i);

  if (external_hash) {
    original_order = (*order) + tsp->taxlabel->nstrings;
    for (i=0; i < tree->nleaves; i++) original_order[i] = tree->leaflist[i]->id;
    for (i=0; i < tree->nleaves; i++) tree->leaflist[i]->id = (*order)[ original_order[i] ];
  }
  /* at this point nexus_tree is ready to be copied to topology */
  if (tsp->has_branch_lengths) topology_malloc_blength (topol); /* then it will copy length values from nexus_tree */
  copy_topology_from_nexus_tree (topol, tree);
  topol->taxlabel = tsp->taxlabel; /* taxlabel is shared among all topologies */
  tsp->taxlabel->ref_counter++;    /* since it is shared, it cannot be deleted if still in use */

  /* branch lenghts scaled to one, but original tree length stored in tsp->tlen */
  if (tsp->has_branch_lengths) {
    treelength= 0.;
    for (i = 0; i < topol->nnodes; i++)  treelength += topol->blength[i];
    for (i = 0; i < topol->nnodes; i++)  topol->blength[i] /= treelength;
    if (!tsp->ndistinct) tsp->tlen = (double*) biomcmc_malloc (3 * sizeof (double)); /* so that we can realloc() later */
  }
  /* look if same topology is already present */
  int found_id = -1;

  tsp->tree = (topology*) biomcmc_realloc ((topology*) tsp->tree, sizeof (topology) * (tsp->ntrees+1));

  split = new_splitset (topol->nleaves);
  /* distinct trees might have same unrooted info, but rooted calculation is faster */
  for (i=0; (i < tsp->ndistinct) && (found_id < 0); i++) if (topology_is_equal (topol, tsp->distinct[i])) found_id = i;
  if ((!use_root_location) && (found_id < 0)) /* if they look distinct (different root), then do slower unrooted calculation */ 
    for (i=0; (i < tsp->ndistinct) && (found_id < 0); i++) if (topology_is_equal_unrooted (topol, tsp->distinct[i], split, i)) found_id = i;
  del_splitset (split);
  // FIXME: better solution is to have vector of splits for each tree, since slowest part is to copy/order bipartitions 

  if (found_id >= 0) {
    tsp->tree[tsp->ntrees] = tsp->distinct[found_id];
    if (tsp->has_branch_lengths) {
      for (i = 0; i < topol->nnodes; i++) { /* first branches (which are already scaled to one) */
        tsp->distinct[found_id]->blength[i] += topol->blength[i]; /* mean length */
        if (tsp->distinct[found_id]->blength[i+  topol->nnodes] > topol->blength[i]) tsp->distinct[found_id]->blength[i+ topol->nnodes] = topol->blength[i]; /* min value */
        if (tsp->distinct[found_id]->blength[i+2*topol->nnodes] < topol->blength[i]) tsp->distinct[found_id]->blength[i+2*topol->nnodes] = topol->blength[i]; /*max value */
      }
      /* then total tree length (sum of unscaled branches) */
      tsp->tlen[3 * found_id] += treelength; /* has 3 consecutive elements: mean, min, max */
      if (tsp->tlen[3 * found_id + 1] > treelength) tsp->tlen[3 * found_id + 1] = treelength; /* min */
      if (tsp->tlen[3 * found_id + 2] < treelength) tsp->tlen[3 * found_id + 2] = treelength; /* max */
    }
    del_topology (topol);
  }
  else {
    topol->id = tsp->ndistinct++;
    if      (!(tsp->ndistinct%10000)) fprintf (stderr, "+");
    else if (!(tsp->ndistinct%1000))  fprintf (stderr, "."); /* the "else" is to avoid printing both */
    fflush (stdout); /* exponentially slower; maybe use hash values to find in linear time? */
    
    if (tsp->has_branch_lengths) {/* mean, min and max are same value, corresponding to the one in nexus_tree */
      for (i = 0; i < topol->nnodes; i++) topol->blength[i + topol->nnodes] = topol->blength[i + 2 * topol->nnodes] = topol->blength[i];
      tsp->tlen = (double*) biomcmc_realloc ((double*) tsp->tlen, 3 * sizeof (double) * (tsp->ndistinct));
      tsp->tlen[3 * topol->id] = tsp->tlen[3 * topol->id + 1] = tsp->tlen[3 * topol->id + 2] = treelength; 
    }

    tsp->distinct = (topology*) biomcmc_realloc ((topology*) tsp->distinct, sizeof (topology) * (tsp->ndistinct));
    tsp->distinct[tsp->ndistinct-1] = topol;
    tsp->tree[tsp->ntrees] = tsp->distinct[tsp->ndistinct-1];
    if (tsp->ndistinct > 1) for (i=0; i < tree->nleaves; i++) {
      /* the leaf bipartitions never change, so can be shared among all topologies */
      del_bipartition (tsp->distinct[tsp->ndistinct-1]->nodelist[i]->split);
      tsp->distinct[tsp->ndistinct-1]->nodelist[i]->split = tsp->distinct[0]->nodelist[i]->split;
      tsp->distinct[0]->nodelist[i]->split->ref_counter++;
    }
  }

  tsp->ntrees++;
  del_nexus_tree (tree);
}

/* TODO: create unroot_topol_space() */
void
translate_taxa_topology_space (topology_space tsp, char *string, hashtable external_hash) 
{
  int i, index;
  char *c, *s, *last, *comma_position, seqname[MAX_NAME_LENGTH]="";
  bool good_scan;

  /* the file may have the first token<->taxon_name at the same line as the "TRANSLATE" command. */
  if  ( (s = strcasestr (string, "TRANSLATE")) ) s += strlen ("TRANSLATE");
  else s = (string);
  last = s + strlen (s) +1;

  if (!tsp->taxlabel) tsp->taxlabel = new_char_vector (1); /* we don't know beforehand how many taxa */

  /* one or more "token1 taxon_id1, token2 taxon_id2," */
  while ((comma_position = strchr (s, ',')) && (s<last)) {
    if (strchr (s, '\"')) good_scan = (sscanf (s, " %d \"%[^\"]\",", &index, seqname) == 2); /* double quotes in name */
    else good_scan = (sscanf (s, " %d %[^,]", &index, seqname) == 2); /* leaf name is everything before a comma */ 
    if (!good_scan) biomcmc_error ("could not scan leaf info in TRANSLATE command");
    
    index--; /* in nexus we have index = 1...NTAX */
    if ((!strlen (seqname)) || (index<0)) biomcmc_error ( "unexpected leaf name/location in TRANSLATE command\n"); 
    char_vector_add_string_at_position (tsp->taxlabel, seqname, index);
    
    s = comma_position+1;
  }

  /* maybe only one "token taxon_id" (the last line, e.g.) */
  if (strchr (s, '\"')) good_scan = (sscanf (s, " %d \"%[^\"]\" ", &index, seqname) == 2); /* double quotes in name */
  else good_scan = (sscanf (s, " %d %[^,;] ", &index, seqname) == 2); /* leaf name is everything before a comma */
  if (good_scan) {
    good_scan = false; /* recicling the variable (use it to see if we found the end of translation) */
    while ( (c = strrchr (seqname, ';')) ) { 
      *c = '\0';  /* remove possible semicolon at end (even a freaky case of several ones!) */
      good_scan = true; /* A semicolon (even within seqname) means that translation table is completely read */
    }
    
    index--; /* in nexus we have index = 1...NTAX */
    if ((!strlen (seqname)) || (index<0)) biomcmc_error ( "unexpected leaf name/location in TRANSLATE command\n"); 
    char_vector_add_string_at_position (tsp->taxlabel, seqname, index);
  }
  
  /* when we finished reading the leaf names (semicolon separated (line below) or not (good_scan) by space) */ 
  if ((strchr (string, ';') || good_scan) && !external_hash) { /* if external_has is present, do not create local one */
    tsp->taxlabel_hash = new_hashtable (tsp->taxlabel->nstrings);
    for (i=0; i<tsp->taxlabel->nstrings; i++) insert_hashtable (tsp->taxlabel_hash, tsp->taxlabel->string[i], i);
  }
}

void
copy_topology_from_nexus_tree (topology tree, nexus_tree nxs_tree)
{
  int i, id, node_id;

  for (i = 0; i < tree->nnodes; i++) {
    node_id = nxs_tree->nodelist[i]->id;
    tree->nodelist[node_id]->mid[0] = tree->nodelist[node_id]->mid[1] = tree->nodelist[node_id]->id = node_id;
    if (tree->blength) tree->blength[node_id] = nxs_tree->nodelist[i]->branch_length; /* should be malloc'ed beforehand */

    if (nxs_tree->nodelist[i]->up) {
      id = nxs_tree->nodelist[i]->up->id;
      tree->nodelist[node_id]->up = tree->nodelist[id];
    }
    else tree->nodelist[node_id]->up = NULL;
    
    if (nxs_tree->nodelist[i]->left) {
      id = nxs_tree->nodelist[i]->left->id;
      tree->nodelist[node_id]->left = tree->nodelist[id];
    }
    else tree->nodelist[node_id]->left = NULL; 
    
    if (nxs_tree->nodelist[i]->right) {
      id = nxs_tree->nodelist[i]->right->id;
      tree->nodelist[node_id]->right = tree->nodelist[id];
    }
    else tree->nodelist[node_id]->right = NULL;
  } // for (nnodes)

  tree->root = tree->nodelist[nxs_tree->root->id];
  for (i = 0; i < tree->nleaves; i++) { tree->nodelist[i]->u_done = false; tree->nodelist[i]->d_done = true; }
  for (i = tree->nleaves; i < tree->nnodes; i++) tree->nodelist[i]->u_done = tree->nodelist[i]->d_done = false;

  update_topology_sisters (tree);
  update_topology_traversal (tree);
}

nexus_tree
new_nexus_tree_from_string (char **string) 
{
  char *lptr, *rptr;
  int id, nleaves;
  nexus_tree T; 
  
  *string = remove_space_from_string (*string);
  nleaves = number_of_leaves_in_newick (string, true); /* true = resolves trifurcation, changing string allocation */
  T = new_nexus_tree (nleaves);
  if (strchr (*string, ':')) T->has_branches = true;
  
  /* begin & end of string */
  lptr = *string;  
  rptr = *string + strlen (*string) - 1;
  
  id = 0; /* This function does the actual creation of the tree */
  T->root = subtree_nexus_tree (T, lptr, rptr, &id, NULL);
  
  id = 0; /* vector of pointers to the tree leaves */
  create_leaflist_nexus_tree (T, T->root, &id); 

  return T;
}

nexus_node
subtree_nexus_tree (nexus_tree tree, char *lsptr, char *rsptr, int *node_id, nexus_node up) 
{
  nexus_node thisnode;
  
  thisnode = tree->nodelist[*node_id];  
  thisnode->up = up;
  thisnode->id = -1; 
  thisnode->branch_length = tree->has_branches ? read_branch_length (rsptr) : 0.;
  thisnode->left = NULL;
  thisnode->right = NULL;
  thisnode->taxlabel = NULL;

  (*node_id)++;

  if (*lsptr == '(') { /* internal node */
    char *newend = rsptr;
    int comma_pos = find_branch_split_newick (lsptr, rsptr);

    thisnode->left = subtree_nexus_tree (tree, lsptr+1, lsptr+comma_pos-1, node_id, thisnode);
    while ((newend != lsptr) && (newend != NULL) && (*newend != ')')) newend--;
    if (newend == lsptr) newend = rsptr;
    thisnode->right = subtree_nexus_tree (tree, lsptr+comma_pos+1, newend-1, node_id, thisnode);
  }

  else thisnode->taxlabel = read_taxlabel (lsptr, rsptr); /* leaf */

  return thisnode;
}

char *
read_taxlabel ( const char *name_start, const char *name_end) 
{
  size_t seqsize;
  size_t i;
  char *tmp, *label=NULL;
  tmp = (char*) name_start;
  while ((tmp <= name_end) && (*tmp != ',') && (*tmp != ')') && (*tmp != ':')) tmp++;
  seqsize = tmp - name_start;
  i = sizeof (char)*(seqsize+1);
  label = (char*) biomcmc_malloc (i);
  label[0] = '\0';
  strncat (label, name_start, seqsize);
  return label;
}

void
create_leaflist_nexus_tree (nexus_tree tree, nexus_node this, int *id) 
{
  if (this->taxlabel != NULL) tree->leaflist[(*id)++] = this;
  else {
    if (this->left)  create_leaflist_nexus_tree (tree, this->left, id);
    if (this->right) create_leaflist_nexus_tree (tree, this->right, id);
  }
}

void
create_node_id_nexus_tree (nexus_node this, int *id) 
{
  this->id = (*id)++;
  if ( (this->left)  && (this->left->id < 0)  ) create_node_id_nexus_tree (this->left, id);
  if ( (this->right) && (this->right->id < 0) ) create_node_id_nexus_tree (this->right, id);
}

double
read_branch_length (char *right_string_ptr) 
{
   char *backwards = right_string_ptr;
   double branch=0.;
  
  if ((*backwards == ')') || (*backwards == ',')) return DEFAULTBLENGTH;
  while (*backwards != ':') backwards--;
  sscanf (backwards, ": %lf", &branch);
  if (branch < 0.) return 0.;
  else return branch;
}

int
find_branch_split_newick (char *left_string_ptr, char *right_string_ptr) 
{
  int nLevel = 0;
  int comma_position = 0;
  char *treeptr = left_string_ptr;
  
  while ((*treeptr != ',') || (nLevel != 1)) { /* stop when *treeptr == ',' AND nLevel == 1 */
    if (treeptr == right_string_ptr) biomcmc_error ("unbalanced tree: couldn't find innermost comma for subtree");
    if (*treeptr == '(') nLevel++;
    if (*treeptr == ')') nLevel--;
    treeptr++;
    comma_position++;
  }
  return comma_position;
}

int
number_of_leaves_in_newick (char **string, bool resolve_trifurcation)
{
  int nopen = 0, nclose = 0, ncommas = 0, i, nsplit = 0;
  int len = strlen (*string);
  int has_branches = 0; /* could be a bool, but I want to debug number of branches */
  char current;

  if (*(*string + len - 1) == ';') *(*string + len - 1) = '\0';
  for (i = 0; i < len; i++) {
    current = (*string)[i];
    if (current == ',' && (nopen - nclose) == 1) {
      if (!nsplit) nsplit = i;
      ncommas++;
    }
    else if (current == '(') nopen++;
    else if (current == ')') nclose++;
    else if (current == ':') has_branches++;
  }
 //FIXME:: bug with notebook 010 
  if (nopen != nclose || ncommas > 2 || ncommas < 1) biomcmc_error ( "Invalid tree structure: %s", *string);
  if (ncommas == 2) nopen++; /* trifurcation (unrooted tree, hopefully) */
  if (!resolve_trifurcation) return nopen + 1; /* do not fix trifurcation */
  
  if (ncommas == 2) { /* try to root an unrooted tree (assuming that the trifurcation means an unrooted tree) */
    char *tstring = NULL;
    int add = 3;
    if (has_branches) add = 7;
    *string = (char *) biomcmc_realloc ((char *) *string, sizeof (char) * (len + add + 1));
    tstring = (char *) biomcmc_malloc (sizeof (char) * (len + add + 1));
    bzero (tstring, sizeof (char) * (len + add + 1));
    tstring = strncat (tstring, *string, nsplit);
    tstring = strcat (tstring, ",("); // replaces "," for ",("
    tstring = strncat (tstring, *string + nsplit + 1, len - nsplit - 1);

    if (has_branches) tstring = strcat (tstring, "):0.0");
    else tstring = strcat (tstring, ")");

    strcpy (*string, tstring);
    free (tstring);
  }
  return nopen + 1;
}

void
store_filename_in_topology_space (topology_space tre, char *filename)
{
  char *first_char = NULL, *last_char = NULL;
  size_t size;

  last_char = strrchr (filename, '.'); /* point to last occurence of "." */
  if (last_char) size = last_char - filename + 1; /* plus one for terminating EOL */
  else           size = strlen (filename) + 1;
  first_char = strrchr (filename, '/'); /* point to last occurence of "/" -- the "beginning" of file name removes all directories */
  if (first_char) size -= (first_char - filename); 
  else            first_char = filename;

  tre->filename = biomcmc_malloc (size * sizeof (char));
  memcpy (tre->filename, first_char, size);
  memcpy (tre->filename + size - 1, "\0", 1);
}

