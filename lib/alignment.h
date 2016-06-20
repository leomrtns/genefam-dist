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

/*! \file alignment.h 
 *  \brief File handling functions and calculation of distances for sequence data in nexus format.
 *
 *  Reading of sequence data in nexus format (sequencial or interleaved) and fasta format. For fasta format the
 *  sequences don't need to be aligned, but for all formats if the sequences are aligned a data compression is used so
 *  that we keep only the distinct site (column) patterns and a mapping between original and compressed site columns.
 *  Based on the sequence pairs we can also calculate the matrix of distances between sequences.
 */

#ifndef _biomcmc_alignment_h_
#define _biomcmc_alignment_h_

#include "hashtable.h"
#include "empirical_frequency.h"

typedef struct alignment_struct* alignment;

/*! \brief Data from alignment file. */
struct alignment_struct
{
  int ntax, nchar, npat;   /*! \brief Number of species, sites and patterns according to sequence file. */
  char_vector character;   /*! \brief Vector with aligned sequence for each taxon. */
  char_vector taxlabel;    /*! \brief Taxon names from file. */
  char_vector taxshort;    /*! \brief Alias (short version) for taxon names that can be used in newick trees. */
  hashtable taxlabel_hash; /*! \brief Lookup table with taxon names. */
  int n_charset;                    /*! \brief Number of gene segments (ASSUMPTIONS BLOCK). */
  int *charset_start, *charset_end; /*! \brief Start and end of each gene segment (from 1...NCHAR) (ASSUMPTIONS ). */
  bool is_aligned;         /*! \brief FASTA files don't need to be aligned; NEXUS files do. */
  int *site_pattern;       /*! \brief pattern, in alignment_struct::character, to which original site belongs. */
  int *pattern_freq;       /*! \brief if sequences are aligned, this is the frequency of each pattern. */
  char *filename;          /*! \brief name of the original file, with extension removed */
};


/*! \brief General function that stores file content into char_vector_struct, removing shell-type comments */
char_vector new_char_vector_from_file (char *filename);
/*! \brief Order the elements of char_vector_struct from longer string to smaller; can be used after calling
 * new_char_vector_from_file() but not on topology-associated char_vectors since other structures may depend on current
 * ordering (like alignment or tree leaves) */
void char_vector_longer_first_order (char_vector vec);

/*! \brief Reads DNA alignment (guess format between FASTA and NEXUS) from file and store info in alignment_struct. */
alignment read_alignment_from_file (char *seqfilename);
/*! \brief Reads DNA FASTA alignment from file and store info in alignment_struct. */
alignment read_fasta_alignment_from_file (char *seqfilename);
/*! \brief Reads DNA NEXUS alignment from file and store info in alignment_struct. */
alignment read_nexus_alignment_from_file (char *seqfilename);
/*! \brief Prints alignment to FILE stream in FASTA format (debug purposes). */
void print_alignment_in_fasta_format (alignment align, FILE *stream);
/*! \brief Frees memory from alignment_struct. */
void del_alignment (alignment align);

/*! \brief new matrix of pairwise distance by simply excluding original elements not present in valid[] */
distance_matrix new_distance_matrix_from_valid_matrix_elems (distance_matrix original, int *valid, int n_valid);
/*! \brief creates and calculates matrix of pairwise distances based on alignment */
distance_matrix new_distance_matrix_from_alignment (alignment align);

/*! \brief transform aligned sequence into likelihood for terminal taxa (e.g. A -> 0001, C-> 0010 etc) (e.g. A -> 0001,
 * C-> 0010 etc) (e.g. A -> 0001, C-> 0010 etc) (e.g. A -> 0001, C-> 0010 etc) (e.g. A -> 0001, C-> 0010 etc) (e.g. A ->
 * 0001, C-> 0010 etc) (e.g. A -> 0001, C-> 0010 etc) (e.g. A -> 0001, C-> 0010 etc) */
void store_likelihood_info_at_leaf (double **l, char *align, int n_pat, int n_state);

#endif

