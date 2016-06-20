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

#include "alignment.h"
#include "nexus_common.h"

#define EPSLON 1.e-3
int char2bit[256][2] = {{0xffff}}; /* DNA base to bitpattern translation, with 1st element set to arbitrary value */
int pairdist[15][15][2]; /* pairwise distance table with #matches, #transitions and #transversions for each pair */

/*! \brief Allocates space for new alignment struct. */
alignment new_alignment (int ntax, int nchar);
/*! \brief Reads one line of alignment (sequence data for one taxa) assuming the taxon name precedes alignment. */
void read_interleaved_nexus_sequence (char *line, alignment align);
/*! \brief Reads one line of alignment (sequence data for one taxa) untill all sites for taxon are read. */
void read_sequential_nexus_sequence (char *line, alignment align);
/*! \brief Reads partition information from alignment. */
void read_assumptions_nexus_sequence (char *line, alignment align);
/*! \brief Reduce site columns to only those with distinct patterns; original sites can be found through
 * alignment_struct::site_pattern */
void alignment_create_sitepattern (alignment align);
/*! \brief If taxon labels have spaces or characters with special meaning for newick trees then
 * alignment_struct::taxshort will have their legal versions; otherwise it will link to taxlabel */
void alignment_shorten_taxa_names (alignment align);
/*! \brief low-level function that copies only allowed characters from taxlabel to taxshort */
void copy_taxlabel_to_shortname (char *big, char *small, size_t size);
/*! \brief store file name (without extension) in alignment_struct->filename */
void store_filename_in_alignment (alignment align, char *seqfilename);

/*! \brief calculates empirical equilibrium site frequencies (site counts) */
void calc_empirical_equilibrium_freqs (char *seq, int *pfreq, int nsites, double *result);
/*! \brief uses Kimura's two-parameter model to calculate distance and ti/tv rate ratio between sequencies */
void calc_pairwise_distance_K2P (char *s1, char *s2, int *w, int nsites, double *result);
/*! \brief initializes char2bit vector (local to this file) A->0001 C->0010 G->0100 T->1000 */
void initialize_char2bit_table (void);

char_vector
new_char_vector_from_file (char *filename)
{
  FILE *infile;
  char *line = NULL, *line_read = NULL, *start = NULL;
  size_t linelength = 0;
  int i;
  char_vector vec = new_char_vector (1);

  /* start reading file */
  infile = biomcmc_fopen (filename, "r");
  while (biomcmc_getline (&line_read, &linelength, infile) != -1) {
    line = remove_nexus_comments (&line_read, &linelength, infile);
    if (nonempty_fasta_line (line)) {
      for (start = line; (*start != '\0') && isspace (*start); start++); /* skip leading spaces */
      if (*start == '\0') biomcmc_error ("found EOL while reading non-empty line in file \"%s\"\n", filename);
      char_vector_add_string (vec, start);
    }
  }
  fclose (infile);
  if (line) free (line);

  for (i = 0; i < vec->nstrings; i++) {
    linelength = strcspn (vec->string[i], "#;"); /* string length until beginning of comments (if any) */
    if (linelength < vec->nchars[i]) *(vec->string[i]+linelength) = '\0'; /* EOL here; space released below */
  }

  /* this function will remove comments (by releasing the memory after EOL) and update nchars[] */
  if (char_vector_remove_empty_strings (vec)) biomcmc_error ("illegal empty lines in file \"%s\"", filename);

  return vec;
}

void
char_vector_longer_first_order (char_vector vec) 
{
  empfreq ef;
  int i, *order;

  order = (int*) biomcmc_malloc (vec->nstrings * sizeof (int));
  ef = new_empfreq_sort_decreasing (vec->nchars, vec->nstrings, 1); /*1 = size_t (0=char, 2=int) */
  
  for (i = 0; i < vec->nstrings; i++) order[i] = ef->i[i].idx;

  char_vector_reorder_strings (vec, order);

  del_empfreq (ef);
  if (order) free (order);
}

alignment
read_alignment_from_file (char *seqfilename)
{
  FILE *seqfile;
  char *line = NULL, *line_read = NULL;
  size_t linelength = 0;
  int is_nexus = 0, i;

  seqfile = biomcmc_fopen (seqfilename, "r");

  /* Search for evidence of nexus file; if obligatory syntax "#NEXUS" is found, see if it is an alignment */
  for (i=0; (is_nexus < 4) && (i < 256) && (biomcmc_getline (&line_read, &linelength, seqfile) != -1);) {
    line = line_read; /* the variable *line_read should point always to the same value (no line++ or alike) */
    if (nonempty_string (line)) { 
      if ((!is_nexus) && strcasestr (line, "#NEXUS")) is_nexus++; 
      else if ((is_nexus == 1) && strcasestr (line, "BEGIN") && strcasestr (line, "DATA")) is_nexus++; 
      else if ((is_nexus == 2) && strcasestr (line, "DIMENSIONS")) is_nexus++; 
      else if ((is_nexus == 3) && strcasestr (line, "MATRIX")) is_nexus++; 
    }
    if (!is_nexus) i++; /* give up if "#NEXUS" is not found in first 256 lines; otherwise keep looking */
  }

  fclose (seqfile);
  if (line_read) free (line_read);

  if (is_nexus == 4) return read_nexus_alignment_from_file (seqfilename);
  else return read_fasta_alignment_from_file (seqfilename);
}

alignment
read_fasta_alignment_from_file (char *seqfilename)
{
  alignment align = NULL; /* we will create a new alignment by hand (no call to new_alignment() */
  FILE *seqfile;
  char *line = NULL, *line_read = NULL, *delim = NULL;
  size_t linelength = 0;
  int i, max_length = 0;

  /* equivalent to call to new_alignment(), but hashtable will be initialized afterwards */
  align = (alignment) biomcmc_malloc (sizeof (struct alignment_struct));
  align->n_charset = 0; /* FASTA files don't support partitioning */ 
  align->charset_start = align->charset_end = NULL;
  align->taxlabel      = new_char_vector (1); /* will increase dynamically */ 
  align->character     = new_char_vector (1); /* will increase dynamically */ 
  align->is_aligned    = true; /* we will check after reading file */
  align->site_pattern = NULL; 
  align->pattern_freq = NULL;
  align->npat = 0; /* number of site patterns only make sense if aligned (checked by alignment_create_sitepattern()) */

  /* start reading file */
  seqfile = biomcmc_fopen (seqfilename, "r");
  while (biomcmc_getline (&line_read, &linelength, seqfile) != -1) {
    line = line_read; /* the variable *line_read should point always to the same value (no line++ or alike) */
    if (nonempty_fasta_line (line)) { /* each line can be either the sequence or its name, on a strict order */
      /* sequence description (in FASTA jargon); the sequence name */
      if ((delim = strchr (line, '>'))) char_vector_add_string (align->taxlabel, ++delim);
      /* the sequence itself, which may span several lines */
      else {
        line = remove_space_from_string (line);
        line = uppercase_string (line);
        char_vector_append_string_at_position (align->character, line, align->taxlabel->next_avail-1);
      }
    }
  }
  fclose (seqfile);

  if (line_read) free (line_read);


  if (align->taxlabel->nstrings != align->character->nstrings)
    biomcmc_error ("number of sequences and number of sequence names disagree in FASTA file \"%s\"\n", seqfilename);
  if (char_vector_remove_empty_strings (align->taxlabel)) /* store sequence length info in taxlabel->nchars[] */
    biomcmc_error ("problem (empty string) reading sequence name for FASTA file \"%s\"\n", seqfilename);
  if (char_vector_remove_empty_strings (align->character))/* store sequence length info in character->nchars[] */
    biomcmc_error ("problem (empty string) reading sequence data for FASTA file \"%s\"\n", seqfilename);

  for (i=0; i < align->character->nstrings; i++) { /* check if is aligned and store size of largest sequence */
    if (max_length && (max_length != (int) align->character->nchars[i])) align->is_aligned = false;
    if (max_length < (int) align->character->nchars[i]) max_length = align->character->nchars[i];
  }

  align->ntax  = align->taxlabel->nstrings;
  align->nchar = max_length;

  align->taxlabel_hash = new_hashtable (align->ntax);
  for (i=0; i < align->ntax; i++) insert_hashtable (align->taxlabel_hash, align->taxlabel->string[i], i);

  if (align->is_aligned) alignment_create_sitepattern (align); /* compact align->character to site patterns */
  if (char2bit[0][0] == 0xffff) initialize_char2bit_table (); /* translation table between ACGT to 1248 */
  alignment_shorten_taxa_names (align);
  store_filename_in_alignment (align, seqfilename);

  return align;
}

alignment
read_nexus_alignment_from_file (char *seqfilename)
{
  alignment align=NULL;
  FILE *seqfile;
  char *line = NULL, *line_read = NULL, *needle_tip;
  bool option_begin_data   = false, 
       option_begin_matrix = false, 
       option_end_matrix   = false, 
       option_interleave   = false,
       option_begin_assumptions  = false,
       option_end_assumptions = false,
       option_read_complete = false;
  size_t linelength = 0;
  int ntax=0, nchar=0;

  seqfile = biomcmc_fopen (seqfilename, "r");

  /* the variable *line_read should point always to the same value (no line++ or alike) */
  biomcmc_getline (&line_read, &linelength, seqfile);
  line = remove_nexus_comments (&line_read, &linelength, seqfile); /* skip lines until all comments are closed */
  while (!nonempty_string (line)) {
    /* skip (possibly not following NEXUS format) initial comments and blank lines */
    if (biomcmc_getline (&line_read, &linelength, seqfile) < 0) 
      biomcmc_error ("Premature end of NEXUS alignment file %s\n", seqfilename);
    line = remove_nexus_comments (&line_read, &linelength, seqfile);
  }
  if (!strcasestr (line, "NEXUS")) 
    biomcmc_error ( "%s is not ot a Nexus sequence file (first line should be \"#NEXUS\")\n", seqfilename);

  while (biomcmc_getline (&line_read, &linelength, seqfile) != -1) {
    line = remove_nexus_comments (&line_read, &linelength, seqfile);
    if (nonempty_string (line)) { 
      /* do not do anything until 'BEGIN DATA' block */
      if ((!option_begin_data) && (strcasestr (line, "BEGIN DATA"))) option_begin_data = true;

      else if (!option_begin_matrix) {
        /* read ntax and nchar info, and wait for 'MATRIX' block */

        if (strcasestr (line, "MATRIX")) option_begin_matrix = true;
        else if (strcasestr (line, "INTERLEAVE")) option_interleave = true;
        else if ( (needle_tip = strcasestr (line, "DIMENSIONS")) ) {
          /* needle_tip will point to first character of 'DIMENSIONS', removing indenting spaces */
          needle_tip += strlen("DIMENSIONS");
          needle_tip = uppercase_string (needle_tip);
          sscanf (needle_tip, " NTAX = %d NCHAR = %d ", &ntax, &nchar);
          align = new_alignment (ntax, nchar);
        } // else if (needle_tip)
      }

      else if (!option_read_complete){
        /* entering 'MATRIX' block with sequence data */
        if ( (needle_tip = strstr (line, ";")) ) {
          /* semi-colon may be at the last line of seq data, or at new one */
          option_read_complete = true;
          if (needle_tip > line) {
            if (option_interleave) read_interleaved_nexus_sequence (line, align);
            else read_sequential_nexus_sequence (line, align);
          }
        } // if (needle_tip = strstr())
        else { 
          if (option_interleave) read_interleaved_nexus_sequence (line, align);
          else read_sequential_nexus_sequence (line, align);
        }
      }

      /* after finish reading matrix data, we need to wait for END keyword */
      else if ((!option_end_matrix) && (strcasestr (line, "END"))) option_end_matrix = true;

      /* wait for assumptions block */
      else if ((!option_begin_assumptions) && (strcasestr (line, "BEGIN ASSUMPTION"))) option_begin_assumptions = true; 

      else if (!option_end_assumptions) {
        /* once in assumptions block, just read CHARSET (gene segments) and order segments after END*/ 
        if (strcasestr (line, "CHARSET") && (needle_tip = strcasestr (line, "="))) {
          needle_tip++;
          read_assumptions_nexus_sequence (needle_tip, align);
        }
        else if (strcasestr (line, "END")) {
          option_end_assumptions = true;
          qsort (align->charset_start, align->n_charset, sizeof (int), compare_int_increasing);
          qsort (align->charset_end, align->n_charset, sizeof (int), compare_int_increasing);
        }
      }

      /* other commands, at the end of file, come here (empty now) */

    } // if (line)
  } //while (biomcmc_getline)
  fclose (seqfile);
  if (line_read) free (line_read);

  if (align->is_aligned) alignment_create_sitepattern (align); /* compact align->character to site patterns */
  if (char2bit[0][0] == 0xffff) initialize_char2bit_table (); /* translation table between ACGT to 1248 */
  alignment_shorten_taxa_names (align);
  store_filename_in_alignment (align, seqfilename);

  return align;
}

void
print_alignment_in_fasta_format (alignment align, FILE *stream)
{
  int i, j, k;
  static int columns = 120;

  if (align->is_aligned) for (i = 0; i < align->ntax; i++) {
    fprintf (stream, ">%s\n", align->taxlabel->string[i]);
    for (j = 0; j < align->nchar; j+= columns) {
      for (k = 0; (k < columns) && ((j+k) < align->nchar); k++) 
        fprintf (stream, "%c", align->character->string[i][ align->site_pattern[j+k] ]);
      fprintf (stream, "\n");
    }
  }

  else for (i = 0; i < align->ntax; i++) {
    fprintf (stream, ">%s\n", align->taxlabel->string[i]);
    for (j = 0; j < (int) align->character->nchars[i]; j+= columns) {
      for (k = 0; (k < columns) && ((j+k) < (int) align->character->nchars[i]); k++) 
        fprintf (stream, "%c", align->character->string[i][j+k]);
      fprintf (stream, "\n");
    }
  }
}

alignment
new_alignment (int ntax, int nchar)
{
  alignment align;

  align = (alignment) biomcmc_malloc (sizeof (struct alignment_struct));
  align->ntax  = ntax;
  align->nchar = nchar;
  align->npat = 0; /* only if sequences are aligned we calculate this value (number of site patterns) */
  align->is_aligned = true; /* default to NEXUS; FASTA files don't use this function new_alignment() */

  align->site_pattern = NULL; 
  align->pattern_freq = NULL;
  align->filename = NULL;

  align->n_charset = 0;
  align->charset_start = NULL; /* (re)allocated by read_assumptions() */
  align->charset_end   = NULL; /* (re)allocated by read_assumptions() */

  align->taxlabel      = new_char_vector (ntax); 
  align->taxshort      = NULL; /* allocated (or points to taxlabel) by alignment_shorten_taxa_names() */ 
  align->character     = new_char_vector_fixed_length (ntax, nchar); 
  align->taxlabel_hash = new_hashtable (ntax);

  return align;
}

void 
del_alignment (alignment align)
{
  if (!align) return;
  del_char_vector (align->character); 
  del_char_vector (align->taxlabel); 
  del_char_vector (align->taxshort); 
  del_hashtable   (align->taxlabel_hash);
  if (align->charset_start) free (align->charset_start);
  if (align->charset_end)   free (align->charset_end);
  if (align->site_pattern)  free (align->site_pattern);
  if (align->pattern_freq)  free (align->pattern_freq);
  if (align->filename)      free (align->filename);
  free (align);
}

void
read_interleaved_nexus_sequence (char *line, alignment align)
{
  char seqname[MAX_NAME_LENGTH]="", *last = NULL;
  size_t namelength; 
  int position;
  if (strchr (line, '\"')) { /* taxon name has double quotes (in which case may have spaces) */
    namelength = strcspn (line+1, "\""); /* how many chars until last char before double quotes */ 
    if (namelength < MAX_NAME_LENGTH) strncpy (seqname, line+1, namelength);
    else biomcmc_error ( "Taxon name too long (more than %d characters)\n", MAX_NAME_LENGTH);
    namelength+= 2; /* how many chars until first char after the double quotes */
  }
  else { /* simple taxon name without spaces */
    namelength = strcspn (line, " \t\n");
    if (namelength < MAX_NAME_LENGTH) strncpy (seqname, line, namelength);
    else biomcmc_error ( "Taxon name too long (more than %d characters)\n", MAX_NAME_LENGTH);
  }
  position = lookup_hashtable (align->taxlabel_hash, seqname);
  if (position < 0) { /* new taxon */
    position = align->taxlabel->next_avail;
    char_vector_add_string (align->taxlabel, seqname); /* creates new string and updates taxlabel->next_avail */
    insert_hashtable (align->taxlabel_hash, seqname, position);
  }
  line = remove_space_from_string (line+namelength); /* modifies (*line), but we still have the original (*line_read) */
  line = uppercase_string (line);
  if ((last = strchr (line, ';'))) (*last) = '\0'; /* read up to this position */
  char_vector_append_string_at_position (align->character, line, position);
}

void
read_sequential_nexus_sequence (char *line, alignment align)
{
  char *last = NULL;
  int position = align->character->next_avail;
  size_t namelength = 0;

  /* sometimes called once more after finishing reading - because of ending ";" */
  if (position >= align->character->nstrings) return;

  if (position == align->taxlabel->next_avail) { /* new taxon */
    char seqname[MAX_NAME_LENGTH]="";
    /* read taxon name */
    if (strchr (line, '\"')) { /* taxon name has double quotes (in which case may have spaces) */
      namelength = strcspn (line+1, "\""); /* strcspn() will give us (the size until) the last element before '"' */
      if (namelength < MAX_NAME_LENGTH) strncpy (seqname, line+1, namelength);
      else biomcmc_error ( "Taxon name too long (more than %d characters)\n", MAX_NAME_LENGTH);
      namelength+= 2; /* how many chars until first char after the double quotes */
    }
    else { /* simple taxon name without spaces */
      namelength = strcspn (line, " \t\n");
      if (namelength < MAX_NAME_LENGTH) strncpy (seqname, line, namelength);
      else biomcmc_error ( "Taxon name too long (more than %d characters)\n", MAX_NAME_LENGTH);
    }
    insert_hashtable (align->taxlabel_hash, seqname, position);
    char_vector_add_string (align->taxlabel, seqname); /* creates new string and updates taxlabel->next_avail */
  }

  line = remove_space_from_string (line+namelength); /* modifies (*line), but we still have (*line_read) */
  line = uppercase_string (line);
  if ((last = strchr (line, ';'))) (*last) = '\0'; /* read up to this position */

  char_vector_append_string_at_position (align->character, line, position);
  if (strlen (align->character->string[position]) == align->character->nchars[position]) 
    align->character->next_avail++; /* finished reading alignment for this taxon */
}

void
read_assumptions_nexus_sequence (char *line, alignment align)
{
  int start, end;

  align->n_charset++;
  align->charset_start = (int*) biomcmc_realloc ((int*) align->charset_start, sizeof (int) * align->n_charset);
  align->charset_end   = (int*) biomcmc_realloc ((int*) align->charset_end,   sizeof (int) * align->n_charset);

  if (sscanf (line, " %d - %d ;", &start, &end) != 2) 
    biomcmc_error ( "problem in ASSUMPTIONS block\n");

  start--; end--;
  /* error handling */
  if (start < 0)             biomcmc_error ( "CHARSET start < 1\n");
  if (start >= align->nchar) biomcmc_error ( "CHARSET start > NCHAR\n");
  if (end <= start)          biomcmc_error ( "CHARSET end <= CHARSET start\n");
  if (end >= align->nchar)   biomcmc_error ( "CHARSET end > NCHAR\n");

  align->charset_start[align->n_charset-1] = start;
  align->charset_end[align->n_charset-1] = end;
}

void
//alignment_create_sitepattern_OLD_implementation (alignment align)
alignment_create_sitepattern (alignment align)
{
  /* faster version: whenever we find a duplicate, we swap the column with the last one; even with this extra copying we
   * have O(nchar * ln(nchar)) comparisions */
  int s1, s2, seq, *index, nchar = align->nchar;
  bool equal;

  /* vector char_vector::nchars is redundant for aligned sequences (all have same value) */
  if (align->character->nchars) free (align->character->nchars);
  align->character->nchars = NULL;

  /* only aligned sequences have vector site_pattern (if one needs original pattern at position */
  align->site_pattern = (int *) biomcmc_malloc (nchar * sizeof (int));
  index = (int *) biomcmc_malloc (nchar * sizeof (int)); 

  for (s1=0; s1 < nchar; s1++) index[s1] = align->site_pattern[s1] = s1;

  for (s1 = 0; s1 < nchar- 1; s1++) {
    for (s2 = s1 + 1; s2 < nchar; s2++) {
      /* compare site columns s1 and s2 for same pattern (=same state for all taxa) */
      for (equal = true, seq = 0; (equal == true) && (seq < align->ntax); seq++) 
        if (align->character->string[seq][s1] != align->character->string[seq][s2]) equal = false;

      if (equal == true) { /* s1 and s2 have identical patterns */
        align->site_pattern[ index[s2] ] = s1; /* map between compressed (just patterns) and original columns */
        nchar--;
        if (s2 < nchar) { /* replace pattern s2 by pattern of last site under analysis */
          index[s2] = index[nchar];
          for (seq = 0; seq < align->ntax; seq++) 
            align->character->string[seq][s2] = align->character->string[seq][nchar];
        }
        s2--;/* compare again, since character[][s2] now contains last site */
      }// if (equalpattern)
    }// for (s2)
    /* in some cases (next)idx[s1] = (prev)idx[s2] = idx[nchar] */
    if (align->site_pattern[ index[s1] ] > s1) align->site_pattern[ index[s1] ] = s1; 
  }// for (s1)
  /* in some cases (next)idx[s1] = (prev)idx[s2] = idx[nchar] */
  if (align->site_pattern[ index[s1] ] > s1) align->site_pattern[ index[s1] ] = s1; 

  align->npat = nchar; /* since char_vector::nchars is obsolete we need to store the number of patterns here */

  if (nchar < align->nchar) { /* we have more than one site column with same pattern */
    for (seq = 0; seq < align->character->nstrings; seq++) { /* compress char_vector to store only patterns (unique) */
      align->character->string[seq] = 
      (char *) biomcmc_realloc ((char*) align->character->string[seq], (nchar + 1) * sizeof (char));
      align->character->string[seq][nchar] = '\0';
    }
  }

  /* calculate frequency of each pattern */
  align->pattern_freq = (int *) biomcmc_malloc (nchar * sizeof (int)); 
  for (s1=0; s1 < nchar; s1++) align->pattern_freq[s1] = 0;
  for (s1=0; s1 < align->nchar; s1++) align->pattern_freq[ align->site_pattern[s1] ]++;

  if (index) free (index);
}

void
alignment_create_sitepattern_NEW_implementation (alignment align)
//alignment_create_sitepattern (alignment align)
{ 
  /* alternative implementation: valid[] vector s.t. columns copying is done only at the end, slower since we scan
   * nchar^2 elements (even though the duplicate ones we skip) */
  int s1, s2, seq, npat = 0, *valid;
  bool equal;

  /* valid (unique) sites */
  valid = (int *) biomcmc_malloc (align->nchar * sizeof (int));

  /* vector char_vector::nchars is redundant for aligned sequences (all have same value) */
  if (align->character->nchars) free (align->character->nchars);
  align->character->nchars = NULL;

  /* only aligned sequences have vector site_pattern (if one needs original pattern at position */
  align->site_pattern = (int *) biomcmc_malloc (align->nchar * sizeof (int));

  for (s1=0; s1 < align->nchar; s1++) align->site_pattern[s1] = s1;

  /* fill align->site_pattern[] with id of distinct patterns */
  for (s1 = 0; s1 < align->nchar - 1; s1++) if (align->site_pattern[s1] == s1) {
    valid[npat] = s1;
    align->site_pattern[s1] = npat++; /* number of distinct patterns */
    for (s2 = s1 + 1; s2 < align->nchar; s2++) if (align->site_pattern[s2] == s2) {
      /* compare site columns s1 and s2 for same pattern (=same state for all taxa) */
      for (equal = true, seq = 0; (equal == true) && (seq < align->ntax); seq++) 
        if (align->character->string[seq][s1] != align->character->string[seq][s2]) equal = false;

      if (equal == true) align->site_pattern[s2] = align->site_pattern[s1]; /* s1 and s2 have identical patterns */
    }
  }
  if (align->site_pattern[s1] == s1) align->site_pattern[s1] = npat++; /* last column may be unique */

  align->npat = npat; /* since char_vector::nchars is obsolete we need to store the number of patterns here */

  if (npat < align->nchar) { /* we have more than one site column with same pattern */
    /* copy columns: character->string[] will store only distinct site columns */
    for (s1 = 0; s1 < npat; s1++) for (seq = 0; seq < align->character->nstrings; seq++)
      align->character->string[seq][s1] = align->character->string[seq][ valid[s1] ];
    /* compress char_vector to store only patterns (unique) */
    for (seq = 0; seq < align->character->nstrings; seq++) {
      align->character->string[seq] = 
      (char *) biomcmc_realloc ((char*) align->character->string[seq], (npat + 1) * sizeof (char));
      align->character->string[seq][npat] = '\0';
    }
  }

  if (valid) free (valid);

  /* calculate frequency of each pattern */
  align->pattern_freq = (int *) biomcmc_malloc (npat * sizeof (int)); 
  for (s1=0; s1 < npat; s1++) align->pattern_freq[s1] = 0;
  for (s1=0; s1 < align->nchar; s1++) align->pattern_freq[ align->site_pattern[s1] ]++;
}

void
alignment_shorten_taxa_names (alignment align)
{
  /* Specially in fasta files, the short version of a sequence name (e.g. "Ecoli") is followed by a longer description
   * (like "Escherichia coli W3110"); in such cases we must have both versions available: the short for printing the
   * tree in newick format and the long for finding the species tree (substring match) or printing the alignment. */
  int i, j, same_size = 1;
  size_t *short_size, max_size = 0;
  char *shortname;

  short_size = (size_t*) biomcmc_malloc (align->taxlabel->nstrings * sizeof (size_t));

  for (i=0; i < align->taxlabel->nstrings; i++) { 

    /* skip possible leading spaces (in FASTA files usually, since NEXUS taxnames are already free of leading spcs */
    j = strspn (align->taxlabel->string[i], " \t"); /* strspn returns the size_t chars belonging to ' ' or '\t' */
    // for (j = 0; (j < (int) align->taxlabel->nchars[i]) && (align->taxlabel->string[i][j] == ' '); j++); /* OLD */

    /* calc label size until a space (or semicolon etc.), with leading spaces neglected */
    short_size[i] = strcspn (align->taxlabel->string[i]+j, " \t;|") + j;
    /* same_size is zero if at least one label has space */
    same_size = same_size && (short_size[i] == align->taxlabel->nchars[i]);
    /* same_size is zero if at least one label has characters forbidden in newick-formated trees */
    same_size = same_size && (short_size[i] >  strcspn (align->taxlabel->string[i], " (),:#"));

    if (short_size[i] > max_size) max_size = short_size[i];
  }

  if (same_size) { /* no spaces found on sequences */
    if (short_size) free (short_size);
    align->taxshort = align->taxlabel;
    align->taxlabel->ref_counter++;
    return;
  }

  /* space found - taxshort will be a shortened version of taxlabel */
  align->taxshort = new_char_vector (align->taxlabel->nstrings);

  if (!max_size) biomcmc_error ("malformed taxon label (maybe starting with more than one space or \";\" )");
  shortname = (char*) biomcmc_malloc ((max_size+1) * sizeof (char));

  for (i=0; i < align->taxlabel->nstrings; i++) {
    copy_taxlabel_to_shortname (align->taxlabel->string[i], shortname, short_size[i]);
    char_vector_add_string_at_position (align->taxshort, shortname, i);
  }

  if (shortname)  free (shortname);
  if (short_size) free (short_size);
}

void
copy_taxlabel_to_shortname (char *big, char *small, size_t size)
{
  int i, j = 0;
  for (i=0; i < (int) size; i++) /* skip symbols forbidden in newick trees (main use of taxshort, after all) */ 
    if ((big[i] != '(') && (big[i] != ')') && (big[i] != ' ') && 
        (big[i] != ':') && (big[i] != '#') && (big[i] != ',')) small[j++] = big[i];

  small[j] = '\0';  
}

void
store_filename_in_alignment (alignment align, char *seqfilename)
{
  char *last_char = NULL;
  size_t size, stringlength = strlen (seqfilename);
  bool needfix = false; /* weird case where alignment file has extension ".tre" */

  last_char = strrchr (seqfilename, '.'); /* point to last occurence of "." */
  if (last_char) {
    size = last_char - seqfilename + 1; /* plus one for terminating EOL */
    if ((stringlength - size == 3) && strstr (last_char, "tre")) { size++; needfix = true; }
  }
  else size = stringlength + 1;
 
  align->filename = biomcmc_malloc (size * sizeof (char));
  memcpy (align->filename, seqfilename, size);
  if (needfix) memcpy (align->filename + size - 1, "_\0", 2); /* filename has extension ".tre" which is created by program */
  else         memcpy (align->filename + size - 1, "\0", 1); 
} /* FIXME: the solution to a guenomu-specific problem (name overlap) should NOT be fixed here */

distance_matrix
new_distance_matrix_from_valid_matrix_elems (distance_matrix original, int *valid, int n_valid)
{
  int i, j;
  distance_matrix new = new_distance_matrix (n_valid);
  for (i = 0; i < n_valid; i++) for (j = 0; j < n_valid ; j++) new->d[i][j] = original->d[ valid[i] ][ valid[j] ];

  return new;
}

distance_matrix
new_distance_matrix_from_alignment (alignment align)
{
  distance_matrix dist;
  int i, j;
  double result[4], s1 = 0., s2 = 0., jc_proportion, 
         count = 0., this_r, this_d, delta_d, delta_r; /* online calculation of mean and var */

  if (!align->is_aligned) biomcmc_error ("pairwise distances can be calculated only for aligned sequences");
  if (align->character->nstrings < 2) biomcmc_error ("must have at least two sequences to calculate distances");
  dist = new_distance_matrix (align->character->nstrings);

  /*     Count empirical base frequencies and initialize site weights */

  for (i=0; i < 4; i++) result[i] = 0.;
  for (i=0; i < dist->size; i++) /* result[] will accumulate weighted site counts */
    calc_empirical_equilibrium_freqs (align->character->string[i], align->pattern_freq, align->npat, result);
  for (i=0; i < 4; i++)  s1 += result[i]; 
  for (i=0; i < 4; i++) dist->freq[i] = result[i]/s1; /* we may have roundoff errors */


  /*    Kimura's two-parameter distance; result[] will hold partial counts (which will also be used in JC distance) */

  for (i=1; i < dist->size; i++) for (j=0; j < i; j++) {
    calc_pairwise_distance_K2P (align->character->string[i], align->character->string[j], align->pattern_freq,
                                align->npat, result);

    jc_proportion = result[0] + result[1]; /* total proportion of differences (ti+tv) used in Jukes-Cantor formula */

    if (jc_proportion) { /* sequences are not identical */
      /* K2P distance must return real (not complex) numbers */
      if (result[1] >= (0.5 - EPSLON)) result[1] = 0.5 - EPSLON; /* Q < 1/2 */
      if ((2*result[0] +result[1]) >= (1. - EPSLON)) result[0] = 0.5 * (1. - result[1] - EPSLON); /* 2P+Q < 1 */

      /* calculation of distance using K2P (also called K80) formula (JMolecEvol 1999, p274; Felsenstein book 2004) */
      s1 = log (1.- 2.*result[0] - result[1]);
      s2 = log (1. - 2.* result[1]);
      dist->d[j][i] = -0.5 * s1 - 0.25 * s2; /* upper triangular matrix will hold K2P pairwise distance */

      /* online mean and variance calculation */
      count += 1.;
      this_d = dist->d[j][i];       /* average K2P distance over all sequences */
      if (s2) this_r = s1/s2 - 0.5; /* average value of ti/tv rate ratio (=alpha/(2*beta) */
      else this_r = 40;             /* large number (if all subst. are transitions) */
      delta_d = this_d - dist->mean_K2P_dist; 
      dist->mean_K2P_dist += delta_d / count; 
      dist->var_K2P_dist += delta_d * (this_d - dist->mean_K2P_dist); 
      delta_r = this_r - dist->mean_R;
      dist->mean_R += delta_r / count;
      dist->var_R += delta_r * (this_r - dist->mean_R);

      /* calculation of distance using Jukes-Cantor formula (Z. Yang book 2006; Felsenstein book 2004) */
      if (jc_proportion >= (0.75 - EPSLON)) jc_proportion = 0.75 - EPSLON;
      s1 = 1. - (4. * jc_proportion)/3.;
      dist->d[i][j] = -0.75 * log (s1);      /* lower triangular matrix will hold JC pairwise distance*/
      dist->mean_JC_dist += dist->d[i][j];   /* average JC distance over all sequences */ 
    }

    else { /* sequences are identical */
      dist->d[j][i] = dist->d[i][j] = 0.;
    }
  }

  dist->mean_JC_dist  *=  2./(double)(dist->size * (dist->size-1));
  dist->var_K2P_dist /= count - 1.; /* online algortihm ... */ 
  dist->var_R        /= count - 1.; /* ... means are already calculated */

  return dist;
}

void
calc_empirical_equilibrium_freqs (char *seq, int *pfreq, int nsites, double *result)
{
  int i, j, bit, ns;

  for (i=0; i < nsites; i++) {
    bit = char2bit[(int)seq[i]][0]; /* integer representation of char (A=0001 C=0010 G=0100 T=1000) */
    ns  = char2bit[(int)seq[i]][1]; /* how many ambiguous states (e.g. A->1, R(=AG)=2, N=4 etc.) */
    /* increase count of result[] weighted by site_pattern/(#ambiguous) for non-indels (indel: ns == 0)*/
    for (j=0; j < 4; j++) if (ns && ((bit >> j) & 1U)) result[j] += (double)(pfreq[i])/(double)(ns); 
  }
}

void
calc_pairwise_distance_K2P (char *s1, char *s2, int *w, int nsites, double *result)
{
  int i, b1, b2;
  double weight = 1., degeneracy, valid_sites = 0.;

  result[0] = result[1] = 0.;
  for (i=0; i < nsites; i++) {
    /* integer (bit) representation of site states */
    b1  = char2bit[ (int)s1[i] ][0]; b2  = char2bit[ (int)s2[i] ][0];
    degeneracy  = (double) (char2bit[ (int)s1[i] ][1] *  char2bit[ (int)s2[i] ][1]);

    /* exclude indels (char2bit[][0]=15 and char2bit[][1]=0) and eventually illegal chars */
    if (!degeneracy) continue; 

    /* frequency of pattern divided by number of ambiguous states */
    weight = (double)(w[i]) / degeneracy;
    valid_sites += (double)(w[i]); /* without indels this should be equal to align->nchar */

    result[0] += (double)(pairdist[b1-1][b2-1][0]) * weight; /* number of transitions */
    result[1] += (double)(pairdist[b1-1][b2-1][1]) * weight; /* number of transversions */
  }
  if (valid_sites) {
    result[0] /= valid_sites; /* result[] now has total counts (per sequence) but we want */ 
    result[1] /= valid_sites; /* fraction per site (between zero and one) */
    if (!result[1]) result[1] = 0.5/valid_sites; /* must be larger than zero to avoid NaN */
  }
  else result[0] = result[1] = 1.;
}

void
initialize_char2bit_table (void)
{
  int i, j, k, l;
  for (i = 0; i < 256; i++) char2bit[i][0] = char2bit[i][1] = 0;
  /* The ACGT is PAUP convention (and maybe DNAml, fastDNAml); PAML uses TCAG ordering */
  char2bit['A'][0] = 1;   char2bit['A'][1] = 1;  /* .   A */ /* 0001 */
  char2bit['B'][0] = 14;  char2bit['B'][1] = 3;  /* .TGC  */ /* 1110 */
  char2bit['C'][0] = 2;   char2bit['C'][1] = 1;  /* .  C  */ /* 0010 */
  char2bit['D'][0] = 13;  char2bit['D'][1] = 3;  /* .TG A */ /* 1101 */
  char2bit['G'][0] = 4;   char2bit['G'][1] = 1;  /* . G   */ /* 0100 etc.*/
  char2bit['H'][0] = 11;  char2bit['H'][1] = 3;  /* .T CA */ 
  char2bit['K'][0] = 12;  char2bit['K'][1] = 2;  /* .TG   */
  char2bit['M'][0] = 3;   char2bit['M'][1] = 2;  /* .  CA */
  char2bit['N'][0] = 15;  char2bit['N'][1] = 4;  /* .TGCA */
  char2bit['O'][0] = 15;  char2bit['O'][1] = 4;  /* .TGCA */
  char2bit['R'][0] = 5;   char2bit['R'][1] = 2;  /* . G A */
  char2bit['S'][0] = 6;   char2bit['S'][1] = 2;  /* . GC  */
  char2bit['T'][0] = 8;   char2bit['T'][1] = 1;  /* .T    */
  char2bit['U'][0] = 8;   char2bit['U'][1] = 1;  /* .T    */
  char2bit['V'][0] = 7;   char2bit['V'][1] = 3;  /* . GCA */
  char2bit['W'][0] = 9;   char2bit['W'][1] = 2;  /* .T  A */
  char2bit['X'][0] = 15;  char2bit['X'][1] = 4;  /* .TGCA */
  char2bit['Y'][0] = 10;  char2bit['Y'][1] = 2;  /* .T C  */
  char2bit['?'][0] = 15;  char2bit['?'][1] = 4;  /* .TGCA */
  char2bit['-'][0] = 15;  char2bit['-'][1] = 0;  /* .TGCA */ /* the "size" of zero makes it skip distance calcs */

  /* distance between patterns (i+1) and (j+1) (e.g. (i+1)=0001=A and (j+1)=0010=C )
   * "static" table to speed up calculations */
  for (i=1; i < 15; i++) for (j=0; j <=i; j++) { /* 16x16 matrix with all bases (ambiguous included) */
    pairdist[i][j][0] = pairdist[i][j][1] = pairdist[j][i][0] = pairdist[j][i][1] = 0;

    for (k=0; k<4; k++) for (l=0;l<4;l++) /* 4x4 matrix with all bit pairs */
      if ( (k!=l) && (((i+1)>>k)&1U) && (((j+1)>>l)&1U) ) { /* both bits are sets between (i+1) and (j+1) */
        /* AC=(0,1) GC=(2,1) AT=(0,3) GT=(2,3) are transversions */
        if ((k+l)%2) { pairdist[i][j][1]++; if (i!=j) pairdist[j][i][1]++; } 
        /* AG=(0,2) CT=(1,3) are transitions. Note that (2,0) and (3,1) are also possible */
        else { pairdist[i][j][0]++; if (i!=j) pairdist[j][i][0]++; } 
        /* OBS: the number of matches (when l==k) can be found by (pseudocode)
         * charbit[seq1][1]*charbit[seq2][1] - pairdist[0] - pairdist[1]
         * since all combinations (size1*size2) must be OR a match OR a transition OR a transversion */
      }
  }
}

void 
store_likelihood_info_at_leaf (double **l, char *align, int n_pat, int n_state)
{
  int i, j; /*the calling function should check if char2bit is initialized or not... */ 
  for (j = 0; j < n_pat; j++) for (i=0; i < n_state; i++) l[j][i] = 0.;
  for (j = 0; j < n_pat; j++) for (i=0; i < n_state; i++)
    if (char2bit[ (int)align[j] ][0] & (1 << i)) l[j][i] = 1.;
}
