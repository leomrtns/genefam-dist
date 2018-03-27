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

/*! \file lowlevel.c 
 *  \brief Lowest level basic functions, that should be available to all other modules. 
 */

#include "lowlevel.h"

/* error-safe memory allocation functions */
void *
biomcmc_malloc (size_t size)
{
  void *value = malloc (size);
  if (value == NULL) biomcmc_error ( "biomcmc_malloc error allocating %d bites", size);
  return value;
}

void *
biomcmc_realloc (void *ptr, size_t size)
{
  void *value = (void *) realloc ((void *)ptr, size);
  if (value == NULL) biomcmc_error ( "biomcmc_realloc error on pointer 0x%08X of %d bites\n", ptr, size);
  return value;
}

FILE *
biomcmc_fopen (const char *path, const char *mode)
{
  FILE *fp = fopen (path, mode);
  if (fp == NULL) {
    fprintf (stderr, "Please check if path is correct, if there are non-ASCII characters in file name,\n");
    fprintf (stderr, "if you have enough permissions (to read/write). Remember that paths are relative to\n");
    fprintf (stderr, "where this program is being called\n");
    biomcmc_error ( "problem opening file \"%s\" with mode \"%s\"", path, mode);
  }
  return fp;
}

void
biomcmc_error (const char *template, ...)
{
  va_list ap;

  fprintf (stderr, "genefam error: ");
  va_start (ap, template);
  vfprintf (stderr, template, ap);
  va_end (ap);
  fprintf (stderr, "\n");
  fprintf (stderr, "[note to developers] If you want to debug, set a breakpoint on function biomcmc_error()\n");
  fflush (stderr);
  exit (EXIT_FAILURE);
}

int
compare_int_increasing (const void *a, const void *b)
{
  return (*(int *) a - *(int *) b);
}

int
compare_int_decreasing (const void *a, const void *b)
{
  return (*(int *) b - *(int *) a);
}

int
compare_double_increasing (const void *a, const void *b)
{
  if ((double *) a > (double *) b) return 1;
  if ((double *) a < (double *) b) return -1;
  return 0;
}

int
compare_double_decreasing (const void *a, const void *b)
{
  if ((double *) a < (double *) b) return 1;
  if ((double *) a > (double *) b) return -1;
  return 0;
}

/* \brief size, in bytes, when extending the buffer of biomcmc_getline() */
#define MIN_CHUNK 64
int
biomcmc_getline (char **lineptr, size_t *n, FILE *stream)
{
  int nchars_avail;    /* Allocated but unused chars in *LINEPTR.  */
  char *read_pos;      /* Where we're reading into *LINEPTR. */

  if (!lineptr) biomcmc_error ("NULL pointer sent to biomcmc_getline() as target string");
  if (!n)       biomcmc_error ("string length unavailable to biomcmc_getline()");
  if (!stream)  biomcmc_error ("lack of input file in biomcmc_getline()");

  if (!(*lineptr)) {
    *n = MIN_CHUNK;
    *lineptr = (char *) biomcmc_malloc (*n);
  }
  nchars_avail = *n;
  read_pos = *lineptr;

  for (;;) {
    register int c = getc (stream);

    /* We always want at least one char left in the buffer, since we always (unless we get an error while reading the 
     * first char) NUL-terminate the line buffer.  
     */
    if ((*lineptr + *n) != (read_pos + nchars_avail)) 
      biomcmc_error ("problem setting string size in biomcmc_getline()");
    if (nchars_avail < 2) {
      if (*n > MIN_CHUNK) (*n) *= 2;
      else (*n) += MIN_CHUNK;

      nchars_avail = *n + *lineptr - read_pos;
      *lineptr = (char *) biomcmc_realloc ((char*) *lineptr, *n);
      read_pos = *n - nchars_avail + *lineptr;
      if ((*lineptr + *n) != (read_pos + nchars_avail)) 
        biomcmc_error ("problem setting string size in biomcmc_getline()");
    }

    if (ferror (stream)) return -1;

    if (c == EOF) {
      /* Return partial line, if any.  */
      if (read_pos == *lineptr) return -1;
      else break;
    }

    if (c == '\r') c = '\n';
    *read_pos++ = c; 
    nchars_avail--;
    if (c == '\n') break;
  }

  /* Done - NUL terminate and return the number of chars read.  */
  *read_pos = '\0';
  return (read_pos - (*lineptr));
}

char_vector
new_char_vector (int nstrings)
{
  char_vector vec;
  int i;

  if (nstrings < 1) biomcmc_error ("Vector of strings should have at least one string");

  vec = (char_vector) biomcmc_malloc (sizeof (struct char_vector_struct));
  vec->nstrings = nstrings;
  vec->ref_counter = 1;
  vec->next_avail = 0;

  vec->string = (char**)  biomcmc_malloc (nstrings * sizeof (char*));
  vec->nchars = (size_t*) biomcmc_malloc (nstrings * sizeof (size_t));

  for (i=0; i < nstrings; i++) {
    vec->string[i] = (char*) biomcmc_malloc (sizeof (char));
    vec->string[i][0] = '\0';
    vec->nchars[i] = 0;
  }
  return vec;
}

char_vector
new_char_vector_fixed_length (int nstrings, int nchars)
{
  char_vector vec;
  int i;

  if (nstrings < 1) biomcmc_error ("Vector of strings should have at least one string");

  vec = (char_vector) biomcmc_malloc (sizeof (struct char_vector_struct));
  vec->nstrings = nstrings;
  vec->ref_counter = 1;
  vec->next_avail = 0;

  vec->string = (char**)  biomcmc_malloc (nstrings * sizeof (char*));
  vec->nchars = (size_t*) biomcmc_malloc (nstrings * sizeof (size_t));

  for (i=0; i < nstrings; i++) {
    vec->string[i] = (char*) biomcmc_malloc ((nchars+1) * sizeof (char)); //ending '\0'
    vec->string[i][0] = '\0';
    vec->nchars[i] = (size_t) nchars;
  }
  return vec;
}

void
del_char_vector (char_vector vec)
{
  int i;
  if (!vec) return;
  if (--vec->ref_counter) return;
  if (vec->string) {
    for (i=vec->nstrings-1; i >=0; i--) if (vec->string[i]) free (vec->string[i]);
    free (vec->string);
  }
  if (vec->nchars) free (vec->nchars);
  free (vec);
}

void
char_vector_link_string_at_position (char_vector vec, char *string, int position)
{
  if (position >= vec->nstrings) char_vector_expand_nstrings (vec, position+1);
  if (vec->nchars[position]) free (vec->string[position]);

  vec->nchars[position] = strlen (string); /* Actually alloc'ed memory may be larger than this (next line fix it) */
  string = (char*) biomcmc_realloc ((char*)string, (vec->nchars[position]+1) * sizeof (char));
  vec->string[position] = string;
  
  // TODO: next_avail may be before position; should we assume char_vector is always increasing?
  vec->next_avail = position+1;
}

void
char_vector_add_string_at_position (char_vector vec, char *string, int position)
{
  size_t l;

  string = string + strspn (string, " \t"); /* skip leading spaces */
  l = strlen (string);

  if (!l) return; /* do nothing with empty strings - like last line of some nexus alignments */
  if (position >= vec->nstrings) char_vector_expand_nstrings (vec, position+1);

  if (l > vec->nchars[position]) {
    vec->nchars[position] = l;
    vec->string[position] = (char*) biomcmc_realloc ((char*)vec->string[position], (l+1) * sizeof (char));
  }
  strncpy (vec->string[position], string, l+1); /* l+1 will insert the ending null */
  vec->next_avail = position+1;
}

void
char_vector_add_string (char_vector vec, char *string)
{
  char_vector_add_string_at_position (vec, string, vec->next_avail);
}

void
char_vector_append_string_at_position (char_vector vec, char *string, int position)
{
  size_t l, this_l;
  
  string = string + strspn (string, " \t"); /* skip leading spaces */
  l = strlen (string);

  if (!l) return; /* do nothing with empty strings - like last line of some nexus alignments */
  if (position < 0) position = 0; /* vec->next_avail works for add_string() but not here... */

  if (position >= vec->nstrings) char_vector_expand_nstrings (vec, position+1);
  this_l = strlen (vec->string[position]); /* we assume there is an ending null at vec->string  */

  if ((l+this_l) > vec->nchars[position]) {
    vec->nchars[position] = l + this_l;
    vec->string[position] = (char*) biomcmc_realloc ((char*)vec->string[position], (l + this_l + 1) * sizeof (char));
  }
  strncpy (vec->string[position] + this_l, string, l + 1); /* l+1 will insert the ending null */
}

void
char_vector_append_string (char_vector vec, char *string)
{
  char_vector_append_string_at_position (vec, string, vec->next_avail-1);
}

void
char_vector_expand_nstrings (char_vector vec, int new_size)
{
  int i;

  if (new_size < vec->nstrings) biomcmc_error ("I refuse to reduce char_vector size. This is a bug/feature.");

  vec->string = (char**)  biomcmc_realloc ((char**)  vec->string, new_size * sizeof (char*));
  vec->nchars = (size_t*) biomcmc_realloc ((size_t*) vec->nchars, new_size * sizeof (size_t));

  for (i=vec->nstrings; i < new_size; i++) {
    vec->string[i] = (char*) biomcmc_malloc (sizeof (char));
    vec->string[i][0] = '\0';
    vec->nchars[i] = 0;
  }
  
  vec->nstrings = new_size;
}

void
char_vector_reorder_strings (char_vector vec, int *order)
{
  char  **ptr_c; // tmp pointer to keep track of labels (otherwise once t[2]=t[1] we loose old t[2]...)
  size_t *ptr_i; // tmp pointer to keep track of allocated memory inside char_vector_struct::
  int i;

  if (!vec->next_avail) return; /* if char_vector doesn't have any elements yet */

  ptr_c = (char**)  biomcmc_malloc (vec->nstrings * sizeof (char*));
  ptr_i = (size_t*) biomcmc_malloc (vec->nstrings * sizeof (size_t));

  /* Vector version of tmp = a; a = b; b = tmp; [geek curiosity: how would a "XOR swap" work in this case?] */
  for (i=0; i < vec->nstrings; i++) { 
    ptr_c[i] = vec->string[i]; 
    ptr_i[i] = vec->nchars[i];
  }
  for (i=0; i < vec->nstrings; i++) { 
    vec->string[i] = ptr_c[order[i]]; 
    vec->nchars[i] = ptr_i[order[i]];
  }
  if (ptr_c) free (ptr_c);
  if (ptr_i) free (ptr_i);
}

int
char_vector_remove_empty_strings (char_vector vec)
{
  int i, n_invalid = 0, n_valid=0, *valid;
  size_t length;

  if (!vec->next_avail) return 0; /* if char_vector doesn't have any elements yet */

  /* vector with indexes of valid (non-empty) strings */
  valid = (int*) biomcmc_malloc (vec->nstrings * sizeof (int));
  for (i=0; i < vec->nstrings; i++) { /* find non-empty strings and realloc them */
    length = strlen (vec->string[i]);
    if (!length) { /* empty string - this shouldn't happen and will be reported through the return value */
      n_invalid++;
      if (vec->string[i]) free (vec->string[i]); /* the condition is always true... */
    }
    else { /* valid string */
      valid[n_valid++] = i;
      if (length < vec->nchars[i]) { /* we can save some space, since nchars gives us the allocated memory */
        vec->nchars[i] = length;
        vec->string[i] = (char*) biomcmc_realloc ((char*) vec->string[i], (length+1) * sizeof (char));
      }
    }
  }

  if (!n_invalid) { if (valid) free (valid); return n_invalid; } /* all strings are valid: return zero */

  /* now we can reduce the vector of strings */
  char_vector_reduce_to_valid_strings (vec, valid, n_valid);
  
  if (valid) free (valid);
  return n_invalid;
}

int
char_vector_remove_duplicate_strings (char_vector vec)
{
  int i, j, k, *valid, n_valid = 0;
  bool equal;

  valid = (int*) biomcmc_malloc (vec->nstrings * sizeof (int));

  if (!vec->next_avail) return 0; /* if char_vector doesn't have any elements yet */

  for (i=0; i < vec->nstrings - 1; i++) if (vec->string[i]) {
    valid[n_valid++] = i; /* valid since it is leftmost, and not NULLified before */
    for (j = i + 1; j < vec->nstrings; j++) if ((vec->string[j]) && (vec->nchars[i] == vec->nchars[j])) {
      for (equal = true, k = 0; (equal == true) && (k < (int) vec->nchars[i]); k++) /* scan both strings */
        if (vec->string[i][k] != vec->string[j][k]) equal = false; /* if they are different, do nothing */
      if (equal == true) { /* if they are equal, we free the memory and NULLify superfluous element */
        if (vec->string[j]) free (vec->string[j]); /* the condition is always true... */
        vec->string[j] = NULL;
      }
    }
  }
  if (vec->string[i]) valid[n_valid++] = i; /* check if last string is unique or not (my most common mistake) */

  if (n_valid == vec->nstrings) { if (valid) free (valid); return 0; }

  i = vec->nstrings;
  /* TODO: for alignments (seq AND seqname) we want both, the original and the reduced */
  char_vector_reduce_to_valid_strings (vec, valid, n_valid);

  if (valid) free (valid);
  return (i - n_valid); /* number of chopped elements (= (old vec->nstrings) - (new vec->nstrings) ) */
}

void
char_vector_reduce_to_valid_strings (char_vector vec, int *valid, int n_valid)
{
  int i;

  for (i=0; i < n_valid; i++) {
    vec->string[i] = vec->string[ valid[i] ];
    vec->nchars[i] = vec->nchars[ valid[i] ];
  }
  vec->nstrings = n_valid;

  vec->string = (char**)  biomcmc_realloc ((char**)  vec->string, n_valid * sizeof (char*));
  vec->nchars = (size_t*) biomcmc_realloc ((size_t*) vec->nchars, n_valid * sizeof (size_t));
}  

/* The hungarian method below is copied from http://www.informatik.uni-freiburg.de/~stachnis/misc.html
 * The (edited) original message follows:
 *
 ** libhungarian by Cyrill Stachniss, 2004  Solving the Minimum Assignment Problem using the 
 ** Hungarian Method.         ** This file may be freely copied and distributed! **
 **
 ** Parts of the used code was originally provided by the "Stanford GraphGase", but I made changes to this code.
 ** As asked by  the copyright node of the "Stanford GraphGase", I hereby proclaim that this file are *NOT* part of the
 ** "Stanford GraphGase" distrubition! */

void
hungarian_reset (hungarian p)
{
  int i, j;

  for (i = 0; i < p->size; i++) {
    p->col_mate[i] = p->unchosen_row[i] = p->row_dec[i] = p->slack_row[i] = p->row_mate[i] = p->parent_row[i] = p->col_inc[i] = p->slack[i] = 0;
    for (j = 0; j < p->size; j++) p->cost[i][j] = 0;
  }
  p->final_cost = 0;
}

hungarian
new_hungarian (int size)
{
  int i;
  hungarian p;

  p = (hungarian) biomcmc_malloc (sizeof (struct hungarian_struct)); 
  p->size = size; /* n_rows = n_columns; if it's not, fill with zeroes (no cost) */
  p->cost = (int**) biomcmc_malloc (size * sizeof (int*));
  for (i = 0; i < p->size; i++)
    p->cost[i] = (int*) biomcmc_malloc (size * sizeof (int));
  /* edges would be assignment_matrix[ i * ncols + col_mate[i] ] = true; and other elems "false" (but we don't use the matrix notation) */
  p->col_mate     = (int*) biomcmc_malloc (size * sizeof (int)); /* for a given row node, col_mate[row] is the assigned col node */
  p->unchosen_row = (int*) biomcmc_malloc (size * sizeof (int));
  p->row_dec      = (int*) biomcmc_malloc (size * sizeof (int));
  p->slack_row    = (int*) biomcmc_malloc (size * sizeof (int));
  p->row_mate     = (int*) biomcmc_malloc (size * sizeof (int));
  p->parent_row   = (int*) biomcmc_malloc (size * sizeof (int));
  p->col_inc      = (int*) biomcmc_malloc (size * sizeof (int));
  p->slack        = (int*) biomcmc_malloc (size * sizeof (int));

  hungarian_reset (p);
  return p;
}

void
hungarian_update_cost (hungarian p, int row, int col, int cost)
{
  if (row >= p->size) return;
  if (col >= p->size) return;
  p->cost[row][col] = cost;
}

void 
del_hungarian (hungarian p)
{
  int i;
  if (!p) return;
  if (p->cost) {
    for (i = p->size - 1; i >= 0; i--) if (p->cost[i]) free (p->cost[i]);
    free (p->cost);
  }
  free (p->col_mate); /* this is the important one, with i assigned to col_mate[i] */
  free (p->slack);
  free (p->col_inc);
  free (p->parent_row);
  free (p->row_mate);
  free (p->slack_row);
  free (p->row_dec);
  free (p->unchosen_row);
  free (p);
}

void 
hungarian_solve (hungarian p, int this_size)
{
  int i, j, nrows = this_size, ncols = this_size, k, l, s, t, q, unmatched;
  p->final_cost = p->initial_cost = 0;

  if (this_size > p->size) { p->final_cost = -1; return; } /* we don't call biomcmc_error(), but it *is* an error! */

  for (l = 0; l < ncols; l++) { // Begin subtract column minima in order to start with lots of zeroes 12
    s = p->cost[0][l];
    for (k = 1; k < nrows; k++) if (p->cost[k][l] < s) s = p->cost[k][l];
    p->initial_cost += s; /* this should be added to final_cost to have classical assignment cost; here we distinguish them */
    if (s!=0)	for (k = 0; k < nrows; k++) p->cost[k][l] -= s;
  } // End subtract column minima in order to start with lots of zeroes 12
 
  /*for (i=0;i<nrows; i++) {
    for (j=0;j<ncols; j++) printf ("%4d ", p->cost[i][j]);
    printf (" DEBUGcost\n");
  }  // DEBUG */ 

  // Begin initial state 16
  t=0;
  for (l = 0; l < ncols; l++)  { // n => num_cols
    p->row_mate[l]= -1;
    p->parent_row[l]= -1;
    p->col_inc[l]=0;
    p->slack[l]= 0x7FFFFFFF;
  }
  for (k = 0; k < nrows; k++) { // m => num_rows
    s = p->cost[k][0];
    for (l = 1; l < ncols; l++) if (p->cost[k][l] < s) s = p->cost[k][l];
    p->row_dec[k]=s;
    for (l = 0; l < ncols; l++) if ((s==p->cost[k][l]) && (p->row_mate[l] < 0)) {
      p->col_mate[k] = l;
      p->row_mate[l] = k;  // fprintf(stderr, "matching col %d==row %d\n",l,k);
      goto row_done;
    }
    p->col_mate[k] = -1;  // fprintf(stderr, "node %d: unmatched row %d\n",t,k);
    p->unchosen_row[t++] = k;
row_done:
    ;
  }
  // End initial state 16

  // Begin Hungarian algorithm 18
  if (t==0)    goto done;
  unmatched=t;
  while (1) {
    q=0; // fprintf(stderr, "Matched %d rows.\n",m-t);
    while (1)	{
      while (q<t) {
         { // Begin explore node q of the forest 19
          k = p->unchosen_row[q];
          s=p->row_dec[k];
          for (l=0;l<ncols;l++) if (p->slack[l]) {
            int del;
            del = p->cost[k][l] - s + p->col_inc[l];
            if (del < p->slack[l]) {
              if (del==0) {
                if (p->row_mate[l]<0)  goto breakthru;
                p->slack[l]=0;
                p->parent_row[l]=k; // fprintf(stderr, "node %d: row %d==col %d--row %d\n", t,row_mate[l],l,k);
                p->unchosen_row[t++]=p->row_mate[l];
              }
              else { p->slack[l]=del; p->slack_row[l]=k; }
            }
          }
         } // End explore node q of the forest 19
        q++;
      }

      // Begin introduce a new zero into the matrix 21
      s = 0x7FFFFFFF;
      for (l = 0;l < ncols; l++) if (p->slack[l] && p->slack[l] < s) s = p->slack[l];
      for (q = 0; q < t; q++) p->row_dec[ p->unchosen_row[q] ] += s;
      for (l = 0; l < ncols; l++) if (p->slack[l]) {
        p->slack[l]-=s;
        if (p->slack[l]==0) {  // Begin look at a new zero 22
          k = p->slack_row[l]; // fprintf(stderr, "Decreasing uncovered elements by %d produces zero at [%d,%d]\n", s,k,l);
          if (p->row_mate[l]<0)  {
            for (j=l+1;j<ncols;j++)  if (p->slack[j]==0) p->col_inc[j]+=s;
            goto breakthru;
          }
          else {
            p->parent_row[l]=k; // fprintf(stderr, "node %d: row %d==col %d--row %d\n",t,row_mate[l],l,k);
            p->unchosen_row[t++]=p->row_mate[l];
          }
        } // End look at a new zero 22

      }
      else  p->col_inc[l]+=s;
      // End introduce a new zero into the matrix 21
    }
breakthru:
    // fprintf(stderr, "Breakthrough at node %d of %d!\n",q,t);
    while (1)	{    // Begin update the matching 20
      j=p->col_mate[k];
      p->col_mate[k]=l;
      p->row_mate[l]=k; // fprintf(stderr, "rematching col %d==row %d\n",l,k);
      if (j<0)   break;
      k=p->parent_row[j];
      l=j;
    }    // End update the matching 20
    if (--unmatched==0)	goto done;
    // Begin get ready for another stage 17
    t=0;
    for (l=0;l<ncols;l++) {
      p->parent_row[l]= -1;
      p->slack[l]=0x7FFFFFFF;
    }
    for (k=0;k<nrows;k++) if (p->col_mate[k]<0) p->unchosen_row[t++]=k; // fprintf(stderr, "node %d: unmatched row %d\n",t,k);
    // End get ready for another stage 17
  }
done:
  // Begin doublecheck the solution 23
  for (k = 0; k < nrows; k++) for (l=0;l<ncols;l++) if (p->cost[k][l] < p->row_dec[k] - p->col_inc[l]) { p->final_cost = -1; printf ("\n**\n"); return;}
  for (k = 0; k < nrows; k++) {
    l=p->col_mate[k];
    if ((l < 0) || (p->cost[k][l] != p->row_dec[k] - p->col_inc[l])) { p->final_cost = -1; return; }
  }
  k=0;
  for (l=0;l<ncols;l++) if (p->col_inc[l])  k++;
  if (k>nrows) { p->final_cost = -1; return; }
  // End doublecheck the solution 23
  // End Hungarian algorithm 18

  for (k = 0; k < nrows; ++k) for (l = 0; l < ncols; ++l) p->cost[k][l] = p->cost[k][l] - p->row_dec[k] + p->col_inc[l];
  for (i = 0; i < nrows; i++) p->final_cost += p->row_dec[i];
  for (i = 0; i < ncols; i++) p->final_cost -= p->col_inc[i]; // fprintf(stderr, "Cost is %d\n",cost);
}
