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

/*! \file empirical_frequency.c  
 *  \brief histogram of vectors, ordered by frequency. Also calculates MAP (modal) values.
 *
 * Sorts a vector of integers by their frequencies, preserving their original indexes. It is a simple extension to qsort
 * where the original order can be reconstructed, or still a key/value sorting.
 */

#include "empirical_frequency.h"

int compare_empfreq_element_decreasing (const void *a, const void *b);
int compare_empfreq_element_increasing (const void *a, const void *b);
int compare_empfreq_double_element_decreasing (const void *a, const void *b);
int compare_empfreq_double_element_increasing (const void *a, const void *b);
empfreq create_empfreq_from_value_sorted_empfreq (empfreq e_idx);

int
compare_empfreq_element_decreasing (const void *a, const void *b)
{
  int result = ((empfreq_element *)b)->freq - ((empfreq_element *)a)->freq;
  if (result) return result;
  /* break ties by position (in decreasing order) */
  return ((empfreq_element *)b)->idx - ((empfreq_element *)a)->idx;
}

int
compare_empfreq_element_increasing (const void *a, const void *b)
{
  int result = ((empfreq_element *)a)->freq - ((empfreq_element *)b)->freq;
  if (result) return result;
  /* break ties by position (in increasing order) */
  return ((empfreq_element *)a)->idx - ((empfreq_element *)b)->idx;
}

int
compare_empfreq_double_element_decreasing (const void *a, const void *b)
{
  if ( ((empfreq_double_element *)b)->freq > ((empfreq_double_element *)a)->freq ) return 1;
  else return -1;
}

int
compare_empfreq_double_element_increasing (const void *a, const void *b)
{
  if ( ((empfreq_double_element *)b)->freq < ((empfreq_double_element *)a)->freq ) return 1;
  else return -1;
}

void
sort_empfreq_decreasing (empfreq ef)
{
  qsort (ef->i, ef->n, sizeof (empfreq_element), compare_empfreq_element_decreasing);
}

void
sort_empfreq_increasing (empfreq ef)
{
  qsort (ef->i, ef->n, sizeof (empfreq_element), compare_empfreq_element_increasing);
}

void
sort_empfreq_double_decreasing (empfreq_double efd)
{
  qsort (efd->d, efd->n, sizeof (empfreq_double_element), compare_empfreq_double_element_decreasing);
}

void
sort_empfreq_double_increasing (empfreq_double efd)
{
  qsort (efd->d, efd->n, sizeof (empfreq_double_element), compare_empfreq_double_element_increasing);
}

empfreq
new_empfreq (int n_elements)
{
  empfreq ef;
  int i;
  ef =  (empfreq) biomcmc_malloc (sizeof (struct empfreq_struct));
  ef->n = n_elements;
  ef->min = 0;
  ef->max = n_elements - 1;
  ef->i = (empfreq_element*) biomcmc_malloc (n_elements * sizeof (empfreq_element));
  for (i=0; i< n_elements; i++) { ef->i[i].freq = 0; ef->i[i].idx  = i; }
  return ef;
}

void
del_empfreq (empfreq ef)
{
  if (!ef) return;
  if (ef->i) free (ef->i);
  free (ef);
}

empfreq_double
new_empfreq_double (int n_elements)
{
  empfreq_double efd;
  int i;
  efd =  (empfreq_double) biomcmc_malloc (sizeof (struct empfreq_double_struct));
  efd->n = n_elements;
  efd->min = 0.;
  efd->max = (double)(n_elements) - 1.;
  efd->d = (empfreq_double_element*) biomcmc_malloc (n_elements * sizeof (empfreq_double_element));
  for (i=0; i< n_elements; i++) { efd->d[i].freq = 0.; efd->d[i].idx  = i; }
  return efd;
}

void
del_empfreq_double (empfreq_double efd)
{
  if (!efd) return;
  if (efd->d) free (efd->d);
  free (efd);
}

empfreq
new_empfreq_sort_decreasing (void *vec, int n, char type)
{
  int i;
  empfreq e_idx = new_empfreq (n);

  /* handles char (weird but possible), size_t (order char_vector[] by seq length) and int (type=0,1,2) */
  if (type == 0) for (i=0; i < n; i++) e_idx->i[i].freq = (int) ((char*)vec)[i];
  if (type == 1) for (i=0; i < n; i++) e_idx->i[i].freq = (int) ((size_t*)vec)[i]; 
  if (type == 2) for (i=0; i < n; i++) e_idx->i[i].freq = (int) ((int*)vec)[i]; 

  sort_empfreq_decreasing (e_idx); /* equiv. to qsort (vec, n, sizeof (int), compare_int) but preserving weights. */
  return e_idx;
}

empfreq
new_empfreq_sort_increasing (void *vec, int n, char type)
{
  int i;
  empfreq e_idx = new_empfreq (n);

  if (type == 0) for (i=0; i < n; i++) e_idx->i[i].freq = (int) ((char*)vec)[i];
  if (type == 1) for (i=0; i < n; i++) e_idx->i[i].freq = (int) ((size_t*)vec)[i]; 
  if (type == 2) for (i=0; i < n; i++) e_idx->i[i].freq = (int) ((int*)vec)[i]; 

  sort_empfreq_increasing (e_idx); /* equiv. to qsort (vec, n, sizeof (int), compare_int) but preserving weights. */
  return e_idx;
}

empfreq_double
new_empfreq_double_sort_decreasing (double *vec, int n)
{
  int i;
  empfreq_double e_idx = new_empfreq_double (n);
  for (i=0; i < n; i++) e_idx->d[i].freq = vec[i];
  sort_empfreq_double_decreasing (e_idx); /* equiv. to qsort (vec, n, sizeof (double), compare_double) but preserving weights. */
  return e_idx;
}

empfreq_double
new_empfreq_double_sort_increasing (double *vec, int n)
{
  int i;
  empfreq_double e_idx = new_empfreq_double (n);
  for (i=0; i < n; i++) e_idx->d[i].freq = vec[i];
  sort_empfreq_double_increasing (e_idx); /* equiv. to qsort (vec, n, sizeof (double), compare_double) but preserving weights. */
  return e_idx;
}

empfreq
new_empfreq_from_int (int *vec, int n)
{
  int i;
  empfreq e_idx, e_count;

  e_idx = new_empfreq (n);
  for (i=0; i < n; i++) {
    /* trick to make it sort by value ( vec[] ) and not by frequency (the "real" sort comes below) */
    e_idx->i[i].idx  = 1; 
    e_idx->i[i].freq = vec[i];
  }

  sort_empfreq_increasing (e_idx); /* equiv. to qsort (vec, n, sizeof (int), compare_int) but preserving weights. */
  e_count = create_empfreq_from_value_sorted_empfreq (e_idx); /* create e_count with freq = e_idx->idx */

  del_empfreq (e_idx);
  return e_count;
}

empfreq
new_empfreq_from_int_weighted (int *vec, int n, int *weight)
{
  int i, new_n = 0, *nonzero;
  empfreq e_idx, e_count;

  nonzero = (int*) biomcmc_malloc (n * sizeof (int));

  /* excludes elements of vec with weight=0 */
  for (i=0; i < n; i++) if (weight[i]) { nonzero[new_n++] = i; }

  if (!new_n) biomcmc_error ("vector of empirical frequencies with all freqs=0");

  e_idx = new_empfreq (new_n);
  for (i=0; i < new_n; i++) {
    /* trick to make it sort by value and not by frequency (the "real" sort comes below) */
    e_idx->i[i].idx  = weight[ nonzero[i] ];
    e_idx->i[i].freq =    vec[ nonzero[i] ];
  }

  sort_empfreq_increasing (e_idx); /* equiv. to qsort (vec, n, sizeof (int), compare_int) but preserving weights. */
  e_count = create_empfreq_from_value_sorted_empfreq (e_idx); /* create e_count with freq = e_idx->idx */

  del_empfreq (e_idx);
  if (nonzero) free (nonzero);
  return e_count;
}

empfreq
create_empfreq_from_value_sorted_empfreq (empfreq e_idx)
{
  int i, distinct_values = 1;
  empfreq e_count;

  for (i=1; i < e_idx->n; i++) if (e_idx->i[i].freq != e_idx->i[i-1].freq) distinct_values++;

  e_count = new_empfreq (distinct_values);

  for (distinct_values = 0, i = 0; i < e_idx->n - 1; i++) {
    e_count->i[distinct_values].freq += e_idx->i[i].idx;
    e_count->i[distinct_values].idx = e_idx->i[i].freq;
    if (e_idx->i[i].freq != e_idx->i[i+1].freq) distinct_values++;
  }
  e_count->i[distinct_values].freq += e_idx->i[i].idx;
  e_count->i[distinct_values].idx = e_idx->i[i].freq;

  e_count->min = e_count->i[0].idx;
  e_count->max = e_count->i[distinct_values].idx;

  sort_empfreq_decreasing (e_count); /* Now it will sort from largest to smallest freq (sum of weights) */
  return e_count;
}

int
find_mode_int (int *vec, int n)
{
  empfreq e_count;
  int map;
  e_count = new_empfreq_from_int (vec, n);
  map = e_count->i[0].idx;
  del_empfreq (e_count);
  return map;
}

int
find_mode_int_weighted (int *vec, int n, int *weight)
{
  empfreq e_count;
  int map;
  e_count = new_empfreq_from_int_weighted (vec, n, weight);
  map = e_count->i[0].idx;
  del_empfreq (e_count);
  return map;
}
