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

/*! \file empirical_frequency.h 
 *  \brief Creates a histogram of a vector, ordered by frequency. 
 *
 * Sorts a vector of integers by their frequencies, preserving their original indexes. It is a simple extension to qsort
 * where the original order can be reconstructed, or still a key/value sorting.
 */

#ifndef _biomcmc_empirical_frequency_h_
#define _biomcmc_empirical_frequency_h_

#include "lowlevel.h"

typedef struct empfreq_struct* empfreq;
typedef struct empfreq_double_struct* empfreq_double;

typedef struct
{
	int freq;
	int idx;
} empfreq_element;

struct empfreq_struct
{
	empfreq_element *i;
	int n;
  /*! \brief Min value for index. */
  int min; 
  /*! \brief Max value for index. */
  int max;
};

struct empfreq_double_struct
{
  /*! \brief normalized frequencies in log scale, s.t. sum up to one (after exponentiating...). */
  double *freq;
  /*! \brief normalized idx values (used in ndups, nloss) */
  double *idx;
	int n;
  /*! \brief Min value for index. */
  double min; 
  /*! \brief Max value for index. */
  double max;
};

void sort_empfreq_decreasing (empfreq ef);
void sort_empfreq_increasing (empfreq ef);

empfreq new_empfreq (int n_elements);
void del_empfreq (empfreq ef);

empfreq_double new_empfreq_double (empfreq ef); /* works only with existing empfreq; limited functionality... */
void del_empfreq_double (empfreq_double efd);

empfreq new_empfreq_sort_decreasing (void *vec, int n, char type);
empfreq new_empfreq_sort_increasing (void *vec, int n, char type);
empfreq new_empfreq_from_int (int *vec, int n);
empfreq new_empfreq_from_int_weighted (int *vec, int n, int *weight);
int find_mode_int (int *vec, int n);
int find_mode_int_weighted (int *vec, int n, int *weight);


#endif
