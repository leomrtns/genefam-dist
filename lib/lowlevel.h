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

/*! \file lowlevel.h 
 *  \brief Lowest level header file. Header file for lowlevel.c
 */

#ifndef _biomcmc_lowlevel_h_
#define _biomcmc_lowlevel_h_

#include "config.h"

#include <stdio.h>      /* [ANSI C C89] */
#include <stdlib.h>     /* random number, searching, sorting, EXIT_SUCCESS [ANSI C C89] */
#include <string.h>     /* string manipulation [ANSI C C89] */
#include <stdarg.h>     /* Access to va_arg, va_list [ANSI C C89] */
#include <stdint.h>     /* standard integer types (int32_t typedef etc.) [C99]*/
#include <ctype.h>      /* char operation functions (e.g. isspace() ), case convertion [ANSI C C89] */
#include <math.h>       /* standard math functions (e.g. exp() ) [ANSI C C89] */
#include <float.h>      /* DBL_MAX_EXP, DBL_EPSILON constants (to avoid underflow etc) */
#include <time.h>       /* speed profiling(e.g. clock(), clock_gettime(), struct timespec ) [ANSI C C89] */
#include <unistd.h>     /* system values checking at runtime (e.g. sysconf() ) [POSIX C] */
#include <sys/time.h>   /* random seed (e.g. gettimeofday(), struct timeval) [POSIX C] */
#include <sys/times.h>  /* speed profiling in clock ticks (e.g. times() ) [POSIX.1 (or GNU extension?)] */ 
#include <sys/types.h>  /* pid_t for process ID, used together with unistd.h (may not be necessary) */
#include <unistd.h>     /* getpid() function, used together with sys/types.h (may not be necessary) */

#include <libgen.h> /* standard XPG basename() - the one provided by string.h is a GNU extension, fails on macOSX*/


#ifdef _OPENMP
#include <omp.h>         /* OpenMP parallel threading library when available */
#endif

/* Global constants */

#define EXP_1 2.71828182845904523536028747135266 /* Euler's number */

#define true  1U /*!< Boolean TRUE  */
#define false 0U /*!< Boolean FALSE */

#define MIN(x,y) (((x)<(y)) ? (x) : (y))
#define MAX(x,y) (((x)>(y)) ? (x) : (y))
#define MOD(a)   (((a)>0)   ? (a) :(-a))


/*! \brief Mnemonic for boolean (char is smaller than int) */
typedef unsigned char bool;

typedef struct char_vector_struct* char_vector;
typedef struct hungarian_struct* hungarian;

/*! \brief vector of strings (char vectors) of variable length */
struct char_vector_struct
{
  char **string;  /*! \brief vector of strings */
  int nstrings;   /*! \brief how many strings */
  size_t *nchars; /*! \brief length of allocated memory for each string excluding the ending '\0' (the actual size in 
                    use needs strlen() or a call to char_vector_compress() over the structure )*/
  int ref_counter;/*! \brief how many times this char_vector_struct is being used */
  int next_avail; /*! \brief next available position (empty string) */
};

struct hungarian_struct
{
  int **cost, *col_mate; /*! \brief cost matrix, and col_mate[row] with column match for row */
  int size,  /*! \brief assignment size. Cost is a square matrix, so size should be an overestimate where "missing" nodes are added w/ cost zero */
      initial_cost, /*! \brief sum of lowest input cost values for each column. The hungarian method rescales them so that minimum per column is zero */
      final_cost;   /*! \brief our final cost is on rescaled cost matrix, therefore to restore the "classical" optimal cost one should sum it with initial_cost */
  int *unchosen_row, *row_dec, *slack_row, *row_mate, *parent_row, *col_inc, *slack; /* aux vectors */
};

/*! \brief Memory-safe malloc() function.
 *
 *  Allocates size bytes and returns a pointer to the allocated memory. An error message is thrown in case of failure.
 *  \param[in] size allocated size, in bytes
 *  \return pointer to newly allocated memory */
void *biomcmc_malloc (size_t size);

/*! \brief Memory-safe realloc() function.
 *
 * Changes the size of the memory block pointed to by ptr to size bytes. An error message is thrown in case of failure.
 * \param[in] size allocated size, in bytes
 * \param[in,out] ptr pointer to previously allocated memory
 * \return pointer to newly allocated memory */
void *biomcmc_realloc (void *ptr, size_t size);

/*! \brief Memory-safe fopen() function.
 *
 * Opens the file whose name is the string pointed to by path and associates a stream with it. An error message is 
 * thrown in case of failure.
 * \param[in] path file name 
 * \param[in] mode opening mode ("r" for reading, "w" for writing, etc)
 * \result pointer to file stream */
FILE *biomcmc_fopen (const char *path, const char *mode);

/*! \brief Prints error message and quits program.
 *
 * similar to fprintf (stderr, ...), but exits after printing the message
 * \param[in] template va_list following same format as printf()
 * \result exits program */
void biomcmc_error (const char *template, ...);

/*! \brief Comparison between two integers used by qsort() in crescent order */
int compare_int_increasing (const void *a, const void *b);

/*! \brief Comparison between two integers used by qsort() in decrescent order */
int compare_int_decreasing (const void *a, const void *b);

/*! \brief Comparison between two double floats used by qsort() in crescent order */
int compare_double_increasing (const void *a, const void *b);

/*! \brief Comparison between two double floats used by qsort() in decrescent order */
int compare_double_decreasing (const void *a, const void *b);

/*! \brief read file line-by-line (like homonymous function from GNU C library)
 *
 * This implementation is originally from the CvsGui project (http://www.wincvs.org/). The explanation from the 
 * original file adapted to our system  follows:
 * \verbatim 
   Read up to (and including) a newline ("\n") from STREAM into *LINEPTR and null-terminate it. *LINEPTR is a pointer 
   returned from malloc (or NULL), pointing to *N characters of space.  It is realloc'd as necessary.  Return the 
   number of characters read (not including the null terminator), or -1 on error or EOF. \endverbatim */
int biomcmc_getline (char **lineptr, size_t *n, FILE *stream);


/*! \brief Create a vector of strings with initial size for each string of zero */
char_vector new_char_vector (int nstrings);

/*! \brief Create a vector of strings where each string is assigned an initial value of nchars */
char_vector new_char_vector_fixed_length (int nstrings, int nchars);

/*! \brief Delete vector of strings only after nobody is using it */
void del_char_vector (char_vector vec);

/*! \brief Link a previously allocated string (to avoid copying all characters) */
void char_vector_link_string_at_position (char_vector vec, char *string, int position);

/*! \brief Add a new string (vector of characters) at specific location */
void char_vector_add_string_at_position (char_vector vec, char *string, int position);

/*! \brief Add a new string (vector of characters) at next available location */
void char_vector_add_string (char_vector vec, char *string);

/*! \brief Append string at the end of existing string at location */
void char_vector_append_string_at_position (char_vector vec, char *string, int position);

/*! \brief Append string at the end of existing string at most recently used location */
void char_vector_append_string (char_vector vec, char *string);

/*! \brief Increase size of vector of strings (called automatically by other functions) */
void char_vector_expand_nstrings (char_vector vec, int new_size);

/*! \brief update order of strings in vector based on a vector of new positions */
void char_vector_reorder_strings (char_vector vec, int *order);

/*! \brief Reduce size of vector of strings by removing empty strings (returns number of empty strings) */
int char_vector_remove_empty_strings (char_vector vec);
/*! \brief Remove identical strings and resizes char_vector_struct */
int char_vector_remove_duplicate_strings (char_vector vec);
/*! \brief reduce char_string_struct to only those elements indexed by valid[] */
void char_vector_reduce_to_valid_strings (char_vector vec, int *valid, int n_valid);

/* Hungarian method for bipartite matching (assignment) */
hungarian new_hungarian (int size);
void hungarian_reset (hungarian p);
void hungarian_update_cost (hungarian p, int row, int col, int cost);
void del_hungarian (hungarian p);
void hungarian_solve (hungarian p, int this_size);

#endif
