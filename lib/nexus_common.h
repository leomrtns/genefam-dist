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

/*! \file nexus_common.h 
 *  \brief File handling functions for nexus format in general.
 */

#ifndef _biomcmc_nexus_common_h_
#define _biomcmc_nexus_common_h_

#include "lowlevel.h"

/*! \brief maximum name length for taxa (alignment and tree files). */
#define MAX_NAME_LENGTH 4096 

/*! \brief Removes (possible nested/multiline) nexus comments of the form [] (brackets). */
char* remove_nexus_comments (char **string, size_t *stringsize, FILE *stream);
/*! \brief Changes uppercase characters by lowercase versions. */
char* lowercase_string (char *string);
/*! \brief Changes lowercase characters by uppercase versions. */
char* uppercase_string (char *string);
/*! \brief Removes spaces, tabs from string. */
char* remove_space_from_string (char *string);
/*! \brief returns bool::false if string is composed only of space characters (' ', '\n', '\r', '\t', etc). */
bool nonempty_string (char *string);
/*! \brief returns bool::false if first nonspace character of string is ';' (FASTA comment) or '#' (HUPO extension) */
bool nonempty_fasta_line (char *string);

#endif
