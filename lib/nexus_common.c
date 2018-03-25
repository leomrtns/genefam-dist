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

#include "nexus_common.h"

/*! \brief Count balance of nexus comments (number of '['s minus number of ']'s). */
int   count_nested_nexus_comments (char *string);
/*! \brief Remove comments that start on this line (we need to guarantee that they also finish on this line). */
char* remove_oneline_nexus_comments (char *string);

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

char*
remove_nexus_comments (char **string, size_t *stringsize, FILE *stream)
{ /* shouldn't change address **string points to (like string++ or similar), since it is reallocated by readline */
  int nested_comment = count_nested_nexus_comments (*string);
  while (nested_comment > 0) {
    if (biomcmc_getline (string, stringsize, stream) == -1)  
      biomcmc_error ( "Premature end of file while removing nexus comments (unbalanced ['s and ]'s ?)\n");
    nested_comment += count_nested_nexus_comments (*string);
  }
  *string = remove_oneline_nexus_comments (*string);
  return (*string + strspn (*string, " \t")); /* remove leading spaces */
}

int 
count_nested_nexus_comments (char *string)
{
  char *s, *last = string + strlen (string);
  int count = 0;
  for (s = string; s <= last; s++) {
    if (*s == '[') count++;
    else if (*s == ']') count--;
  }
 return count;
} 

char* 
remove_oneline_nexus_comments (char *string) 
{
  char *s, *first, *last = string+strlen (string);
  int count = 1; /* we only use it if we find a "]" */
  for (s = last; (s >= string) && (*s != ']'); s--); /* find last "]", backwards */
  if (s < string) return (string); 

  first = s-1; /* s points to "]", first will find "[" now */ 
  s++;         /* first position after "]" */
  while ((count > 0) && (first >= string)) {
    if      ( *first == '[') count--;
    else if ( *first == ']') count++;
    first--;
  }
  *last = '\0'; memmove (first+1, s, last - s + 1);

  return remove_oneline_nexus_comments (string); /* recursive, to remove other comments before */
}

char* 
remove_oneline_nexus_comments_OLD (char *string) 
{
  char *s, *first, *last = string+strlen (string);
  int count = 1; /* we only use it if we find a "[" */
  for (s = string; (s <= last) && (*s != '['); s++);
  if (s>last) { 
    /*DEBUG printf ("{%s}", string);*/  /* FIXME: unbalanced "]" in multiline comments */
    return (string); 
  }
  first = s++;
  while ((count > 0) && (s <= last)) {
    if ( *s == '[') count++;
    else if (*s == ']') count--;
    s++;
  }
  if (s < last) { *last = '\0'; memmove (first, s, last - s + 1); }
  else return NULL;

  return string;
}

char*
lowercase_string (char *string)
{
  char *s, *last = string+strlen (string); 
  for (s = string; s <= last; s++) if (isupper (*s)) *s = tolower (*s);
  return string;
}

char*
uppercase_string (char *string)
{
  char *s, *last = string+strlen (string); 
  for (s = string; s <= last; s++) if (islower (*s)) *s = toupper (*s);
  return string;
}

char*
remove_space_from_string (char *string)
{
  char *last = strchr (string, '\0'), *s; 
  for (s = last-1; s >= string; s--) if (isspace (*s) || (*s == '>')) { memmove (s, s + 1, last - s); last--; }
  return string;
}

bool
nonempty_string (char *string)
{
  char *last = strchr (string, '\0'), *s;

  if (*(last-1) == '\n') *(--last) = '\0'; // chomp final "\n"

  if (!(last-string)) return false;
  for (s = string; s < last; s++) if (!isspace (*s)) return true;
  return false;
}

bool
nonempty_fasta_line (char *string)
{ /* nonempty and not a comment */
  char *last = strchr (string, '\0'), *s;

  if (*(last-1) == '\n') *(--last) = '\0'; // chomp final "\n"

  if ((last-string) < 2) return false; /* just one character may be skipped */
  for (s = string; (s < last) && isspace (*s); s++); /* walk until find a nonempty char (or EOL) */
  if (s == last) return false;
  if ((s<last) && ((*s == '#') || (*s == ';'))) return false; /*first char is comment delimiter ('#' is an extension)*/
  return true;
}

