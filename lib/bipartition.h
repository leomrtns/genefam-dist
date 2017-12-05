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

/*! \file bipartition.h 
 *  \brief Unary/binary operators on arbitrarily-sized bitstrings (strings of zeros and ones) like split bipartitions
 *  
 */

#ifndef _biomcmc_bipartition_h_
#define _biomcmc_bipartition_h_


#include "hashtable.h" // includes  "lowlevel.h"

typedef struct bipartition_struct* bipartition;
typedef struct bipsize_struct* bipsize;
typedef struct bip_hashtable_struct* bip_hashtable;
typedef struct bip_hashitem_struct* bip_hashitem;
typedef bipartition* tripartition; /* just a vector of size 3 */

/*! \brief Bit-string representation of splits. */
struct bipartition_struct 
{
  /*! \brief Representation of a bipartition by a vector of integers (bitstrings). */
  uint64_t *bs;
  /*! \brief Counter (number of "one"s) */
  int n_ones;
  /*! \brief number of bits (leaves), vector size and mask */
  bipsize n;
  /*! \brief How many times this struct is being referenced */
  int ref_counter;
};

struct bipsize_struct
{
  /*! \brief mask to make sure we consider only active positions (of last bitstring) */
  uint64_t mask;
  /*! \brief Vector size and total number of elements (n_ints = n_bits/(8*sizeof(long long)) +1). */
  int ints, bits, original_size;
  /*! \brief How many times this struct is being referenced */
  int ref_counter;
};

/*! \brief Hash table of bipartitions (see hashtable.h for original version, with string keys and integer values) */
struct bip_hashtable_struct 
{ 
  int size; /*! \brief Table size. */
  int probelength;  /*! \brief Number of collisions before empty slot is found. */
  int maxfreq; /*! \brief frequency (integer) of most frequent bipartition */
  uint32_t h; /*! \brief Value set by hash(). Used in hash1() and hash2() to avoid calling hash() again. */
  uint32_t a1, a2, b1, b2, P; /*!< \brief Random values used in hash functions. */
  bip_hashitem* table; /*! \brief Vector with key/value pairs. */
  int ref_counter;  /*! \brief Counter of how many external references (structures sharing this hashtable) to avoid deletion */
};

/*! \brief key (bipartition) and value (frequency) pair for hash table of bipartitions */
struct bip_hashitem_struct {
  bipartition key; /*! \brief pointer to bipartition (must update ref_counter) */
  int count;  /*! \brief frequency of bipartition (counter, but can be scaled to max count so far) */
};

/*! \brief create a new bipartition (bitstring) capable of storing an arbitrary number of bits and initialize it to zero
 * \param[in] size number of bits of desired bipartition
 * \returns bipartition (opaquely a vector of long long ints) */
bipartition new_bipartition (int size);
/*! \brief Create a new bipsize, which controls some bipartition sizes */
bipsize new_bipsize (int size);
/*! \brief Create a new bipartition (allocate memory) and initialize it from another bipartition */
bipartition new_bipartition_copy_from (const bipartition from);
/*! \brief create new bipartition that will share bipsize -- useful for bipartition vectors */
bipartition new_bipartition_from_bipsize (bipsize n);
/*! \brief free memory allocated by bipartition */
void del_bipartition (bipartition bip);
/*! \brief free memory allocated by bipsize */
void del_bipsize (bipsize n);
/*! \brief update the valid number of bits and mask -- e.g. when replacing subtrees by leaves in reduced trees */
void bipsize_resize (bipsize n, int nbits);
/*! \brief set all bits to zero except the one at position-th bit */
void bipartition_initialize (bipartition bip, int position);
/*! \brief set all bits to zero  */
void bipartition_zero (bipartition bip);
/*! \brief simply set the bit at "position" to one, irrespective of other bits */
void bipartition_set (bipartition bip, int position);
void bipartition_set_lowlevel (bipartition bip, int i, int j);
/*! \brief simply unset the bit at "position" (set to zero), irrespective of other bits */
void bipartition_unset (bipartition bip, int position);
void bipartition_unset_lowlevel (bipartition bip, int i, int j);
/*! \brief Copy contents from one bipartition to another */
void bipartition_copy (bipartition to, const bipartition from);
/*! \brief Binary logical OR ("|") between b1 and b2, where update_count should be true if you need to know the
 * resulting size (slow) or false if you don't care or if b1 and b2 are disjoint (no common elements) */
void bipartition_OR (bipartition result, const bipartition b1, const bipartition b2, bool update_count);
/*! \brief Binary logical AND ("&") between b1 and b2, update_count should be set to false  only if you <b>really</b>
 * don't need to know the number of active bits (e.g. sorting, bipartition comparison) */
void bipartition_AND (bipartition result, const bipartition b1, const bipartition b2, bool update_count);
/*! \brief Binary logical AND ("&") between b1 and ~b2 (NOT b2), that is, apply mask b1 on the inverse of b2 */
void bipartition_ANDNOT (bipartition result, const bipartition b1, const bipartition b2, bool update_count);
/*! \brief Binary logical eXclusive OR ("^") between b1 and b2, update_count should be set to false  only if you <b>really</b>
 * don't need to know the number of active bits (e.g. sorting, bipartition comparison) */
void bipartition_XOR (bipartition result, const bipartition b1, const bipartition b2, bool update_count);
/*! \brief Binary logical eXclusive OR ("^") between b1 and complement of b2 (that is, NOT b2: b1 ^ ~b2). Used when
 *finding best disagreement (that in this case erases the complement --other side -- of agreement edge) */
void bipartition_XORNOT (bipartition result, const bipartition b1, const bipartition b2, bool update_count);
/*! \brief Unary complement ("~") of bipartition. Use with caution, since there is no mask for unused padded bits */ 
void bipartition_NOT (bipartition result, const bipartition bip);
/*! \brief Count the number of active bits (equal to one). Used by bipartition_AND() and bipartition_XOR() when 
 * update_count = true. Please use it parsimoniously since it is as slow as without bitstring representation. */
void bipartition_count_n_ones (const bipartition bip);
/*! \brief fill vector id[] with positions of set bits, up to vecsize bits set */
void bipartition_to_int_vector (const bipartition b, int *id, int vecsize);
/*! \brief Compare equality of two bipartitions */
bool bipartition_is_equal (const bipartition b1, const bipartition b2);
/*! \brief Compare if two bipartitions represent the same splits (or they are equal or one is the complement of the other) */
bool bipartition_is_equal_bothsides (const bipartition b1, const bipartition b2);
/*! \brief Compare sizes of two bipartitions, by number of active bits with ties broken by actual bitstrings */
bool bipartition_is_larger (const bipartition b1, const bipartition b2);
/*! \brief invert ones and zeroes in loco when necessary to assure bipartition has more zeroes than ones */
void bipartition_flip_to_smaller_set (bipartition bip);
/*! \brief Check if position-th bit is equal to one or not */
bool bipartition_is_bit_set (const bipartition bip, int position);
/*! \brief Check if first bipartition contains all elements of second bipartition (b2 is a subset of b1) */
bool bipartition_contains_bits (const bipartition b1, const bipartition b2);
/*! \brief Print to screen a bit representation of the bipartition (with number of ones at the end) */
void bipartition_print_to_stdout (const bipartition b1);
/*! \brief replace bit info, copying 'from' one position 'to' another; bool "update" indicates if afterwards size will be reduced */
void bipartition_replace_bit_in_vector (bipartition *bvec, int n_b, int to, int from, bool reduce);
/*! \brief apply mask to last element (useful after manipulations) and count number of bits */
void bipartition_resize_vector (bipartition *bvec, int n_b);

/*! \brief 32bits hash value for bipartition */
uint32_t bipartition_hash (bipartition bip);
/*! \brief Create new hashtable of size bipartitions. */
bip_hashtable new_bip_hashtable (int size);
/*! \brief Free bipartition hashtable space. */
void  del_bip_hashtable (bip_hashtable ht);
/*! \brief Insert key (bipartition) into bipartition hashtable, adding one to its count (freq). */
void bip_hashtable_insert (bip_hashtable ht, bipartition key);
/*! \brief Return frequency of bipartition (count/maxfreq) or zero if not found. */
double bip_hashtable_get_frequency (bip_hashtable ht, bipartition key);

/*! \brief tripartition of a node (a vector with 3 bipartitions, that should not be 'flipped' to smaller set, however) */
tripartition new_tripartition (int nleaves);
/*! \brief free tripartition space (just 3 bipartitions) */
void del_tripartition (tripartition trip);
/*! \brief from node, create tripartition from node->left and node->right (assuming bipartitions were not 'flipped' yet) */
void store_tripartition_from_bipartitions (tripartition tri, bipartition b1, bipartition b2);
/*! \brief sort order of bipartitions s.t. smallest is first */
void sort_tripartition (tripartition tri);
/*! \brief  match bipartitions between two nodes and return optimal score (min disagreement) */
int align_tripartitions (tripartition tp1, tripartition tp2, hungarian h);
/*! \brief  assuming tripartitions are ordered, check if nodes (represented by tripartitions) are the same */
bool tripartition_is_equal (tripartition tp1, tripartition tp2);


#endif
