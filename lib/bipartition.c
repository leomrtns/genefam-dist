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

#include "bipartition.h"

int BitStringSize = 8 * sizeof (uint64_t);

bipartition
new_bipartition (int size)
{
  bipartition bip;
  int i;

  bip = (bipartition) biomcmc_malloc (sizeof (struct bipartition_struct));
  bip->n = new_bipsize (size);
  bip->n_ones = 0;
  bip->ref_counter = 1;

  bip->bs = (uint64_t*) biomcmc_malloc (bip->n->ints * sizeof (uint64_t));
  for (i=0; i < bip->n->ints; i++) bip->bs[i] = 0LL;

  return bip;
}

bipsize
new_bipsize (int size)
{
  bipsize n;
  int i;

  n = (bipsize) biomcmc_malloc (sizeof (struct bipsize_struct));
  n->bits = n->original_size = size;
  n->ref_counter = 1;
  n->ints = size/BitStringSize + 1;
  n->mask = 0LL;
  for (i=0; i < n->bits%BitStringSize; i++) n->mask |= (1LL << i); /* disregard other bits */

  return n;
}

bipartition
new_bipartition_copy_from (const bipartition from)
{
  bipartition bip;
  int i;

  bip = (bipartition) biomcmc_malloc (sizeof (struct bipartition_struct));
  bip->n = new_bipsize (from->n->bits);
  bip->n_ones = from->n_ones;
  bip->ref_counter = 1;

  bip->bs = (uint64_t*) biomcmc_malloc (bip->n->ints * sizeof (uint64_t));
  for (i=0; i < bip->n->ints; i++) bip->bs[i] = from->bs[i];

  return bip;
}

bipartition
new_bipartition_from_bipsize (bipsize n)
{
  bipartition bip;
  int i;

  bip = (bipartition) biomcmc_malloc (sizeof (struct bipartition_struct));
  bip->n = n;
  bip->n->ref_counter++;
  bip->n_ones = 0;
  bip->ref_counter = 1;

  bip->bs = (uint64_t*) biomcmc_malloc (bip->n->ints * sizeof (uint64_t));
  for (i=0; i < bip->n->ints; i++) bip->bs[i] = 0LL;

  return bip;
}

void
del_bipartition (bipartition bip)
{
  if (bip) {
    if (--bip->ref_counter) return;
    if (bip->bs) free (bip->bs); 
    del_bipsize (bip->n);
    free (bip); 
  }
}

void
del_bipsize (bipsize n)
{
  if (n) {
    if (--n->ref_counter) return;
    free (n);
  }
}

void
bipsize_resize (bipsize n, int nbits)
{
  int i;
  n->bits = nbits;
  n->ints = nbits/BitStringSize + 1; // might be smaller than original bs size 
  n->mask = 0LL;
  for (i=0; i < nbits%BitStringSize; i++) n->mask |= (1LL << i); /* disregard other bits */
}

void
bipartition_initialize (bipartition bip, int position)
{
  int i, j;
  for (i=0; i < bip->n->ints; i++) bip->bs[i] = 0LL;
  j = position%BitStringSize; 
  i = position/BitStringSize;

  bip->bs[i] = (1LL << j);
  bip->n_ones = 1;
}

void
bipartition_zero (bipartition bip)
{
  int i;
  for (i=0; i < bip->n->ints; i++) bip->bs[i] = 0LL;
  bip->n_ones = 0;
}

void
bipartition_set (bipartition bip, int position)
{
  bipartition_set_lowlevel (bip, position/BitStringSize, position%BitStringSize);
}

void
bipartition_set_lowlevel (bipartition bip, int i, int j)
{
  if (bip->bs[i] & (1LL << j)) return; // bit already set
  bip->bs[i] |= (1LL << j);
  bip->n_ones++; /* doesn't work if we reduce space later (check replace_int_in_vector() ) */
}

void
bipartition_unset (bipartition bip, int position)
{
  bipartition_unset_lowlevel (bip, position/BitStringSize, position%BitStringSize);
}

void
bipartition_unset_lowlevel (bipartition bip, int i, int j)
{
  if (!(bip->bs[i] & (1LL << j))) return; // bit already unset
  bip->bs[i] &= ~(1LL << j);
  bip->n_ones--;
}

void
bipartition_copy (bipartition to, const bipartition from)
{
  int i;
  for (i=0; i < to->n->ints; i++) to->bs[i] = from->bs[i];
  to->n_ones = from->n_ones;
}

void
bipartition_OR (bipartition result, const bipartition b1, const bipartition b2, bool update_count)
{
  int i;
  for (i=0; i < result->n->ints; i++) result->bs[i] = b1->bs[i] | b2->bs[i];
  result->bs[i-1] &= b1->n->mask; /* do not change last bits (do not belong to bipartition) */
  if (update_count) bipartition_count_n_ones (result);
  else result->n_ones = b1->n_ones + b2->n_ones; // works on topologies where b1 and b2 are disjoint 
}

void
bipartition_AND (bipartition result, const bipartition b1, const bipartition b2, bool update_count)
{
  int i;
  for (i=0; i < result->n->ints; i++) result->bs[i] = b1->bs[i] & b2->bs[i];
  result->bs[i-1] &= b1->n->mask; /* do not change last bits (do not belong to bipartition) */
  if (update_count) bipartition_count_n_ones (result);
  else result->n_ones = 0;// update_count = false should be used only when you don't care about this value (temp var)
}

void
bipartition_ANDNOT (bipartition result, const bipartition b1, const bipartition b2, bool update_count)
{
  int i;
  for (i=0; i < result->n->ints; i++) result->bs[i] = b1->bs[i] & (~b2->bs[i]);
  result->bs[i-1] &= b1->n->mask; /* do not change last bits (do not belong to bipartition) */
  if (update_count) bipartition_count_n_ones (result);
  else result->n_ones = 0;// update_count = false should be used only when you don't care about this value (temp var)
}

void
bipartition_XOR (bipartition result, const bipartition b1, const bipartition b2, bool update_count)
{
  int i;
  for (i=0; i < result->n->ints; i++) result->bs[i] = b1->bs[i] ^ b2->bs[i];
  result->bs[i-1] &= b1->n->mask; /* do not change last bits (do not belong to bipartition) */
  if (update_count) bipartition_count_n_ones (result);
  else result->n_ones = 0;// update_count = false should be used only when you don't care about this value (temp var)
}

void
bipartition_XORNOT (bipartition result, const bipartition b1, const bipartition b2, bool update_count)
{ /* equivalent to XOR followed by NOT */
  int i;
  for (i=0; i < result->n->ints; i++) result->bs[i] = b1->bs[i] ^ (~b2->bs[i]);
  result->bs[i-1] &= b1->n->mask; /* do not change last bits (do not belong to bipartition) */
  if (update_count) bipartition_count_n_ones (result);
  else result->n_ones = 0;// update_count = false should be used only when you don't care about this value (temp var)
}

void
bipartition_NOT (bipartition result, const bipartition bip)
{
  int i;
  for (i=0; i < result->n->ints; i++) result->bs[i] = ~bip->bs[i];
  result->bs[i-1] &= bip->n->mask; /* do not invert last bits (do not belong to bipartition) */
  result->n_ones = bip->n->bits - bip->n_ones;
}

void
bipartition_count_n_ones (const bipartition bip)
{
  int i;
  uint64_t j;
  bip->n_ones = 0;
/* // Naive approach 
  for (i=0; i < bip->n_ints - 1; i++) for (j=0; j < BitStringSize; j++) bip->n_ones += ((bip->bs[i] >> j) & 1LL);
  for (j=0; j < bip->n_bits%BitStringSize; j++) bip->n_ones += ((bip->bs[i] >> j) & 1LL);
 */
  // clear the least significant bit set per iteration (Peter Wegner in CACM 3 (1960), 322, mentioned in K&R)
  for (i=0; i < bip->n->ints; i++) for (j = bip->bs[i]; j; bip->n_ones++) j &= j - 1LL;
}

bool
bipartition_is_equal (const bipartition b1, const bipartition b2)
{
  int i;
  if (b1->n_ones  != b2->n_ones)  return false;
  if (b1->n->ints != b2->n->ints) return false;
  for (i=0; i < b1->n->ints - 1; i++) if (b1->bs[i] != b2->bs[i]) return false;
  b1->bs[i] &= b1->n->mask; b2->bs[i] &= b2->n->mask; /* apply mask before comparing last elems */
  if (b1->bs[i] != b2->bs[i]) return false; 
  return true;
}

bool
bipartition_is_equal_bothsides (const bipartition b1, const bipartition b2)
{
  int i;
  bool equal = true;
  for (i=0; (i < b1->n->ints - 1) && (equal); i++) if (b1->bs[i] != b2->bs[i]) equal = false;
  if ((equal) && ((b1->bs[i] & b1->n->mask) != (b2->bs[i] & b2->n->mask))) equal = false; 
  if (equal) return true; /* the biparitions are already the same, without flipping the bits */
  /* now we compare one bipartition with the complement of the other */
  for (i=0; (i < b1->n->ints - 1); i++) if (b1->bs[i] != ~b2->bs[i]) return false;
  if ((b1->bs[i] & b1->n->mask) != ((~b2->bs[i]) & b2->n->mask)) return false;
  return true; /* they are the exact complement of one another */
}

bool
bipartition_is_larger (const bipartition b1, const bipartition b2)
{
  int i;
  if (b1->n_ones > b2->n_ones) return true;
  if (b1->n_ones < b2->n_ones) return false;

  for (i = b1->n->ints - 1; (i >= 0) && (b1->bs[i] == b2->bs[i]); i--); /* find position of distinct bipartition elem*/

  if (i < 0) return false; /* identical bipartitions */
  if (b1->bs[i] > b2->bs[i]) return true;
  else return false;
}

void
bipartition_flip_to_smaller_set (bipartition bip)
{
  int i = bip->n->ints - 1; /* most significant position -- consistent with is_larger() above, using OLD algo below */
  if ((2 * bip->n_ones) < bip->n->bits) return; /* it is already the smaller set */
  /* OLD always x is different from ~x, so we just look at last element ("largest digits of number") */
  // if (((2 * bip->n_ones) == bip->n->bits) && (bip->bs[i] < (bip->n->mask & ~bip->bs[i]))) return;
  /* NEW: resolve ties by always showing the same "side" of bipartition, that is, the one having an arbitrary leaf (first one, in our case) */
  if (((2 * bip->n_ones) == bip->n->bits) && (bip->bs[0] & 1LL)) return;

  for (i=0; i < bip->n->ints; i++) bip->bs[i] = ~bip->bs[i]; /* like bipartition_NOT() */
  bip->bs[i-1] &= bip->n->mask; /* do not invert last bits (do not belong to bipartition) */
  bip->n_ones = bip->n->bits - bip->n_ones;
  return;
}

bool
bipartition_is_bit_set (const bipartition bip, int position)
{
  if (bip->bs[(int)(position/BitStringSize)] & (1LL << (int)(position%BitStringSize))) return true; 
  return false;
}

bool
bipartition_contains_bits (const bipartition b1, const bipartition b2)
{ /* generalization of bipartition_is_bit_set(); b1 contains or not b2 */
  int i;
  if (b1->n_ones < b2->n_ones) return false;
  for (i=0; i < b1->n->ints; i++) if ((b2->bs[i]) && (b2->bs[i] != (b1->bs[i] & b2->bs[i]))) return false;
  return true;
}

void
bipartition_to_int_vector (const bipartition b, int *id, int vecsize)
{
  int i, j, k = 0;
  for (i=0; i < b->n->ints; i++) for (j=0; (j < BitStringSize) && (k < vecsize); j++) if ( ((b->bs[i] >> j) & 1LL) ) id[k++] = i * BitStringSize + j;
}

void
bipartition_print_to_stdout (const bipartition b1)
{
  int i, j;
  for (i = 0; i < b1->n->ints - 1; i++) {
    for (j = 0; j < BitStringSize; j++) printf ("%d", (int)((b1->bs[i] >> j) & 1LL));
    printf (".");
  }
  for (j = 0; j < b1->n->bits%BitStringSize; j++) printf ("%d", (int)((b1->bs[i] >> j) & 1LL));
  printf ("[%d] ", b1->n_ones);
}

void
bipartition_replace_bit_in_vector (bipartition *bvec, int n_b, int to, int from, bool reduce)
{ /* copy info from position "from" to position "to" */
  int k, j = from%BitStringSize, i = from/BitStringSize, j2 = to%BitStringSize, i2 = to/BitStringSize;

  /* boolean "reduce" means that bitstring space will be reduced (last bits will be removed), therefore the update of
   * n_ones is different from default bipartition_set() behaviour: it's not an extra "1" (that is, one that did not
   * contribute to n_ones), but an existing "1" that change places. Schematically:
   * from -> to | normal n_ones count | when bitstring is reduced afterwards
   * 0    -> 0  |     0               |  0 
   * 0    -> 1  |    -1               | -1
   * 1    -> 0  |    +1               |  0  (since it's a leaf that belonged to position "from" and now is on position "to")
   * 1    -> 1  |     0               | -1  (in fact one of the two "1"s dissapeared after reducing the bitstring 
   * (the above description is outdated since I rewrote by hand the bit functions -- observe how we must erase 1 values from "from") */
  if (reduce) for (k = 0; k < n_b; k++) { // copy 0 or 1 values, erasing "from" values to avoid problems after reducing space (hanging 1s out of range)
    if      ( ((bvec[k]->bs[i] >> j) & 1LL) && ((bvec[k]->bs[i2] >> j2) & 1LL) )  { bvec[k]->n_ones--; bvec[k]->bs[i] &= ~(1LL << j); }
    else if ( ((bvec[k]->bs[i] >> j) & 1LL) && !((bvec[k]->bs[i2] >> j2) & 1LL) ) { bvec[k]->bs[i2] |=  (1LL << j2); bvec[k]->bs[i] &= ~(1LL << j); }
    else if ( !((bvec[k]->bs[i] >> j) & 1LL) && ((bvec[k]->bs[i2] >> j2) & 1LL) ) { bvec[k]->bs[i2] &= ~(1LL << j2); bvec[k]->n_ones--; }
    /* else do nothing (from zero to zero) */
  }

  else for (k = 0; k < n_b; k++) { // copy 0 or 1 values 
    if ( ((bvec[k]->bs[i] >> j) & 1LL) ) bipartition_set_lowlevel   (bvec[k], i2, j2); // will check if n_ones change or not 
    else                                 bipartition_unset_lowlevel (bvec[k], i2, j2);
  }
}

void
bipartition_resize_vector (bipartition *bvec, int n_b)
{
  int k, i = bvec[0]->n->ints - 1;
  for (k = 0; k < n_b; k++) { bvec[k]->bs[i] &= bvec[0]->n->mask; bipartition_count_n_ones (bvec[k]); }
}

