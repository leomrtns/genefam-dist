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
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULLAR PURPOSE.  See the GNU General Public License for more 
 * details (file "COPYING" or http://www.gnu.org/copyleft/gpl.html).
 */

/* leo@uvigo: In most PRNGs which need a struct (vector and counter) the initial state of the vector is assumed 
 * "in equilibrium" (in a good initial state). The idea from sprng parallel implementation is to replace the "poor"
 * vector initialization (implemented in GSL, for instance) by an explicit sampling from another PRNG, since 
 * distinct vectors (of size 2^14+1 in case of GFSR4 e.g.) will lead to independent streams. Besides the weak point
 * of shift register generators is the possible correlation of vector elements given the single seed. 
 *
 * My idea was to increase the randomness by combining several PRNGs, in the initialization and the higher level
 * functions (in random_number.c). The algorithms are approximately in order of complexity, so that functions at the
 * bottom should not call functions at the top (at the risk of overflowing the stack ;). */

#include "random_number_gen.h"
#include "random_number_aux.h"

void rng_set_marsaglia_constants (uint32_t *m, uint32_t s1, uint32_t s2);

/*
inline uint64_t rng_get_taus (rng_taus_struct *r);
inline uint64_t rng_get_xorshift (rng_xorshift_struct *r);
inline uint64_t rng_get_mt19937 (rng_mt19937_struct *r);
inline uint32_t rng_get_mt19937ar (rng_mt19937ar_struct *r);
inline uint32_t rng_get_gfsr4 (rng_gfsr4_struct *r);
*/
inline uint64_t
rng_get_taus (rng_taus_struct *r)
{ /* r->x[30] (6 vectors of size 4 or 5) */
  int i;
  uint64_t combined = 0ULL, 
           *A  = r->x, /* different (non-overlapping) elements of same vector */
           *q  = r->x + (2 * r->n), /* for this generator only r->n doesn't mean current iteration   */
           *rr = r->x + (3 * r->n), /* it means the number of elements (4 or 5, depending on stream) */
           *C  = r->x + (4 * r->n), 
           *s  = r->x + (5 * r->n);
  /* obs: "r" from the article is called "rr" here since I already called "r" the RNG... */

  for (i=0; i < r->n; i++) { /* order in TAUSWORTHE define: (A,q,r,C,s) where C=2^64 - 2^{64-k} and r=k-s */
    A[i] = ((((A[i] & C[i]) << s[i]) & 0xffffffffffffffffULL) ^ 
                   ((((A[i] << q[i]) & 0xffffffffffffffffULL)^A[i]) >> rr[i]));
    combined ^= A[i];
  }
  return combined;
}

void
rng_set_taus (rng_taus_struct *r, uint64_t seed, int stream)
{ /* r->x[30] (6 vectors of size 4 or 5) */
  int i;
  uint64_t *A, *k, *q; 
 
  rng_set_stream_taus (r, stream);
 
  A = r->x; /* different (non-overlapping) elements of same vector */
  k = r->x +      r->n;
  q = r->x + (2 * r->n);

  if (!seed) seed = (uint64_t) marsaglia_constants[stream % 80U]; /* any number should work */ 

  seed = rng_randomize_array_64bits (r->x, r->n, seed, true); /* initialize vector */
  rng_randomize_array_64bits (r->x, r->n, seed, false);

  for (i=0; i < r->n; i++) /* Initial values should be larger or equal to 2^(64-k) */
    if (A[i] < ( 1ULL << (64 - k[i]))) A[i] += (1ULL << (64 - k[i]));

  for (i=0; i < r->n; i++) /* may not be necessary for our stream parameters */
    A[i] = ((((A[i] <<  q[i]) ^ A[i]) >> k[i]) ^A[i]);

  /* warm up */
  for (i=0; i < 10; i++) rng_get_taus (r);
}

void
rng_set_stream_taus (rng_taus_struct *r, int stream_number)
{
  int i;
  uint64_t *k, *q, *rr, *C, *s;
  
  stream_number %= 150; /* we have 150 distinct streams, 44 with five components and 106 with four */
  if (stream_number < 44) r->n = 5;
  else r->n = 4;

  k  = r->x +       r->n;/* different (non-overlapping) elements of same vector */
  q  = r->x + (2 * r->n);
  rr = r->x + (3 * r->n);
  C  = r->x + (4 * r->n);
  s  = r->x + (5 * r->n);

  /* number of elements of each stream: 20, 24, 4, 2, 92, 8 */
  if (stream_number < 20) { /* table 7 */
    for (i=0; i < r->n; i++) { q[i] = qTable76[0][i]; k[i] = kTable76[0][i]; s[i] = sTable76[stream_number][i]; }
  }
  else if (stream_number < 44) { /* table 6 */
    for (i=0; i < r->n; i++) { q[i] = qTable76[1][i]; k[i] = kTable76[1][i]; s[i] = sTable76[stream_number][i]; }
  }
  else if (stream_number < 48) { /* table 5 (1 to 4) */
    int sn = stream_number - 44;
    for (i=0; i < r->n; i++) { q[i] = qTable543[0][i]; k[i] = kTable543[0][i]; s[i] = sTable543[sn][i]; }
  }
  else if (stream_number < 50) { /* table 5 (5 to 6) */
    int sn = stream_number - 44;
    for (i=0; i < r->n; i++) { q[i] = qTable543[1][i]; k[i] = kTable543[1][i]; s[i] = sTable543[sn][i]; }
  }
  else if (stream_number < 142) { /* table 4 */
    int sn = stream_number - 44;
    for (i=0; i < r->n; i++) { q[i] = qTable543[2][i]; k[i] = kTable543[2][i]; s[i] = sTable543[sn][i]; }
  }
  else { /* (stream_number < 150) table 3 */
    int sn = stream_number - 44;
    for (i=0; i < r->n; i++) { q[i] = qTable543[3][i]; k[i] = kTable543[3][i]; s[i] = sTable543[sn][i]; }
  }

  for (i=0; i < r->n; i++) { /* auxiliary variables, common to all streams ("rr" is original "r" in article) */
    rr[i] = k[i] - s[i];
    C[i] = Cmask[k[i] - 36]; /* we start at k=36 and k < 64 always */
  }
}

inline uint64_t
rng_get_xorshift (rng_xorshift_struct *r)
{ /* x[65]  but last one (r->x[64]) is not part of the stream (is aux variable) */
  int i;
  uint64_t t, v;

  i = r->n = (r->n + 1) & 63;
  t = r->x[i];
  v = r->x[(i+11) & 63];     /* Index is (i-53) mod 64 */
  t ^= t << 33;  t ^= t >> 26; /* (I + L^a)(I + R^b) */
  v ^= v << 27;  v ^= v >> 29; /* (I + L^c)(I + R^d) */
  r->x[i] = (v ^= t) & 0xffffffffffffffffULL; /* Update circular array */
  /* 0x61c8864680b583ebULL = odd approximation to 2**64*(3-sqrt(5))/2. */
  r->x[64] += 0x61c8864680b583ebULL;/* Update Weyl generator */

  return (v + (r->x[64] ^ (r->x[64] >> 27))) & 0xffffffffffffffffULL;
}

void
rng_set_xorshift (rng_xorshift_struct *r, uint64_t seed)
{ /* x[65]  but last one (r->x[64]) is not part of the stream (is aux variable) */
  int i, j;
  uint64_t t, s1 = seed;

  if (!seed) seed = 1064612766ULL; 

  /* Avoid correlations for close seeds; Recurrence has period 2**64-1 */
  for (i=0; i < 64; i++) rng_get_brent_64bits (&seed);

  for (r->x[64] = seed, i = 0; i < 64; i++) { /* Initialise circular array */
    rng_get_brent_64bits (&seed);
    r->x[i] = (seed + (r->x[64] += 0x61c8864680b583ebULL)) & 0xffffffffffffffffULL;
  }

  rng_randomize_array_64bits (r->x, 64, s1, false); /* increase randomness by concatenation */

  for (i = 63, j = 0; j < 256; j++) { /* Discard first 256 results */
    i = (i+1) & 63;
    t = r->x[i];
    seed = r->x[(i+11) & 63];               /* Index is (i-53) mod 64 */
    t ^= t << 33; t ^= t >> 26;             /* (I + L^a)(I + R^b) */
    seed ^= seed << 27; seed ^= seed >> 29; /* (I + L^c)(I + R^d) */
    r->x[i] = (seed ^ t) & 0xffffffffffffffffULL; /* Update circular array */
  }
  r->n = i;
  return;
}

inline uint64_t
rng_get_mt19937 (rng_mt19937_struct *r)
{
  static const uint64_t mag01[2]={ 0ULL, 0xB5026F5AA96619E9ULL}; /* this is magic vector, don't change */
  uint64_t x;

  if (r->n >= 312) { /* generate all 312 words at once */
    int i;
    for (i = 0; i < 156; i++) {
      x = (r->x[i] & 0xFFFFFFFF80000000ULL)| (r->x[i+1] & 0x7FFFFFFFULL);
      r->x[i] = r->x[i+156] ^ (x >> 1) ^ mag01[(int)(x & 1ULL)];
    }
    for (; i < 311; i++) {
      x = (r->x[i] & 0xFFFFFFFF80000000ULL) | (r->x[i+1] & 0x7FFFFFFFULL);
      r->x[i] = r->x[i-156] ^ (x >> 1) ^ mag01[(int)(x & 1ULL)];
    }
    x = (r->x[311] & 0xFFFFFFFF80000000ULL) | (r->x[0] & 0x7FFFFFFFULL);
    r->x[311] = r->x[155] ^ (x >> 1) ^ mag01[(int)(x & 1ULL)];
    r->n = 0; 
  }
  x = r->x[r->n++];

  x ^= (x >> 29) & 0x5555555555555555ULL;
  x ^= (x << 17) & 0x71D67FFFEDA60000ULL;
  x ^= (x << 37) & 0xFFF7EEE000000000ULL;
  x ^= (x >> 43);

  return x;
}

void
rng_set_mt19937 (rng_mt19937_struct *r, uint64_t seed)
{
  int i;
  if (!seed) seed = 548919650218ULL;

  r->n = 313;
  r->x[0] = seed;
  for (i = 1; i < 312; i++) r->x[i] = (6364136223846793005ULL * (r->x[i-1] ^ (r->x[i-1] >> 62)) + i);

  rng_randomize_array_64bits (r->x, 312, seed, false); /* increase randomness by concatenation */
  rng_twist_array_64bits (r->x, 312, seed);
}

inline uint32_t
rng_get_mt19937ar (rng_mt19937ar_struct *r)
{
  static const uint32_t mag01[2]={ 0ULL, 0x9908b0dfUL}; /* this is magic vector, don't change */
  uint32_t y;
  if (r->n >= 624) { /* generate N words at one time */
    int i;

    for (i = 0; i < 227; i++) {
      y = (r->x[i] & 0x80000000UL) | (r->x[i+1] & 0x7fffffffUL);
      r->x[i] = r->x[i+397] ^ (y >> 1) ^ mag01[y & 1UL];
    }
    for (; i < 623; i++) {
      y = (r->x[i] & 0x80000000UL) | (r->x[i+1] & 0x7fffffffUL);
      r->x[i] = r->x[i-277] ^ (y >> 1) ^ mag01[y & 1UL];
    }
    y = (r->x[623] & 0x80000000UL) | (r->x[0] & 0x7fffffffUL);
    r->x[623] = r->x[396] ^ (y >> 1) ^ mag01[y & 1UL];

    r->n = 0;
  }

  y = r->x[r->n++];

  /* Tempering */
  y ^= (y >> 11);
  y ^= (y << 7) & 0x9d2c5680UL;
  y ^= (y << 15) & 0xefc60000UL;
  y ^= (y >> 18);

  return y;
}

void
rng_set_mt19937ar (rng_mt19937ar_struct *r, uint32_t seed)
{
  int i;

  if (!seed) seed = 19650218UL;

  r->n = 625;
  r->x[0]= seed;
  for (i = 1; i < 624; i++) r->x[i] = (1812433253UL * (r->x[i-1] ^ (r->x[i-1] >> 30)) + i);
  
  seed = rng_randomize_array_32bits (r->x, 624, seed, false); /* increase randomness by concatenation */
  rng_twist_array_32bits (r->x, 624, seed);
}

inline uint32_t
rng_get_gfsr4 (rng_gfsr4_struct *r)
{ /* r->x[16384] */
  r->n = ((r->n)+1) & 16383;
  return r->x[r->n] =
  ( r->x[(r->n + 15913) & 16383] ^ r->x[(r->n + 14798) & 16383] ^ 
    r->x[(r->n + 9396)  & 16383] ^ r->x[(r->n + 6695)  & 16383]) & 0xffffffffUL;
}

void
rng_set_gfsr4 (rng_gfsr4_struct *r, uint32_t seed)
{ /* r->x[16384] */
  int i;
  rng_tt800_struct tt;

  if (!seed) seed = 627749721UL;

  rng_set_tt800 (&tt, seed); /* initialize tt800 generator */

  for (i = 0; i < 16384; i++) r->x[i] = (rng_get_tt800 (&tt) ^ rng_get_cong_many (&seed)) & 0xffffffffUL; 
  seed = rng_randomize_array_32bits (r->x, 16384, seed, false); /* increase randomness by concatenation */
  for (i = 0; i < 16384; i++) r->x[i] ^= rng_get_tt800 (&tt);
  r->n = 0;
}

/* * * below, generators whose initial states depend only on simple functions (single-variable PRNGs etc.) */ 

uint32_t
rng_get_diaconis (rng_diaconis_struct *r)
{ /*r->x[128] */
  /* two implicid bits version. Period lenght: (2**R + 1) * 2**31. */
  uint32_t b0, sr, ss, br, bs;
  
  /* fib(n) = fib(n-R) * fib(n-S); with all fib() odd. */
  r->n--;
  br = r->x[(r->n + 127) & 127];
  bs = r->x[(r->n + 30)  & 127];
  sr = br & 1; br ^= sr;
  ss = bs & 1; bs ^= ss;
  b0 = 4 * br * bs;
  if (sr) br *= 3;
  if (ss) bs *= 3;
  b0 += br + bs + sr + ss;
  r->x[r->n & 127] = b0;
  return (b0 + (b0 >> 16)) & 0xffffffffUL; /* low bit improvement */
}

void
rng_set_diaconis (rng_diaconis_struct *r, uint32_t seed)
{ /*r->x[128] */
  int i;
  uint32_t m[4];
  r->n = 0;

  if (!seed) seed = 1372460312UL;

  rng_set_marsaglia (m, seed);
  rng_get_shr (&seed);

  for (i = 0; i < 127; i++) r->x[i] = rng_get_brent (&seed) + rng_get_marsaglia (m);
  rng_randomize_array_32bits (r->x, 127, seed, false); /* increase randomness by concatenation */
}

/* Makoto Matsumoto & Y. Kurita, Twisted GFSR Generators II, ACM Trans. Model. Comput. Simul., 4 (1994) 254-266 */
uint32_t
rng_get_tt800 (rng_tt800_struct *r)
{ /* r->x[25] */
  static const uint32_t mag01[2]={ 0x0, 0x8ebfd028}; /* this is magic vector, don't change */
  uint32_t y;

  if (r->n >= 25) { /* generate N words at one time */
    int i;
    for (i = 0; i < 18; i++) r->x[i] = r->x[i+7]  ^ (r->x[i] >> 1) ^ mag01[r->x[i] % 2];
    for (; i < 25; i++)      r->x[i] = r->x[i-18] ^ (r->x[i] >> 1) ^ mag01[r->x[i] % 2];
    r->n = 0;
  }
  y = r->x[r->n];
  y ^= (y << 7)  & 0x2b5b2500UL; /* s and b, magic vectors */
  y ^= (y << 15) & 0xdb8b0000UL; /* t and c, magic vectors */
  r->n++;
  return (y ^ (y >> 16)) & 0xffffffffUL;
}

void
rng_set_tt800 (rng_tt800_struct *r, uint32_t seed)
{ /* r->x[25] */
  int i;
  uint32_t s1, s2;

  if (!seed) seed = 1529210297UL;
  s1 = seed;
  s2 = rng_get_cong_many (&s1);

  r->n = 26;
  for (i = 0; i < 25; i++) r->x[i] = (rng_get_shr (&s2) ^ rng_get_brent (&seed)) + rng_get_cong (&s1);
  seed = rng_randomize_array_32bits (r->x, 25, seed, false); /* increase randomness by concatenation */
  rng_randomize_array_32bits (r->x, 25, seed, false);
}

uint32_t
rng_get_lfib4 (rng_lfib4_struct *r)
{ /* r->x[256] */  /* set = vector filling */
  r->n = (r->n + 1) & 255;
  return r->x[r->n] = r->x[r->n & 255] + r->x[(r->n + 58) & 255] + r->x[(r->n + 119) & 255] + r->x[(r->n + 178) & 255];
}

void
rng_set_lfib4 (rng_swb_struct *r, uint32_t seed)
{
  int i;
  uint32_t s1, s2, m[4];

  if (!seed) seed = 2642725982UL;
  s2 = seed;

  r->n = 0;
  s1 = biomcmc_hashint_7 (seed); /* arbitrary shuffling */
  for (i = 0; i < 32; i++) { rng_get_cong (&s2); rng_get_cong (&s1); rng_get_brent (&s1); }

  rng_set_marsaglia (m, s1);

  for (i = 0; i < 192; i++)  r->x[i]  = rng_get_marsaglia (m) ^ rng_get_shr (&s1) ^ rng_get_std31 (&seed);
  for (i = 64; i < 256; i++) r->x[i] ^= rng_get_marsaglia (m) ^ rng_get_shr (&s2) ^ rng_get_std31 (&seed);
  rng_randomize_array_32bits (r->x, 256, seed, false); /* increase randomness by concatenation */
}

uint32_t
rng_get_swb (rng_swb_struct *r)
{ /* r->x[258] (x[256] + two aux variables) */
  r->x[256] = r->x[(r->n + 15) & 255]; /* x in original algorithm */
  r->x[(r->n + 237) & 255] = r->x[256] - (r->x[257] = r->x[(r->n+1) & 255] + (r->x[256] < r->x[257]));
  r->n = (r->n + 1) & 255;
  return r->x[r->n];
}

void
rng_set_swb (rng_swb_struct *r, uint32_t seed)
{
  int i;
  uint32_t s1;
  r->n = 0;
  r->x[256] = r->x[257] = 0UL;

  if (!seed) seed = 904977562UL;

  s1 = seed;
  for (i=0; i < 32; i++) rng_get_brent (&s1); 
  
  for (i = 0; i < 256; i++) r->x[i] = rng_get_shr (&s1) ^ rng_get_std31 (&seed);
  rng_randomize_array_32bits (r->x, 256, seed, false); /* increase randomness by concatenation */
}

/* Panneton, L'Ecuyer, and Matsumoto WELLRNG1024a */
uint32_t
rng_get_well1024 (rng_well1024_struct *r)
{ /* r->x[32] */ /* set = vector filling */
  uint32_t z0, z1, z2;
  z0 = r->x[(r->n + 31) & 0x0000001fU];
  z1 = r->x[r->n] ^ (r->x[(r->n + 3) & 0x0000001fU] ^ (r->x[(r->n + 3) & 0x0000001fU] >> 8));
  z2 = ((r->x[(r->n + 24) & 0x0000001fU] ^ (r->x[(r->n + 24) & 0x0000001fU] << 19)) ^
        (r->x[(r->n + 10) & 0x0000001fU] ^ (r->x[(r->n + 10) & 0x0000001fU] << 14)));
  r->x[r->n] = z1 ^ z2;
  r->x[(r->n + 31) & 0x0000001fU] = (z0 ^ (z0 << 11)) ^ (z1 ^ (z1 << 7)) ^ (z2 ^ (z2 << 13));
  r->n = (r->n + 31) & 0x0000001fU;
  return r->x[r->n];
}

void
rng_set_well1024 (rng_well1024_struct *r, uint32_t seed)
{
  int i;
  uint32_t s1, s2;
  r->n = 0;

  if (!seed) seed = 3519793928UL;

  s1 = biomcmc_hashint_8 (seed);
  s2 = rng_get_shr (&s1); 
  for (i = 0; i < 32; i++) r->x[i] = rng_get_shr (&s1) ^ rng_get_std31 (&s2);
  seed = rng_randomize_array_32bits (r->x, 32, seed, false); /* increase randomness by concatenation */
  rng_twist_array_32bits (r->x, 32, seed);
}

/* * * Simple generators (single value or simple vectors of size 2 or 4) */

uint32_t
rng_get_gamerand (uint32_t *game)
{ /* game[2] */
  game[0] = (game[0] << 16) + (game[0] >> 16); game[0] += game[1]; game[1] += game[0];
  return game[0];
}

void
rng_set_gamerand (uint32_t *game, uint32_t seed)
{ /* game[2] */
  int i;

  if (!seed) seed = 7584631UL; 

  game[0] = seed; 
  game[1] = biomcmc_hashint_9 (seed);
  for (i = 0; i < 32; i++) rng_get_brent (game + 1);
}

/* Marsaglia's Super-Duper (two-components multiply-with-carry) algortihm */
uint32_t
rng_get_marsaglia (uint32_t *m)
{ /* m[4] */
  return (((m[0] = m[2] * (m[0] & 0xffffUL) + (m[0] >> 16)) << 16) + 
          ((m[1] = m[3] * (m[1] & 0xffffUL) + (m[1] >> 16 )) & 0xffffUL)) & 0xffffffffUL;
}

void
rng_set_marsaglia (uint32_t *m, uint32_t seed)
{ /* m[4] */
  m[0] = seed; 
  m[1] = 1UL + biomcmc_hashint_9 (seed);
  if (!m[0]) m[0] = 362436069UL;
  if (!m[1]) m[1] = 521288629UL; 
  
  if (m[0] == m[1]) m[1] *= 69069UL;

  rng_set_marsaglia_constants (m, m[0], m[1]);
}

void
rng_set_marsaglia_constants (uint32_t *m, uint32_t s1, uint32_t s2)
{
  /* we must choose two distinct marsaglia_constants[], from a pool of 81. We have 6400 possible streams, since one
   * constant we reserve for the diagonals. */
  s1 = biomcmc_hashint_9 (s1); s1 %= 80;
  s2 = biomcmc_hashint_9 (s2); s2 %= 80;
  if (s1 == s2) s2 = 80;
  m[2] = marsaglia_constants[s1];
  m[3] = marsaglia_constants[s2];
}

uint64_t
rng_get_std61 (uint64_t *x)
{
  (*x) = ((*x) >> 31) + (((*x) << 30) & 0x1fffffffffffffffULL) - ((*x) >> 42) - (((*x) << 19) & 0x1fffffffffffffffULL);
  if ((int64_t) (*x) < 0) (*x) += 0x1fffffffffffffffULL;
  return (*x);
}

/* The Minimal Portable Random Number Generator (32 bits)
   a = 7^5 = 16807   m = 2^31 - 1 = 2147483647 = 0x7fffffff ;   x[n+1] = a * x[n] (mod m) */
uint32_t
rng_get_std31 (uint32_t *x)
{
  uint32_t zh, zl, z;
  z = (*x) << 1; 
  zl = z & 0xffffUL; 
  zh = z >> 16; 
  zl *= 48271; zh *= 48271; 
  zh += zl >> 16; 
  zl = (zl & 0xffffUL) + (zh << 16);
  zh = (zh >> 16) << 1UL; 
  zl += zh;
  if (zh > zl) zl += 2UL;
  (*x) = zl >> 1;
  return (*x);
}

uint32_t
rng_get_shr (uint32_t *x)
{
  (*x) ^= ((*x) << 17); (*x) ^= ((*x) >> 13), (*x) ^= ((*x) << 5);
  return (*x);
}

uint32_t
rng_get_brent (uint32_t *x)
{
  (*x) ^= (*x) << 10; (*x) ^= (*x) >> 15; (*x) ^= (*x) << 4;  (*x) ^= (*x) >> 13;
  return *x;
}

uint64_t
rng_get_brent_64bits (uint64_t *x)
{
  (*x) ^= (*x) << 10; (*x) ^= (*x) >> 15; (*x) ^= (*x) << 4;  (*x) ^= (*x) >> 13;
  return *x;
}

uint32_t
rng_get_cong (uint32_t *x)
{
  uint32_t b1, b2; /* we use only the leading half bits (the lsat ones are too regular */
  b1 = ((*x) = (69069UL * (*x)) + 1234567) >> 16;
  b2 = ((*x) = (69069UL * (*x)) + 1234567) & 0xffff0000UL;
  return (b1 | b2); /* notice that here we don't update x with (b1|b2) */
}

uint32_t
rng_get_cong_many (uint32_t *x)
{ /* quick-and-dirty (well, not so quick)  from GSL */
  uint32_t j, t = 0UL, bit =  0x80000000UL;
  for (j = 0; j < 32; j++) {
    rng_get_cong (x);
    if ((*x) & 0x80000000UL) t |= bit;
    bit >>= 1;
  }
  return t; /* we don't return x, despite it's being being updated */
}

uint32_t
rng_twist_array_32bits (uint32_t *a, uint32_t n_a, uint32_t seed)
{
  uint32_t i, im, mars[4];

  /* initialize Marsaglia's Super-Duper generator */
  mars[0] = biomcmc_hashint_8 (seed) + 1UL; 
  mars[1] = seed;
  mars[1] = rng_get_cong_many (mars + 1);
  rng_set_marsaglia_constants (mars, mars[0], mars[1]);

  if (!a[0]) a[0] = (1ULL << 30) - n_a;
  a[0] += rng_get_std31 (&seed);

  for (i = 0; i < n_a; i++) {
    im = (i + n_a - 1) % n_a; /* im = i - 1 but warping around (so that zero minus one is last element) */ 
    a[i] = (a[i] ^ ((a[im] ^ (a[im] >> 30)) * 1664525UL)) + rng_get_std31 (&seed) + rng_get_marsaglia (mars);
  }

  seed = 69069UL * ((seed >> 6) + ((uint32_t) n_a)); /* change the seed */

  for (i = 0; i < n_a; i++) {
    im = (i + n_a - 2) % n_a; /* im = i - 2 but with warp-around */ 
    a[i] = (a[i] ^ ((a[im] ^ (a[im] >> 30)) * 1566083941UL)) ^ rng_get_std31 (&seed) ^ rng_get_marsaglia (mars);
  }
  if (!a[0]) a[0] = (1ULL << 30) + ((uint32_t) n_a);

  return seed;
}

uint64_t
rng_twist_array_64bits (uint64_t *a, uint32_t n_a, uint64_t seed)
{ /* modified from MT19937 (http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html) */
  uint32_t i, im;
  uint64_t s1 = biomcmc_hashint_64bits (seed) + ((uint64_t) n_a);

  if (!a[0]) a[0] = (1ULL << 63);
  a[0] += rng_get_std61 (&seed);

  for (i = 0; i < n_a; i++) {
    im = (i + n_a - 1) % n_a; /* im = i - 1 but warping around (so that zero minus one is last element) */ 
    a[i] = (a[i] ^ ((a[im] ^ (a[im] >> 62)) * 3935559000370003845ULL)) + rng_get_std61 (&seed) + i;

    for (im = 0; im < 4; im++) rng_get_brent_64bits (&s1);
    a[i] ^= rng_get_std61 (&seed) ^ rng_get_brent_64bits (&s1);
  }

  seed = 69069ULL * ((seed >> 6) + 1ULL); /* change the seed */

  for (i = 0; i < n_a; i++) {
    im = (i + n_a - 2) % n_a; /* im = i - 2 but with warp-around */ 
    a[i] = (a[i] ^ ((a[im] ^ (a[im] >> 62)) * 2862933555777941757ULL)) - i;
    a[i] ^= rng_get_std61 (&seed);
  }
  
  if (!a[0]) a[0] = (1ULL << 63);
  
  return seed; /* so that this function can be called several times */
}

uint32_t
rng_randomize_array_32bits (uint32_t *a, uint32_t n_a, uint32_t seed, bool first_time)
{
  uint32_t i, t, m[4];

  if (!seed) seed = (69069 * n_a) + 69069;

  rng_set_marsaglia (m, seed);
  t = rng_get_marsaglia (m); 
  seed = rng_get_marsaglia (m);
  for (i = 0; i < 16; i++) rng_get_brent (&seed);

  if (first_time) for (i = 0; i < n_a; i++) 
    a[i] = (rng_get_shr (&t) ^ rng_get_marsaglia (m) ^ rng_get_std31 (&t)) ^ rng_get_cong (&seed);

  rng_get_brent (&seed);

  for (i = 0; i < n_a; i++) 
    a[i] ^= (rng_get_shr (&t) ^ rng_get_marsaglia (m) ^ rng_get_std31 (&t)) + rng_get_cong (&seed);

  return seed;
}

uint64_t
rng_randomize_array_64bits (uint64_t *a, uint32_t n_a, uint64_t seed, bool first_time)
{
  uint64_t u1, u2, u3;
  uint32_t i, t, s, m[4];

  if (!seed) seed = (uint32_t) (69069 * n_a);

  m[0] = (uint32_t) (seed >> 32); /* 32 bits component: Marsaglia's multiply-with-carry */
  m[1] = (uint32_t) (seed & 0xffffffffUL);
  rng_set_marsaglia (m, m[0]);
  t = rng_get_marsaglia (m); /* 32 bits component: Marsaglia's 3-shift register */
  s = rng_get_marsaglia (m); /* 32 bits component: Fishman and Moore's std31 */

  if (first_time) for (i = 0; i < n_a; i++) {
    u1 = (((uint64_t) (rng_get_shr (&t))) << 32)      | ((uint64_t) (rng_get_shr (&t))); /* 64 bits filled */
    u2 = (((uint64_t) (rng_get_marsaglia (m))) << 32) | ((uint64_t) (rng_get_marsaglia (m)));
    u3 = (((uint64_t) (rng_get_std31 (&s))) << 32)    | ((uint64_t) (rng_get_std31 (&s)));
    a[i] = u1 ^ u2 ^ u3 ^ rng_get_std61 (&seed);
  }

  m[0] += n_a;

  for (i = 0; i < n_a; i++) {
    u1 = (((uint64_t) (rng_get_shr (&t))) << 32)      | ((uint64_t) (rng_get_shr (&t))); /* 64 bits filled */
    u2 = (((uint64_t) (rng_get_marsaglia (m))) << 32) | ((uint64_t) (rng_get_marsaglia (m)));
    u3 = (((uint64_t) (rng_get_std31 (&s))) << 32)    | ((uint64_t) (rng_get_std31 (&s)));
    a[i] ^= (u1 + u2) ^ u3 ^ rng_get_std61 (&seed);
  }
  return seed;
}


