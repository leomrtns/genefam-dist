/* 
 * This file is part of guenomuULL, a hierarchical Bayesian procedure to estimate the distribution of species trees based
 * on multi-gene families data.
 * Copyright (C) 2009  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 *
 * Guenomu is free software; you can redistribute it and/or modify it under the terms of the GNU General Public 
 * License as published by the Free Software Foundation; either version 2 of the LicenseULL, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be usefulULL, but WITHOUT ANY WARRANTY; without even the implied 
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULLAR PURPOSE.  See the GNU General Public License for more 
 * details (file "COPYING" or http://www.gnu.org/copyleft/gpl.html).
 */

/*! \file random_number_aux.h
 *  \brief Variables and structures local to random_number.c (should be opaque to user) */

#ifndef _biomcmc_random_number_gen_h_
#define _biomcmc_random_number_gen_h_

#include "hashtable.h"

/* Tausworthe linear feedback shift register PRNG */
/* 6 vectors (A,k,q,r,C,s) of size at most 5 (4 or 5, depending on stream) */
typedef struct { uint64_t x[30];    int n; } rng_taus_struct; 

/* \brief XOR-shift (R. P. Brent's xorgens) */
/* the vector itself has only 64 elements but the 65th one is the Weyl generator; */
typedef struct { uint64_t x[65];    int n; } rng_xorshift_struct; 
 
/*! \brief MT19937-64, the Mersenne Twister for 64 bits */
typedef struct { uint64_t x[312];   int n; } rng_mt19937_struct; 

/*! \brief MT19937, the original Mersenne Twister (for 32 bits); the name "ar" comes from "array" */
typedef struct { uint32_t x[624];   int n; } rng_mt19937ar_struct;

/*! \brief GFSR4 implementation from GSL */
typedef struct { uint32_t x[16384]; int n; } rng_gfsr4_struct; /*16384 = 2^14 */

/*! \brief Persi Diaconis' lagged Fibonacci */ 
typedef struct { uint32_t x[128];   int n; } rng_diaconis_struct;

/*! \brief tt800 (small cousin of MT19937, the Mersenne Twister) */
typedef struct { uint32_t x[25];    int n; } rng_tt800_struct;

/*! \brief Marsaglia's LFIB4 lagged Fibonacci using addition */ 
typedef struct { uint32_t x[256];   int n; } rng_lfib4_struct;

/*! \brief Marsaglia's Subtract-with-borrow generator */ 
typedef struct { uint32_t x[258];   int n; } rng_swb_struct;

/* Panneton, L'Ecuyer, and Matsumoto WELLRNG1024a */
typedef struct { uint32_t x[32];    int n; } rng_well1024_struct;

/*! \brief Maximally equidistributed combined Tausworthe (Linear Feedback Shift Register) random number generator.
 * (for more info see file docs/random_number_generation.txt or http://www.gnu.org/software/gsl).  */
/*  https://gcc.gnu.org/onlinedocs/gcc/Inline.html --> inline behaviour changed in GNU11 which is the dafault for GCC5  */
extern uint64_t rng_get_taus (rng_taus_struct *r);
void            rng_set_taus (rng_taus_struct *r, uint64_t seed, int stream);
void            rng_set_stream_taus (rng_taus_struct *r, int stream_number);

extern uint64_t rng_get_xorshift (rng_xorshift_struct *r);
void            rng_set_xorshift (rng_xorshift_struct *r, uint64_t seed);

/*! \brief MT19937-64, the Mersenne Twister for 64 bits */
extern uint64_t rng_get_mt19937 (rng_mt19937_struct *r);
void            rng_set_mt19937 (rng_mt19937_struct *r, uint64_t seed);

/*! \brief MT19937, the Mersenne Twister for 32 bits */
extern uint32_t rng_get_mt19937ar (rng_mt19937ar_struct *r);
void            rng_set_mt19937ar (rng_mt19937ar_struct *r, uint32_t seed);

/*! \brief gfsr lagged-Fibonacci generator (two-tap Generalised Feedback Shift Register). */
extern uint32_t rng_get_gfsr4 (rng_gfsr4_struct *r);
void            rng_set_gfsr4 (rng_gfsr4_struct *r, uint32_t seed);

uint32_t rng_get_diaconis (rng_diaconis_struct *r);
void     rng_set_diaconis (rng_diaconis_struct *r, uint32_t seed);

/* TT800 (period = 2^800)  
 * Makoto Matsumoto & Y. Kurita, Twisted GFSR Generators II, ACM Trans. Model. Comput. Simul., 4 (1994) 254-266 */
uint32_t rng_get_tt800 (rng_tt800_struct *r);
void     rng_set_tt800 (rng_tt800_struct *r, uint32_t seed);

/* Marsaglia's four-lag generator using addition */
uint32_t rng_get_lfib4 (rng_lfib4_struct *r);
void     rng_set_lfib4 (rng_swb_struct *r, uint32_t seed);

/* Marsaglia's subtract-with-borrow generator with long period (2^7578) but fails at birthday tests */
uint32_t rng_get_swb (rng_swb_struct *r);
void     rng_set_swb (rng_swb_struct *r, uint32_t seed);

/* Panneton, L'Ecuyer, and Matsumoto WELLRNG1024a */
uint32_t rng_get_well1024 (rng_well1024_struct *r);
void     rng_set_well1024 (rng_well1024_struct *r, uint32_t seed);

/* simple algorithms (up to four variables, no vector initialization) */

uint32_t rng_get_gamerand (uint32_t *game);
void     rng_set_gamerand (uint32_t *game, uint32_t seed);

/* George Marsaglia's Two-multiply with carry generator (period ~  2^60) */
uint32_t rng_get_marsaglia (uint32_t *m);
void     rng_set_marsaglia (uint32_t *m, uint32_t seed);

/* Multiplicative, Congruential Random-Number Generators with Multipliers +-2^k1 +- 2^k2 and Modulus 2^p - 1,
   ACM Trans. Math. Software 23 (1997) 255-265. Pei-Chi Wu. */ 
uint64_t rng_get_std61 (uint64_t *x);

/* The Minimal Portable Random Number Generator (32 bits)
   G. S. Fishman, L. R. Moore; An statistical exhaustive analysis of multiplicative congruential random number 
   generators with modulus 2^31-1, SIAM J. Sci. Statist. Comput., 7 (1986) 24-45. Erratum, ibid, 7 (1986) 1058 */
uint32_t rng_get_std31 (uint32_t *x);

/* Marsaglia's 3-shift-register generator (period 2^32-1) */
uint32_t rng_get_shr (uint32_t *x);

/* Brent auxiliary function to decrease correlation between seeds */
uint32_t rng_get_brent (uint32_t *x);
uint64_t rng_get_brent_64bits (uint64_t *x);

/* congruential quick-and-dirty generator. */
uint32_t rng_get_cong      (uint32_t *x);
uint32_t rng_get_cong_many (uint32_t *x);

/* randomize an initialized vector following Nishimura and Matsumoto's MT19937 */
uint32_t rng_twist_array_32bits (uint32_t *a, uint32_t n_a, uint32_t seed);
uint64_t rng_twist_array_64bits (uint64_t *a, uint32_t n_a, uint64_t seed);

/* randomize after possibly initializing a vector using XOR concatenation */
uint32_t rng_randomize_array_32bits (uint32_t *a, uint32_t n_a, uint32_t seed, bool first_time);
uint64_t rng_randomize_array_64bits (uint64_t *a, uint32_t n_a, uint64_t seed, bool first_time);

#endif
