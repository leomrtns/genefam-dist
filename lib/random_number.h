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

/*! \file random_number.h
 * \brief Random number generation, with algorithms from the Gnu Scientific Library (GSL) and motivation from the 
 * Scalable Parallel Pseudo Random Number Generators Library (SPRNG). 
 *
 * There are two ways of setting up the pseudo-random number generation: a high-level approach where we simply call one
 * function before start using the generator and another after we finished using it, and a lower-level approach where we
 * must allocate, seed and free the streams by hand. Both approaches are completely incompatible, and there is no check
 * to avoid this mistake so the programmer should use the first, high-level approach unless he understand well the
 * algorithms.
 * The high-level approach is useful when one stream is enough (for instance, a serial program). In this case we use 
 * a (modified) maximally-equidistributed combined Linear Feedback Shift Register (LFSR, or Tausworthe) pseudo-random 
 * number generator (PRNG) whose seed is based on time of day.
 *
 * The low level approach is useful if one needs several uncorrelated streams - for example a parallel program where
 * each node must have its own stream of pseudo-random numbers and a common stream shared by all nodes. For this we
 * implemented a modified <a href="http://www.mcs.anl.gov/~itf/dbpp/text/node119.html">random tree method</a> where
 * a lagged-Fibonacci generator (two-tap Generalised Feedback Shift Register - GFSR4) PRNG is used to generate the 
 * seeds for the (modified) Tausworthe streams. The seed for this GFSR4 generator is given by the time of day and should 
 * be set only once by the program, given its low randomness. The GFSR4 uses a vector of size \f$2^14\f$ which is
 * initialized by an (quick-and-dirty) xorshift randomization of the seed, which makes it sensitive to the choice 
 * of this seed. I included a 
 * <a href="http://ianbullard.squarespace.com/journal/2009/4/28/why-you-should-never-use-rand.html">gamerand 
 * fast randomization</a> over each element of the vector, modified for 64 bits (the algorithm is too simple to require
 * a license and I've seen public domain versions of it; the original disclaimer is GPL-like). A more
 * robust alternative <b>not implemented here</b> would be to use a better PRNG to feed the initial vector of the GFSR4
 * generator. As noticed by <a href="http://www.springerlink.com/content/ek34q1jwla5b6438/">CJK Tan and JAR Blais 
 * (HPCN 2000, LNCS 1823, 127-135)</a> another random tree method could be
 * constructed with parallel GFSR4 streams where the initial states are created by an PRNG - equivalent to our seed
 * generator, but providing not only the seed but all \f$2^{14}\f$ elements to each stream. To ensure independence of
 * streams - since different seeds are just different points of the same stream - we implemented all 150 parameter sets
 * provided by <a href="http://www.iro.umontreal.ca/~lecuyer/myftp/papers/tausme2.ps">
* L'ecuyer (Maths Comput 1999, pp.261)</a>for the Tausworthe generators. So we have at least 150 independent streams, 
* with periods between 
* \f$10^{14}\f$ and \f$10^{35}\f$ - if I understood correctly the interpretation of the number of non-zero solutions of the
* polynomials \f$N_1\f$. 
*
* Our modification to the original Tausworthe (LFSR) algorithms is to combine it with the
* <a href="http://wwwmaths.anu.edu.au/~brent/random.html">generalized Marsaglias's xorshift PRNG called xorgens</a>,
* developed by Richard Brent and available under a GPL license. So the LFSR has one extra component from this
* independent xorshift (combined through a XOR). The seed for the xorshift is the same as for the Tausworthe.
*
* All streams are of 64 bits (dependent on a "long long int" having 64 bits) but 32 bits or even 16 bits are available
* through wrapper functions. If your system does not provide 64 bit integers use the unwrapped versions instead.
* When working with more than one stream at the same time - by the same node, using the example of a 
* parallel program - we must update by hand the variable pointed to by the global variable ::random_number to indicate
* which stream should be used. 
*
* Some original comments from the GSL can be found on "doc/random_number_generation.txt" */

#ifndef _biomcmc_random_number_h_
#define _biomcmc_random_number_h_

#include "random_number_gen.h"

typedef struct biomcmc_rng_struct* biomcmc_rng;

/*! \brief Random number structure (combined Tausworthe algorithm) */
struct biomcmc_rng_struct 
{
  rng_taus_struct taus; /*! \brief Tausworthe linear feedback shift-register from GSL */

  rng_mt19937_struct mt; /*! \brief 64 bits Mersenne Twister from Matsumoto's webpage */

  uint64_t bit32;    /*! \brief temporary values when only 32 bits are necessary */
  bool have_bit32;   /*! \brief when using 32 bits we first check if we have one stored */ 

  double rnorm32, rnorm64;         /*! \brief stored standard normal random values with 32 and 52 bits of precision */
  bool have_rnorm32, have_rnorm64; /*! \brief true if we have normal random variate stored */ 
};

/*! \brief pointer to pseudo-random number generator (should point to real stream, even when there are several) */
extern biomcmc_rng biomcmc_random_number; /* in OSX at least, the linker will bundle all sources into one library,
                                             and thus global variables must be declared "extern" in the headers 
                                             and defined only once in the corresponding .c file (to avoid duplicates).
                                             Functions don't need this because they are implicitly declared "extern" 
                                             Another solution is to declare all global variables as "static" but this 
                                             is dangerous since we may have several static copies. */

/*! \brief High-level setup of a simple random number generator and initialization with a seed (not to be mixed with 
 * other low-level functions that update the seed or allocate memory).
 *
 * The seed may be provided by calling function (mostly for debug) if not zero otherwise it is based on present time 
 * of day, and uses the Tausworthe pseudo-random number generator. 
 * The Tausworthe generator we use has a period of (at least?) \f$10^{35}\f$. This function allocates memory to 
 * global variable ::random_number directly */
void biomcmc_random_number_init (unsigned long long int seed);
/*! \brief High-level finalization (memory release etc.) of the random number environment */ 
void biomcmc_random_number_finalize (void);

/*! \brief Allocate memory for new (Tausworthe + MT19937) generator from a pool of streams */
biomcmc_rng new_biomcmc_rng (unsigned long long int seed, int stream_number);
/*! \brief Release memory occupied by biomcmc_rng:: */
void del_biomcmc_rng (biomcmc_rng r);

/*! \brief Create initial seed based on time, combining time in microseconds and seconds (user-controled seed is not uint64_t) */
unsigned long long int biomcmc_rng_get_initial_seed (void);
/*! \brief Generate a vector of seeds (based on initial one), create and initialize stream with an element of this 
 * vector */
biomcmc_rng new_biomcmc_rng_from_seed (unsigned long long int seed, int stream_number);

/*! \brief Returns a random number from a Standard Normal distribution N(0,1) - prob_distribution.h has general case */
extern double biomcmc_rng_snorm32 (void);
/*! \brief Returns a random number from a Standard Normal distribution with maximum (52 bits) integer precision */
extern double biomcmc_rng_snorm (void);

/*! \brief Returns a random number between 0 and 1 (including 1) with precision \f$\approx 10^{-15}\f$ (52 bits). */
extern double biomcmc_rng_unif (void);
/*! \brief Returns a positive random number between 0 and 1 (excluding 0 and including 1) (52 bits). */
extern double biomcmc_rng_unif_pos (void);
/*! \brief Returns an long integer (64 bits) random number between 0 and n (excluding n), with \f$n < 10^{15}\f$. */
extern uint64_t biomcmc_rng_unif_int64 (uint64_t n);

/*! \brief Returns a random number between 0 and 1 (including 1) with precision \f$\approx 10^{-9}\f$ (32 int bits). */
extern double biomcmc_rng_unif32 (void);
/*! \brief Returns a positive random number between 0 and 1 (including 1) (32 int bits). */
extern double biomcmc_rng_unif_pos32 (void);
/*! \brief Returns an integer (32 bits) random number between 0 and n (excluding n), provided \f$n < 4 10^9\f$ approx. */
extern uint32_t biomcmc_rng_unif_int (uint32_t n);

/* * * low level functions (direct access to generator) */

/*! \brief  new value with 64 random bits */
extern uint64_t biomcmc_rng_get (void);
/*! \brief  new value with 52 random bits as a double precision float */
extern double   biomcmc_rng_get_52 (void);
/*! \brief  new value with 32 random bits */
extern uint32_t biomcmc_rng_get_32 (void);

/*! \brief get current time with maximum precision and soter in vector time[2] */
void biomcmc_get_time (int *time);
/*! \brief returns the floating-point time in seconds elapsed between past[2] and now[2] */
double biomcmc_elapsed_time (int *now, int *past);

#endif
