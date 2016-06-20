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

/*! \file random_number_aux.h */

#ifndef _biomcmc_random_number_aux_h_
#define _biomcmc_random_number_aux_h_

/*! \brief  Five-element streams for L'ecuyer's combined LFSR (Tausworthe) generator */ 
uint64_t sTable76[44][5] = {
  /* table 7 of MathsComput(1999)p261 (from elements 0 to 19 here) - safe for up to \f$10^{35}\f$ samples */
   {9ULL, 34ULL, 5ULL, 26ULL, 18ULL},  {9ULL, 32ULL, 5ULL, 31ULL, 6ULL},   {9ULL, 25ULL, 5ULL, 37ULL, 22ULL},
   {10ULL, 24ULL, 5ULL, 7ULL, 12ULL},  {12ULL, 17ULL, 5ULL, 14ULL, 8ULL},  {12ULL, 40ULL, 5ULL, 16ULL, 22ULL}, 
   {12ULL, 26ULL, 5ULL, 34ULL, 23ULL}, {17ULL, 27ULL, 5ULL, 13ULL, 9ULL},  {17ULL, 8ULL, 5ULL, 37ULL, 19ULL},
   {20ULL, 41ULL, 5ULL, 14ULL, 6ULL},  {22ULL, 40ULL, 5ULL, 4ULL, 18ULL},  {22ULL, 19ULL, 5ULL, 14ULL, 19ULL},
   {22ULL, 41ULL, 5ULL, 16ULL, 6ULL},  {22ULL, 16ULL, 5ULL, 32ULL, 4ULL},  {26ULL, 9ULL, 5ULL, 11ULL, 14ULL},
   {26ULL, 19ULL, 5ULL, 29ULL, 3ULL},  {44ULL, 20ULL, 5ULL, 8ULL, 6ULL},   {44ULL, 31ULL, 5ULL, 22ULL, 14ULL},
   {53ULL, 8ULL, 5ULL, 23ULL, 17ULL},  {53ULL, 12ULL, 5ULL, 31ULL, 18ULL},
   /* table 6 (from elements 20 to 43 here) - safe for up to \f$10^{31}\f$ samples) */
   {10ULL, 5ULL, 29ULL, 23ULL, 8ULL},  {12ULL, 5ULL, 11ULL, 16ULL, 15ULL}, {17ULL, 5ULL, 16ULL, 6ULL, 7ULL},
   {17ULL, 5ULL, 19ULL, 16ULL, 14ULL}, {18ULL, 5ULL, 37ULL, 7ULL, 3ULL},   {19ULL, 5ULL, 31ULL, 15ULL, 13ULL},
   {20ULL, 5ULL, 11ULL, 13ULL, 6ULL},  {22ULL, 5ULL, 17ULL, 10ULL, 11ULL}, {23ULL, 5ULL, 37ULL, 13ULL, 7ULL},
   {24ULL, 5ULL, 7ULL, 16ULL, 8ULL},   {26ULL, 5ULL, 22ULL, 4ULL, 9ULL},   {26ULL, 5ULL, 26ULL, 13ULL, 12ULL},
   {26ULL, 5ULL, 31ULL, 14ULL, 13ULL}, {36ULL, 5ULL, 32ULL, 16ULL, 8ULL},  {36ULL, 5ULL, 32ULL, 21ULL, 8ULL},
   {39ULL, 5ULL, 19ULL, 6ULL, 8ULL},   {43ULL, 5ULL, 14ULL, 20ULL, 15ULL}, {44ULL, 5ULL, 14ULL, 15ULL, 15ULL},
   {44ULL, 5ULL, 29ULL, 6ULL, 13ULL},  {44ULL, 5ULL, 34ULL, 25ULL, 9ULL},  {45ULL, 5ULL, 16ULL, 21ULL, 8ULL},
   {51ULL, 5ULL, 28ULL, 3ULL, 12ULL},  {53ULL, 5ULL, 26ULL, 16ULL, 8ULL},  {54ULL, 5ULL, 28ULL, 13ULL, 3ULL}
};

/*! \brief  Four-element streams for L'ecuyer's combined LFSR (Tausworthe) generator */ 
uint64_t sTable543[106][4] = {
   /* table 5, first 4 lines of MathsComput(1999)p261 (from elements 0 to 3 here) - safe for up to \f$10^{21}\f$ samples */
   {30ULL, 23ULL, 17ULL, 18ULL}, {13ULL, 23ULL, 26ULL, 5ULL}, {17ULL, 38ULL, 23ULL, 24ULL}, {26ULL, 47ULL, 17ULL, 19ULL},
   /* table 5, last 2 lines (from elements 4 to 5 here) - safe for up to \f$10^{21}\f$ samples */
   {26ULL, 34ULL, 20ULL, 17ULL}, {29ULL, 38ULL, 28ULL, 18ULL},
   /* table 4 of (from elements 6 to 97 here) - safe for up to \f$10^{17}\f$ samples */
   {18ULL, 10ULL, 23ULL, 11ULL}, {26ULL, 10ULL, 13ULL, 11ULL}, {48ULL, 17ULL, 30ULL, 11ULL}, {27ULL, 20ULL, 9ULL, 11ULL},  
   {46ULL, 22ULL, 9ULL, 11ULL},  {23ULL, 29ULL, 24ULL, 11ULL}, {25ULL, 29ULL, 13ULL, 11ULL}, {34ULL, 29ULL, 9ULL, 11ULL},  
   {50ULL, 7ULL, 38ULL, 12ULL},  {15ULL, 8ULL, 19ULL, 12ULL},  {44ULL, 22ULL, 16ULL, 12ULL}, {6ULL, 23ULL, 29ULL, 12ULL},
   {16ULL, 5ULL, 22ULL, 13ULL},  {11ULL, 10ULL, 25ULL, 13ULL}, {18ULL, 11ULL, 40ULL, 13ULL}, {19ULL, 16ULL, 30ULL, 13ULL}, 
   {45ULL, 23ULL, 24ULL, 13ULL}, {17ULL, 7ULL, 9ULL, 14ULL},   {52ULL, 11ULL, 20ULL, 14ULL}, {52ULL, 22ULL, 30ULL, 14ULL}, 
   {25ULL, 23ULL, 26ULL, 14ULL}, {27ULL, 7ULL, 19ULL, 15ULL},  {25ULL, 11ULL, 13ULL, 15ULL}, {6ULL, 26ULL, 31ULL, 15ULL},
   {19ULL, 28ULL, 25ULL, 15ULL}, {38ULL, 28ULL, 37ULL, 15ULL}, {53ULL, 28ULL, 18ULL, 15ULL}, {50ULL, 29ULL, 32ULL, 15ULL}, 
   {17ULL, 32ULL, 41ULL, 15ULL}, {39ULL, 8ULL, 12ULL, 16ULL},  {53ULL, 13ULL, 33ULL, 16ULL}, {12ULL, 5ULL, 13ULL, 17ULL},  
   {16ULL, 5ULL, 11ULL, 17ULL},  {25ULL, 7ULL, 32ULL, 17ULL},  {54ULL, 10ULL, 36ULL, 17ULL}, {45ULL, 11ULL, 29ULL, 17ULL},
   {30ULL, 20ULL, 18ULL, 17ULL}, {39ULL, 20ULL, 43ULL, 17ULL}, {19ULL, 22ULL, 22ULL, 17ULL}, {50ULL, 23ULL, 25ULL, 17ULL}, 
   {11ULL, 26ULL, 19ULL, 17ULL}, {19ULL, 26ULL, 11ULL, 17ULL}, {13ULL, 29ULL, 40ULL, 17ULL}, {46ULL, 32ULL, 29ULL, 17ULL}, 
   {20ULL, 4ULL, 31ULL, 18ULL},  {5ULL, 10ULL, 33ULL, 18ULL},  {43ULL, 16ULL, 31ULL, 18ULL}, {38ULL, 23ULL, 37ULL, 18ULL},
   {46ULL, 25ULL, 39ULL, 18ULL}, {47ULL, 4ULL, 26ULL, 19ULL},  {33ULL, 7ULL, 27ULL, 19ULL},  {18ULL, 11ULL, 17ULL, 19ULL}, 
   {43ULL, 11ULL, 37ULL, 19ULL}, {5ULL, 14ULL, 13ULL, 19ULL},  {53ULL, 20ULL, 27ULL, 19ULL}, {24ULL, 25ULL, 25ULL, 19ULL}, 
   {30ULL, 25ULL, 27ULL, 19ULL}, {34ULL, 29ULL, 41ULL, 19ULL}, {18ULL, 5ULL, 36ULL, 20ULL},  {15ULL, 11ULL, 18ULL, 20ULL},
   {52ULL, 11ULL, 34ULL, 20ULL}, {5ULL, 22ULL, 10ULL, 20ULL},  {9ULL, 22ULL, 10ULL, 20ULL},  {16ULL, 23ULL, 38ULL, 20ULL}, 
   {17ULL, 23ULL, 26ULL, 20ULL}, {40ULL, 23ULL, 37ULL, 20ULL}, {46ULL, 23ULL, 5ULL, 20ULL},  {6ULL, 28ULL, 27ULL, 20ULL},  
   {25ULL, 28ULL, 33ULL, 20ULL}, {5ULL, 32ULL, 26ULL, 20ULL},  {13ULL, 7ULL, 37ULL, 21ULL},  {26ULL, 8ULL, 41ULL, 21ULL},
   {37ULL, 10ULL, 43ULL, 21ULL}, {38ULL, 10ULL, 11ULL, 21ULL}, {30ULL, 13ULL, 39ULL, 21ULL}, {38ULL, 16ULL, 43ULL, 21ULL}, 
   {9ULL, 17ULL, 32ULL, 21ULL},  {34ULL, 25ULL, 17ULL, 21ULL}, {38ULL, 26ULL, 41ULL, 21ULL}, {8ULL, 28ULL, 31ULL, 21ULL},  
   {19ULL, 29ULL, 12ULL, 21ULL}, {37ULL, 32ULL, 27ULL, 21ULL}, {27ULL, 8ULL, 5ULL, 22ULL},   {8ULL, 10ULL, 29ULL, 22ULL},
   {41ULL, 10ULL, 25ULL, 22ULL}, {50ULL, 13ULL, 4ULL, 22ULL},  {55ULL, 13ULL, 37ULL, 22ULL}, {50ULL, 17ULL, 36ULL, 22ULL}, 
   {39ULL, 26ULL, 29ULL, 22ULL}, {55ULL, 26ULL, 23ULL, 22ULL}, {13ULL, 28ULL, 16ULL, 22ULL}, {51ULL, 32ULL, 10ULL, 22ULL},
   /* table 3 (from elements 98 to 105 here) - safe for up to \f$10^{14}\f$ samples */
   {18ULL, 28ULL, 7ULL, 8ULL},   {26ULL, 20ULL, 11ULL, 7ULL},  {19ULL, 25ULL, 12ULL, 9ULL},  {18ULL, 31ULL, 13ULL, 6ULL},
   {18ULL, 22ULL, 16ULL, 6ULL},  {30ULL, 28ULL, 17ULL, 9ULL},  {17ULL, 28ULL, 18ULL, 6ULL},  {12ULL, 8ULL, 22ULL, 9ULL}
};

/* \brief q vectors for L'ecuyer's combined Tausworthe generators (table 7 and 6 (five elements)) */
uint64_t qTable76[2][5] = { {1ULL, 7ULL, 24ULL, 3ULL, 5ULL}, {1ULL, 24ULL, 3ULL, 5ULL, 3ULL} };
/* \brief k vectors for L'ecuyer's combined Tausworthe generators (tables 7 and 6 (five elements)) */
uint64_t kTable76[2][5] = { {63ULL, 57ULL, 55ULL, 52ULL, 47ULL}, {63ULL, 55ULL, 52ULL, 47ULL, 41ULL} };

/* \brief q vectors for L'ecuyer's combined Tausworthe generators (table 5[1...4], table 5[5,6], table 4 and table 3 */
uint64_t qTable543[4][4] = {
   {31ULL, 1ULL, 19ULL, 22ULL}, {31ULL, 11ULL, 19ULL, 22ULL}, {1ULL, 19ULL, 7ULL, 24ULL}, {31ULL, 19ULL, 24ULL, 21ULL}
};
/* \brief k vectors for L'ecuyer's combined Tausworthe generators (table 5[1...4], table 5[5,6], table 4 and table 3*/
uint64_t kTable543[4][4] = {
   {63ULL, 60ULL, 58ULL, 57ULL}, {63ULL, 60ULL, 58ULL, 57ULL}, {63ULL, 58ULL, 57ULL, 55ULL}, {63ULL, 58ULL, 55ULL, 47ULL}
};

uint64_t Cmask[28] = {
  0xfffffffff0000000ULL, 0xfffffffff8000000ULL, 0xfffffffffc000000ULL, 0xfffffffffe000000ULL, /* 36-39 */
  0xffffffffff000000ULL, 0xffffffffff800000ULL, 0xffffffffffc00000ULL, 0xffffffffffe00000ULL, /* 40-43 */
  0xfffffffffff00000ULL, 0xfffffffffff80000ULL, 0xfffffffffffc0000ULL, 0xfffffffffffe0000ULL, /* 44-47 */
  0xffffffffffff0000ULL, 0xffffffffffff8000ULL, 0xffffffffffffc000ULL, 0xffffffffffffe000ULL, /* 48-51 */
  0xfffffffffffff000ULL, 0xfffffffffffff800ULL, 0xfffffffffffffc00ULL, 0xfffffffffffffe00ULL, /* 52-55 */
  0xffffffffffffff00ULL, 0xffffffffffffff80ULL, 0xffffffffffffffc0ULL, 0xffffffffffffffe0ULL, /* 56-59 */
  0xfffffffffffffff0ULL, 0xfffffffffffffff8ULL, 0xfffffffffffffffcULL, 0xfffffffffffffffeULL  /* 60-63 */
};

/* http://www.math.uni-bielefeld.de/~sillke/ALGORITHMS/random 
 * some possible 16-bit constants k for which both k*2^16-1 and k*2^15-1 are prime */
uint32_t marsaglia_constants[81] = {
  18000, 18030, 18273, 18513, 18879, 19074, 19098, 19164, 19215, 19584, 19599, 19950, 20088, 20508, 20544, 20664, 20814,
  20970, 21153, 21243, 21423, 21723, 21954, 22125, 22188, 22293, 22860, 22938, 22965, 22974, 23109, 23124, 23163, 23208,
  23508, 23520, 23553, 23658, 23865, 24114, 24219, 24660, 24699, 24864, 24948, 25023, 25308, 25443, 26004, 26088, 26154,
  26550, 26679, 26838, 27183, 27258, 27753, 27795, 27810, 27834, 27960, 28320, 28380, 28689, 28710, 28794, 28854, 28959,
  28980, 29013, 29379, 29889, 30135, 30345, 30459, 30714, 30903, 30963, 31059, 31083, 36969
};

#endif
