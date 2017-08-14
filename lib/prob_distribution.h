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

/*! \file prob_distribution.h
 *  \brief Probability distribution functions and auxiliary mathematical functions from statistical package R 
 *
 * Code derived from the <a href="http://www.r-project.org/">R project for Statistical Computing</a> version 2.9.1, 
 * available under the GPL license. It might be possible to use directly the standalone mathematical library "Rmath.h" 
 * from the R project instead of our implementation. The advantage would be a library updated more often than mine, at 
 * the cost of delegating to the guenomu user the installation and maintenance of the extra libraries (like GSL, for 
 * instance). In Debian this library can be installed through the package "r-mathlib". 
 * The original R library checks for several built-in compiler functions (like log1p(e) for calculating log(e+1) ) but I
 * simply assume the compiler has none and reimplement them.
 * The CDFs always assume the lower tail (upper tails must use 1. - lower tail) or equivalent. 
 *
 * The code for the discrete sampling comes from the GNU Scientific Library version 1.14
 */

#ifndef _biomcmc_prob_distribution_h_
#define _biomcmc_prob_distribution_h_

#include "lowlevel.h"
#include "random_number.h"

typedef struct discrete_sample_struct* discrete_sample;

/* \brief sampling from a discrete distribution in O(1) time (from GSL based on Walker's algorithm) */
struct discrete_sample_struct
{
  size_t K; /* note that here we DO NOT use int, but size_t (which limits somehow the largest vector size, but practical) */
  size_t *A;
  double *F;
};

// reminder: gamma is in E[x] = a/b , but poisson is E[x] = lambda (since lambda is related to alpha, not beta) */

/*! \brief Ziheng Yang's gamma discretization of rates */
void biomcmc_discrete_gamma (double alpha, double beta, double *rate, int nrates);

/*! \brief pdf of discrete truncated exponential (d is discrete, m is maximum value) */
double biomcmc_dexp_dt (double d, double lambda, double m, bool log_p);
/*! \brief cdf of discrete truncated exponential (d is discrete, m is maximum value): calculates P(D <= d) */
double biomcmc_pexp_dt (double d, double lambda, double m, bool log_p);
/*! \brief quantile of discrete truncated exponential, that is, finds d s.t. P(D <= d) >= p */
double biomcmc_qexp_dt (double p, double lambda, double m, bool log_p);

/*! \brief gamma density */
double biomcmc_dgamma (double x, double alpha, double beta, bool log_p);
/*! \brief gamma quantile (inverse CDF) */
double biomcmc_qgamma (double p, double alpha, double beta, bool log_p);
/*! \brief computes the cummulative distribution function for the gamma distribution with shape parameter alpha and rate 
  parameter beta, s.t. \f$ E[X] = \alpha/\beta \f$. The same as the (lower) incomplete gamma function. */
double biomcmc_pgamma (double x, double alpha, double beta, bool log_p);

double biomcmc_dnorm (double x, double mu, double sigma, bool log_p);
double biomcmc_qnorm (double p, double mu, double sigma, bool log_p);
double biomcmc_pnorm (double x, double mu, double sigma, bool log_p);

double biomcmc_dlnorm (double x, double meanlog, double sdlog, bool log_p);
double biomcmc_qlnorm (double p, double meanlog, double sdlog, bool log_p);
double biomcmc_plnorm (double x, double meanlog, double sdlog, bool log_p);

double biomcmc_dpois (double x, double lambda, bool log_p);
/*The quantile function of the Poisson distribution. */
double biomcmc_qpois (double p, double lambda, bool log_p);
double biomcmc_ppois (double x, double lambda, bool log_p);

double biomcmc_rng_gamma (double alpha, double beta);
/*! \brief Returns a random number from a Normal distribution N(mu, sigma^2) using 52 bits of precision */
double biomcmc_rng_norm (double mu, double sigma);
double biomcmc_rng_lnorm (double meanlog, double sdlog);
double biomcmc_rng_pois (double mu);


/* Calculates log|gamma(x)| and returns the sign of the gamma function in sgn if *sgn is not NULL. */
double biomcmc_lgammafn (double x, int *sgn);
/*  This function computes the value of the gamma function. */
double biomcmc_gammafn (double x);

/*! \brief compute the relative error logarithm \f$ \log(1 + x)\f$ (C99 standard) */
double biomcmc_log1p (double x);
/*! \brief accurate calculation of \f$\log(1+x)-x\f$, particularly for small x  */
double biomcmc_log1pmx (double x);
/*! \brief compute \f$ \exp(x) - 1\f$ accurately also when x is close to zero, i.e. \f$|x| \ll 1 \f$ */
double biomcmc_expm1 (double x);

/* create the struct for sampling from an arbitrary empirical frequency vector */
discrete_sample new_discrete_sample_from_frequencies (double *prob, size_t size);
/* free discrete_sample_struct */
void del_discrete_sample (discrete_sample g);
/* sample index from empirical discrete distribution (between zero and g->K which is the original vector size) */
size_t biomcmc_rng_discrete (discrete_sample g);
/* density from empirical frequency vector (that is, original values from vector) */
double biomcmc_discrete_sample_pdf (discrete_sample g, size_t k);

/* Compute log (exp (logx) + exp (logy)) without overflow and without loss of accuracy. */
double biomcmc_logspace_add (double logx, double logy);
/* Compute log (exp (logx) - exp (logy)) without overflow and without loss of accuracy. */
double biomcmc_logspace_sub (double logx, double logy);
/*! \brief check if number is between minus infinity and plus infinity, or NaN */
bool biomcmc_isfinite(double x);

#endif
