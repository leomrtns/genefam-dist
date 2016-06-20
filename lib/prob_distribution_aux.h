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

/*! \file prob_distribution_aux.h
 *  \brief Auxiliary (low level) functions for prob_distribution.c
 */

#ifndef _biomcmc_prob_distribution_aux_h_
#define _biomcmc_prob_distribution_aux_h_

#include "prob_distribution.h"

#define R_LN2          0.693147180559945309417232121458  /* ln(2) */
#define R_PI           3.141592653589793238462643383280  /* pi */
#define R_2PI          6.283185307179586476925286766559  /* 2*pi */
#define R_EXP_M1       0.367879441171442321595523770161  /* exp(-1) = 1/e */
#define R_SQRT_32      5.656854249492380195206754896838  /* sqrt(32) */
#define R_1_SQRT_2PI   0.398942280401432677939946059934  /* 1/sqrt(2*pi) */
#define R_LN_SQRT_2PI  0.918938533204672741780329736406  /* log(sqrt(2*pi)) = log(2*pi)/2 */
#define R_LN_SQRT_PId2 0.225791352644727432363097614947  /* log(sqrt(pi/2)) */

const double pInf = 1./0., mInf = -1./0., NaN = 0./0.;
/* Scalefactor:= (2^32)^8 = 2^256 = 1.157921e+77 */
const double scalefactor = 115792089237316195423570985008687907853269984665640564039457584007913129639936.;
/* If |x| > |k| * M_cutoff,  then  log[ exp(-x) * k^x ]	 =~=  -x */ 
const double M_cutoff = R_LN2 * DBL_MAX_EXP / DBL_EPSILON;/*=3.196577e18*/

/* log gamma correction factor for x >= 10 so that log(gamma(x)) = .5*log(2*pi) + (x-.5)*log(x) -x + lgammacor(x) */  
double lgammacor (double x);
/* determines the number of terms for the double precision orthogonal series "dos" needed to insure
 * the error is no larger than "eta".  Ordinarily eta will be chosen to be one-tenth machine precision. */
int chebyshev_init (const double *dos, int nos, double eta);
/* evaluates the n-term Chebyshev series "a" at "x" */
double chebyshev_eval (double x, const double *a, const int n);
/* This function calculates the minimum and maximum legal bounds for x in biomcmc_gammafn(x). These are not the only 
 * bounds, but they are the only non-trivial ones to calculate. */
void gammalims (double *xmin, double *xmax);
/* Continued fraction for calculation of 1/i + x/(i+d) + x^2/(i+2*d) + x^3/(i+3*d) + ... = sum_{k=0}^Inf x^k/(i+k*d)
 * auxilary in biomcmc_log1pmx() and lgamma1p(); eps = relative tolerance */
double logcf (double x, double i, double d, double eps);
/* Compute  log(gamma(a+1))  accurately also for small a (0 < a < 0.5). */
double lgamma1p (double a);
/* wrapper over poisson density  (== dpois (x_P_1 - 1, lambda, log_p)) */
double dpois_wrap (double x_plus_1, double lambda, bool log_p);
/* lowlevel poisson density (assuming E[X] = lambda, which is the scale and not rate) */
double dpois_raw(double x, double lambda, bool log_p);
/* computes the log of the error term in Stirling's formula. */
double stirlerr(double n);
/* Evaluates the "deviance part" bd0(x,M) =  M * D0(x/M) = M*[ x/M * log(x/M) + 1 - (x/M) ] =  x * log(x/M) + M - x
 *  where M = E[X] = n*p (or = lambda), for   x, M > 0 */
double bd0(double x, double np);
/* auxiliary function for pgamma() */
double pgamma_smallx (double x, double alph, bool log_p);
/* sum of upper terms of series (auxiliar to pgamma()) */
double pd_upper_series (double x, double y, bool log_p);
/* sum of lower terms of series (auxiliar to pgamma()) */
double pd_lower_series (double lambda, double y);
/* Continued fraction for calculation of (i/d) + o(i/d) ( auxiliar to pgamma())*/
double pd_lower_cf (double i, double d);
/* Compute \f$\frac{dnorm (x, 0, 1, FALSE)}{pnorm (x, 0, 1, lower_tail = T, FALSE)} \f$ accurately */
double dpnorm (double x, double lp);
/* Asymptotic expansion to calculate the probability that Poisson variate has value <= x. */
double ppois_asymp (double x, double lambda, bool log_p);
/* lowlevel pgamma() calculation (without variable checking etc) */
double pgamma_raw (double x, double alph, bool log_p);
/* chi-squared approximation to qgamma, g = log Gamma (nu/2) and tol = EPS1 */
double qchisq_appr (double p, double nu, double g, bool log_p, double tol);
/* returns both cummulative points from the normal */
void pnorm_both(double x, double *cum, double *ccum, int i_tail, bool log_p);
/* find quantile for poisson distribution */
double do_poisson_search (double y, double *z, double p, double lambda, double incr);

#endif
