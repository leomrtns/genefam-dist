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

/* Algorithms from R version 2.9.1 
 *
 * prob_distrib_aux.c -link_to->  prob_distrib_aux.h -link_to-> prob_distrib.h 
 * compiled into one object prob_distrib.o */
#include "prob_distribution_aux.c"

typedef struct dsample_stack_struct* dsample_stack;

/* \brief stack structure used in preprocessing of Walker's algorithm (sample from discrete distrib in O(1) time) from GSL */
struct dsample_stack_struct {
  size_t N;  /* max number of elts on stack */
  size_t *v; /* array of values on the stack */
  size_t i;  /* index of top of stack */
};

dsample_stack new_dsample_stack (size_t vector_size);
void         free_dsample_stack (dsample_stack s);


void
biomcmc_discrete_gamma (double alpha, double beta, double *rate, int nrates)
{
  if (nrates == 1) { rate[0] = alpha/beta; return; } /* gross, ill-defined and imoral model */
  else {
    int i;
    double K = nrates, scale = (alpha * K)/beta, prob[nrates-1];

    for (i = 0; i < nrates - 1; i++) prob[i] = biomcmc_qgamma ((double)(i+1)/K, alpha, beta, false);
    for (i = 0; i < nrates - 1; i++) prob[i] = biomcmc_pgamma (prob[i], alpha + 1., beta, false);

    rate[0] = prob[0] * scale;
    rate[nrates-1] = (1. - prob[nrates-2]) * scale;
    for (i = 1; i < nrates - 1; i++) rate[i] = (prob[i] - prob[i-1]) * scale;
    return;
  }
}

/*! \brief pdf of discrete truncated exponential (d is discrete, m is maximum value; E[X]=1/lambda) */
double
biomcmc_dexp_dt (double d, double lambda, double m, bool log_p)
{
  double prob_sum;
  if ( (d==NaN) || (lambda==NaN) || (m==NaN)|| (d < 0.) || (lambda <= 0.) || (m < 1.)) return NaN;
  if (d > m)  return log_p ? mInf : 0.;
  prob_sum = biomcmc_expm1 (-lambda * (m + 1.)) / biomcmc_expm1 (-lambda); /* expm1(x) = (exp(x) - 1) when x ~ 0 */

  return log_p ? ((-lambda * d) - log (prob_sum)) : (exp(-lambda * d)/prob_sum);
}

/*! \brief cdf of discrete truncated exponential (d is discrete, m is maximum value): calculates P(D <= d) */
double
biomcmc_pexp_dt (double d, double lambda, double m, bool log_p)
{
  if ( (d==NaN) || (lambda==NaN) || (m==NaN)|| (d < 0.) || (lambda <= 0.) || (m < 1.)) return NaN;
  if (d >= m) return log_p ? 0. : 1.;

  return log_p ? 
  (log (- biomcmc_expm1 (-lambda * (d + 1.))) - log (- biomcmc_expm1 (-lambda * (m + 1.)))) :
  (biomcmc_expm1 (-lambda * (d + 1.)) / biomcmc_expm1 (-lambda * (m + 1.)));
}

/*! \brief quantile of discrete truncated exponential, that is, finds d s.t. P(D <= d) >= p */
double
biomcmc_qexp_dt (double p, double lambda, double m, bool log_p)
{
  double x;
  if ( (p==NaN) || (lambda==NaN) || (m==NaN)|| (lambda <= 0.) || (m < 1.)) return NaN;

  if (log_p) {         if (p > 0.) return NaN; if (p == 0.) return m; if (p == mInf) return 0.; }
  else { if ((p < 0.) || (p > 1.)) return NaN; if (p == 1.) return m;   if (p == 0.) return 0.; }
  if (log_p) p = exp (p);

  x = biomcmc_log1p (p * biomcmc_expm1 (-lambda * (m + 1.))); /* log1p(x) = log (1+x) */ 

  return ceil ( - x/lambda - 1.);
}

double 
biomcmc_dgamma (double x, double alpha, double beta, bool log_p)
{
  double pr;

  if ( (x==NaN) || (alpha==NaN) || (beta==NaN)|| (alpha < 0.) || (beta <= 0.)) return NaN;

  if (x < 0.) return log_p ? mInf : 0.;
  if (alpha == 0.) return (x == 0.)? pInf : (log_p ? mInf : 0.);
  if (x == 0.) {
    if (alpha < 1.) return pInf; 
    if (alpha > 1.) return log_p ? mInf : 0.;
    return log_p ? log(beta) : beta; /* alpha == 1. */
  }

  if (alpha < 1) {
    pr = dpois_raw (alpha, x * beta, log_p);
    return log_p ? (pr + log(alpha/x)) : (pr * alpha/x);
  }
  pr = dpois_raw (alpha-1, x * beta, log_p);
  return log_p ? (pr + log(beta)) : (pr * beta);
}

/* Compute the quantile function of the gamma distribution. This function is based on the Applied Statistics 
 * Algorithm AS 91 ("ppchi2") and via pgamma(.) AS 239.
 * R core improvements: lower_tail, log_p; non-trivial result for p outside [0.000002, 0.999998]; 
 * p ~ 1 no longer gives +Inf; final Newton step(s)
 * Best, D. J. and D. E. Roberts (1975). Percentage Points of the Chi-Squared Distribution.
 * Applied Statistics 24, page 385.  */
double 
biomcmc_qgamma (double p, double alpha, double beta, bool log_p)
{
  static const int MAXIT = 1000;
  static const double EPS1 = 1e-2, EPS2 = 5e-7, EPS_N = 1e-15; /* final precision after each step */
  static const double i420  = 1./ 420., i2520 = 1./ 2520., i5040 = 1./ 5040, pMIN = 1e-100, pMAX = (1-1e-14);
  double p_, a, b, c, g, ch, ch0, p1, p2, q, s1, s2, s3, s4, s5, s6, t, x;
  int i, max_it_Newton = 1;

  if ( (p==NaN) || (alpha==NaN) || (beta==NaN)|| (alpha < 0.) || (beta <= 0.)) return NaN;

  if (log_p) { if (p > 0.) return NaN; if (!p) return pInf; if (p == mInf) return 0.; }
  else { if ((p < 0.) || (p > 1.)) return NaN; if (!p) return 0.; if (p==1.) return pInf; }

  if (alpha == 0) return 0.; /* all mass at 0 : */

  p_ = log_p ? exp (p) : p;
  g = biomcmc_lgammafn (alpha, NULL);/* log Gamma(v/2) */

  /* Phase I : Starting Approximation */
  ch = qchisq_appr (p, 2*alpha, g, log_p, EPS1); /* 2*alpha = nu = 'df' */
  if(!biomcmc_isfinite (ch)) { max_it_Newton = 0; goto END; }
  if(ch < EPS2) { max_it_Newton = 20; goto END; } /* Corrected according to AS 91; MM, May 25, 1999 */
  if((p_ > pMAX) || (p_ < pMIN)) { max_it_Newton = 20; goto END; } /* did return ML_POSINF or 0.;	much better: */

  /* Phase II: Iteration *	Call pgamma() [AS 239]	and calculate seven term taylor series */
  c = alpha - 1.;
  s6 = (120 + c * (346+127*c)) * i5040; /* used below, is "const" */

  ch0 = ch;/* save initial approx. */
  for (i=1; i <= MAXIT; i++) {
    q = ch;
    p1 = 0.5 * ch;
    p2 = p_ - pgamma_raw (p1, alpha, false);
    if(!biomcmc_isfinite (p2) || (ch <= 0)) { ch = ch0; max_it_Newton = 27; goto END; }

    t = p2 * exp (alpha*R_LN2 + g + p1 - c*log(ch));
    b = t/ch;
    a = 0.5*t - b*c;
    s1 = (210+ a*(140+a*(105+a*(84+a*(70+60*a))))) * i420;
    s2 = (420+ a*(735+a*(966+a*(1141+1278*a)))) * i2520;
    s3 = (210+ a*(462+a*(707+932*a))) * i2520;
    s4 = (252+ a*(672+1182*a) + c*(294+a*(889+1740*a))) * i5040;
    s5 = (84+2264*a + c*(1175+606*a)) * i2520;

    ch += t*(1+0.5*t*s1-b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))));
    if(fabs(q - ch) < (EPS2 * ch)) goto END;
    if(fabs(q - ch) > (0.1 * ch)) { if(ch < q) ch = 0.9 * q; else ch = 1.1 * q;} /* diverging? also forces ch > 0 */
  }

END:
  /* PR# 2214 :	 From: Morten Welinder <terra@diku.dk>, Fri, 25 Oct 2002 16:50 To: R-bugs@biostat.ku.dk
   * With a final Newton step, double accuracy, e.g. for (p= 7e-4; nu= 0.9)
   * Improved (MM): - only if rel.Err > EPS_N (= 1e-15); optionally *iterate* Newton */
  x = 0.5*ch/beta;

  if (max_it_Newton) {
    if (!log_p) { p = log(p); log_p = true; } /* always use log scale */
    p_ = biomcmc_pgamma (x, alpha, beta, log_p);

    for(i = 1; i <= max_it_Newton; i++) {
      p1 = p_ - p;
      if(fabs(p1) < fabs(EPS_N * p)) break;
      if((g = biomcmc_dgamma (x, alpha, beta, log_p)) == mInf) break;
      /* delta x = f(x)/f'(x);
       * if(log_p) f(x) := log P(x) - p; f'(x) = d/dx log P(x) = P' / P
       * ==> f(x)/f'(x) = f*P / P' = f*exp(p_) / P' (since p_ = log P(x)) */
      t = x - (p1 * exp(p_ - g));/* = x - delta x" */
      p_ = biomcmc_pgamma (t, alpha, beta, log_p);
      if (fabs(p_ - p) > fabs(p1) || ((i > 1) && (fabs(p_ - p) == fabs(p1)))) break; /* against flip-flop */

      if(t > (1.1*x)) t = 1.1*x; /* control step length: this could have started at the initial approximation */
      else if(t < (0.9*x)) t = 0.9*x;
      x = t;
    }
  }

  return x;
}

double 
biomcmc_pgamma (double x, double alpha, double beta, bool log_p)
{
  if ( (x==NaN) || (alpha==NaN) || (beta==NaN)|| (alpha < 0.) || (beta <= 0.)) return NaN;
  x *= beta;
  if (x == NaN) return NaN;
  if (!alpha) return (x < 0.)?(log_p? mInf: 0.) : (log_p? 0.: 1.); /* limit case; useful e.g. in pnchisq() */
  return pgamma_raw (x, alpha, log_p);
}

double
biomcmc_dnorm (double x, double mu, double sigma, bool log_p)
{
  if ((x==NaN) || (mu==NaN) || (sigma==NaN) || ((!biomcmc_isfinite (x)) && (mu == x)) || (sigma < 0.)) return NaN;
  if ((!biomcmc_isfinite (sigma)) || ((sigma == 0.) && (x != mu))) return log_p ? mInf : 0.;
  if ((sigma == 0.) && (mu == x)) return pInf; 

  x = (x - mu) / sigma;
  if(!biomcmc_isfinite (x)) return log_p ? mInf : 0.;

  return (log_p? -(R_LN_SQRT_2PI +  0.5 * x * x + log(sigma)) : (R_1_SQRT_2PI * exp(-0.5 * x * x)  /   sigma));
}

/*	Compute the quantile function for the normal distribution. For small to moderate probabilities, algorithm 
 *	referenced below is used to obtain an initial approximation which is polished with a final Newton step.
 *	For very large arguments, an algorithm of Wichura is used.
 *	Beasley, J. D. and S. G. Springer (1977).
 *	Algorithm AS 111: The percentage points of the normal distribution, Applied Statistics, 26, 118-121.
 *  Wichura, M.J. (1988).
 *  Algorithm AS 241: The Percentage Points of the Normal Distribution. Applied Statistics, 37, 477-484. */
double 
biomcmc_qnorm (double p, double mu, double sigma, bool log_p)
{
  double p_, q, r, val;

  if ((p==NaN) || (mu==NaN) || (sigma==NaN) || (sigma < 0.)) return NaN;
  if(sigma == 0)	return mu;

  if (log_p) {         if (p > 0.) return NaN; if (p == 0.) return pInf; if (p == mInf) return mInf; }
  else { if ((p < 0.) || (p > 1.)) return NaN; if (p == 1.) return pInf;   if (p == 0.) return mInf; }

  if (log_p) p_ = exp (p);
  else p_ = p;
  q = p_ - 0.5;

  /* use AS 241: double ppnd16_(double *p, long *ifault) ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3
     Produces the normal deviate Z corresponding to a given lower tail area of P; Z is accurate to about 1 part 
     in 10**16. */
  if (fabs(q) <= .425) {/* 0.075 <= p <= 0.925 */
    r = .180625 - q * q;
    val =
    q * (((((((r * 2509.0809287301226727 + 33430.575583588128105) * r + 67265.770927008700853) * r + 
             45921.953931549871457) * r + 13731.693765509461125) * r + 1971.5909503065514427) * r + 
          133.14166789178437745) * r + 3.387132872796366608)
    / (((((((r * 5226.495278852854561 + 28729.085735721942674) * r + 39307.89580009271061) * r +
           21213.794301586595867) * r + 5394.1960214247511077) * r + 687.1870074920579083) * r + 
        42.313330701600911252) * r + 1.);
  }
  else { /* closer than 0.075 from {0,1} boundary */
    /* r = min(p, 1-p) < 0.075 */
    if (q > 0) r = log_p ? (- biomcmc_expm1 (p)) : (0.5 - p + 0.5); /* 0.5 + 0.5 -p = 1 - p */
    else r = p_; 

    r = sqrt (- ((log_p && (q <= 0)) ? p : log (r))); /* r = sqrt(-log(r))  <==>  min(p, 1-p) = exp( - r^2 ) */

    if (r <= 5.) { /* <==> min(p,1-p) >= exp(-25) ~= 1.3888e-11 */
      r += -1.6;
      val = (((((((r * 7.7454501427834140764e-4 + .0227238449892691845833) * r + .24178072517745061177) *
                 r + 1.27045825245236838258) * r + 3.64784832476320460504) * r + 5.7694972214606914055) *
              r + 4.6303378461565452959) * r + 1.42343711074968357734)
      / (((((((r * 1.05075007164441684324e-9 + 5.475938084995344946e-4) * r + .0151986665636164571966) * r +
             .14810397642748007459) * r + .68976733498510000455) * r + 1.6763848301838038494) * r + 
          2.05319162663775882187) * r + 1.);
    }
    else { /* very close to  0 or 1 */
      r += -5.;
      val = (((((((r * 2.01033439929228813265e-7 + 2.71155556874348757815e-5) * r + .0012426609473880784386) * r + 
                 .026532189526576123093) * r + .29656057182850489123) * r + 1.7848265399172913358) * r + 
              5.4637849111641143699) * r + 6.6579046435011037772)
      / (((((((r * 2.04426310338993978564e-15 + 1.4215117583164458887e-7)* r + 1.8463183175100546818e-5) * r +
             7.868691311456132591e-4) * r + .0148753612908506148525) * r + .13692988092273580531) * r +
          .59983220655588793769) * r + 1.);
    }

    if(q < 0.0) val = -val;
    /* return (q >= 0.)? r : -r ;*/
  }
  return mu + sigma * val;
}

/* The main computation evaluates near-minimax approximations derived from those in "Rational Chebyshev 
 * approximations for the error function" by W. J. Cody, Math. Comp., 1969, 631-637.  This transportable program 
 * uses rational functions that theoretically approximate the normal distribution function to at least 18
 * significant decimal digits.  The accuracy achieved depends on the arithmetic system, the compiler, the intrinsic 
 * functions, and proper selection of the machine-dependent constants.
 * REFERENCE: Cody, W. D. (1993). ALGORITHM 715: SPECFUN - A Portable FORTRAN Package of Special Function Routines 
 * and Test Drivers". ACM Transactions on Mathematical Software. 19, 22-32.*/
double 
biomcmc_pnorm (double x, double mu, double sigma, bool log_p)
{
  double p, cp;

  if ((x==NaN) || (mu==NaN) || (sigma==NaN) || ((!biomcmc_isfinite (x)) && (mu == x)) || (sigma < 0.)) return NaN;
  if (sigma == 0.) return (x < mu) ? (log_p? mInf : 0.) : (log_p ? 0. : 1.);

  p = (x - mu) / sigma;
  if (!biomcmc_isfinite (p)) return (x < mu) ? (log_p? mInf : 0.) : (log_p ? 0. : 1.);

  x = p;
  pnorm_both(x, &p, &cp, 0, log_p); /* the zero stands for lower tail (0,1,2 -> lower, upper, both) */
  return p;
}

double 
biomcmc_dlnorm (double x, double meanlog, double sdlog, bool log_p)
{
  double y;

  if ((x == NaN) || (meanlog == NaN) || (sdlog == NaN) || (sdlog <= 0.)) return NaN;
  if (x <= 0.) return log_p ? mInf : 0.;

  y = (log (x) - meanlog) / sdlog;
  return (log_p ?
          -(R_LN_SQRT_2PI   + 0.5 * y * y + log (x * sdlog)) :
          R_1_SQRT_2PI * exp (-0.5 * y * y)  /  (x * sdlog));
}

double 
biomcmc_qlnorm (double p, double meanlog, double sdlog, bool log_p)
{
  if ((p == NaN) || (meanlog == NaN) || (sdlog == NaN)) return NaN;

  if (log_p) {         if (p > 0.) return NaN; if (p == 0.) return pInf; if (p == mInf) return 0.; }
  else { if ((p < 0.) || (p > 1.)) return NaN; if (p == 1.) return pInf;   if (p == 0.) return 0.; }

  return exp (biomcmc_qnorm (p, meanlog, sdlog, log_p));
}

double 
biomcmc_plnorm (double x, double meanlog, double sdlog, bool log_p)
{
  if ((x == NaN) || (meanlog == NaN) || (sdlog == NaN) || (sdlog <= 0.)) return NaN;
  if (x <= 0.) return log_p ? mInf : 0.;
  return biomcmc_pnorm (log (x), meanlog, sdlog, log_p);
}

double 
biomcmc_dpois (double x, double lambda, bool log_p)
{ /* E[X] = lambda (remember that lambda is related to "alpha" of gamma distrib, NOT "beta" */ 
  if ((x == NaN) || (lambda== NaN) || (lambda < 0.)) return NaN;
  if ((x < 0) || (!biomcmc_isfinite (x)) || (fabs (x - floor (x + 0.5)) > 1e-7)) return log_p ? mInf : 0.;
  x = floor (x + 0.5);
  return dpois_raw (x, lambda, log_p);
}

/*The quantile function of the Poisson distribution. Uses the Cornish-Fisher Expansion to include a skewness
 *correction to a normal approximation.  This gives an initial value which never seems to be off by more than
 *1 or 2. A search is then conducted of values close to this initial start point [ do_poisson_search() ] . */
double 
biomcmc_qpois (double p, double lambda, bool log_p)
{
  double mu, sigma, gamma, z, y;
  if ((p == NaN) || (lambda== NaN) || (lambda < 0.) || (!biomcmc_isfinite (lambda))) return NaN;

  if (log_p) {         if (p > 0.) return NaN; if (p == 0.) return pInf; if (p == mInf) return 0.; }
  else { if ((p < 0.) || (p > 1.)) return NaN; if (p == 1.) return pInf;   if (p == 0.) return 0.; }

  if (lambda == 0) return 0.;

  mu = lambda;
  sigma = sqrt(lambda);
  gamma = 1.0/sigma; /* gamma = sigma; PR#8058 should be kurtosis which is mu^-0.5 */

  /* Note : "same" code in qpois.c, qbinom.c, qnbinom.c - This is far from optimal [cancellation for p ~= 1, etc] */
  if (log_p) { /* need check again (cancellation!): */
    p = exp (p); 
    if (p == 0.) return 0.;
    if (p == 1.) return pInf;
  }
  if (p + 1.01 * DBL_EPSILON >= 1.) return pInf; /* R-mathlib hack, check for better options... */

  /* y = approx.value (Cornish-Fisher expansion) */
  z = biomcmc_qnorm (p, 0., 1., false);
  y = floor (mu + sigma * (z + gamma * (z*z - 1.) / 6.) + 0.5);
  z = biomcmc_ppois (y, lambda, true);
  p *= 1 - 64 * DBL_EPSILON; /* fuzz to ensure left continuity; 1 - 1e-7 may lose too much */

  /* If the mean is not too large a simple search is OK */
  if (lambda < 1e4) return do_poisson_search (y, &z, p, lambda, 1);
   {/* Otherwise be a bit cleverer in the search */
    double incr = floor (y * 0.001), oldincr, new_incr;
    do {
      oldincr = incr;
      y = do_poisson_search (y, &z, p, lambda, incr);
      new_incr = floor (incr / 100.);
      incr = MAX (1, new_incr);
    } while ((oldincr > 1) && (incr > (lambda * 1e-15)));
    return y;
   }
}

double 
biomcmc_ppois (double x, double lambda, bool log_p)
{
  if ((x == NaN) || (lambda== NaN) || (lambda < 0.)) return NaN;
  if (x < 0)               return log_p ? mInf : 0.;
  if (lambda == 0.)        return log_p ? 0. : 1.;
  if (!biomcmc_isfinite (x)) return log_p ? 0. : 1.;
  x = floor (x + 1e-7);

  return 1. - biomcmc_pgamma (lambda, x + 1, 1., log_p); /* upper.tail of gamma = 1. - lower.tail */
}

/* - Shape parameter a >= 1.  Algorithm GD in: Ahrens, J.H. and Dieter, U. (1982). 
 * Generating gamma variates by a modified rejection technique. Comm. ACM, 25, 47-54. 
 * - Shape parameter 0 < a < 1. Algorithm GS in: Ahrens, J.H. and Dieter, U. (1974). 
 * Computer methods for sampling from gamma, beta, poisson and binomial distributions. Computing, 12, 223-246. */
double 
biomcmc_rng_gamma (double alpha, double beta)
{
  /* Coefficients q[k] - for q0 = sum(q[k]*a^(-k))
   * Coefficients a[k] - for q = q0+(t*t/2)*sum(a[k]*v^k)
   * Coefficients e[k] - for exp(q)-1 = sum(e[k]*q^k) */
  static const double q[7] = { 0.04166669, 0.02083148, 0.00801191, 0.00144121, -7.388e-5,  2.4511e-4, 2.424e-4};
  static const double a[7] = { 0.3333333,   -0.250003,  0.2000062, -0.1662921, 0.1423657, -0.1367177, 0.1233795};
  static double aa = 0., aaa = 0.;
  static double s, s2, d;     /* no. 1 (step 1) */
  static double q0, b, si, c; /* no. 2 (step 4) */
  double e, p, q00, r, t, u, v, w, x, ret_val, exp_rand;

  if ( !biomcmc_isfinite (alpha) || !biomcmc_isfinite (beta) || (alpha < 0.) || (beta <= 0.)) return NaN;

  if (alpha < 1.) { /* GS algorithm for parameters a < 1 */
    if(alpha == 0.) return 0.;
    e = 1.0 + R_EXP_M1 * alpha;
    while (1) {
      p = e * biomcmc_rng_unif_pos ();
      exp_rand = - log (biomcmc_rng_unif_pos ());
      if (p >= 1.0) {
        x = -log ((e - p) / alpha);
        if (exp_rand >= (1.0 - alpha) * log (x)) break;
      }
      else {
        x = exp (log (p) / alpha);
        if (exp_rand >= x) break;
      }
    }
    return x / beta;
  }

  /* --- a >= 1 : GD algorithm --- */
  /* Step 1: Recalculations of s2, s, d if alpha has changed */
  if (alpha != aa) {
    aa = alpha;
    s2 = alpha - 0.5;
    s = sqrt (s2);
    d = R_SQRT_32 - s * 12.0;
  }
  /* Step 2: t = standard normal deviate, x = (s,1/2) -normal deviate. */

  /* immediate acceptance (i) */
  t = biomcmc_rng_snorm ();
  x = s + 0.5 * t;
  ret_val = x * x;
  if (t >= 0.) return ret_val / beta;
  /* Step 3: u = 0,1 - uniform sample. squeeze acceptance (s) */
  u = biomcmc_rng_unif ();
  if (d * u <= t * t * t) return ret_val/beta;
  /* Step 4: recalculations of q0, b, si, c if necessary */
  if (alpha != aaa) {
    aaa = alpha;
    r = 1.0 / alpha;
    q0 = ((((((q[6] * r + q[5]) * r + q[4]) * r + q[3]) * r + q[2]) * r + q[1]) * r + q[0]) * r;

    /* Approximation depending on size of parameter alpha; */
    /* The constants in the expressions for b, si and c were established by numerical experiments */
    if (alpha <= 3.686) { b = 0.463 + s + 0.178 * s2; si = 1.235; c = 0.195 / s - 0.079 + 0.16 * s; } 
    else if (alpha <= 13.022) { b = 1.654 + 0.0076 * s2; si = 1.68 / s + 0.275; c = 0.062 / s + 0.024; } 
    else { b = 1.77; si = 0.75; c = 0.1515 / s; }
  }
  /* Step 5: no quotient test if x not positive */
  if (x > 0.0) {
    /* Step 6: calculation of v and quotient q */
    v = t / (s + s);
    if (fabs (v) <= 0.25)
      q00 = q0 + 0.5 * t * t * ((((((a[6] * v + a[5]) * v + a[4]) * v + a[3]) * v + a[2]) * v + a[1]) * v + a[0]) * v;
    else
      q00 = q0 - s * t + 0.25 * t * t + (s2 + s2) * log(1.0 + v);
    /* Step 7: quotient acceptance (q) */
    if (log(1.0 - u) <= q00) return ret_val / beta;
  }

  while (1) {
    /* Step 8: e = standard exponential deviate
     *	       u =  0,1 -uniform deviate
     *	       t = (b,si)-double exponential (laplace) sample */
    u = biomcmc_rng_unif_pos ();
    e = - log (biomcmc_rng_unif_pos ());
    u = u + u - 1.0;
    if (u < 0.0) t = b - si * e;
    else t = b + si * e;
    /* Step	 9:  rejection if t < tau(1) = -0.71874483771719 */
    if (t >= -0.71874483771719) {
      /* Step 10:	 calculation of v and quotient q */
      v = t / (s + s);
      if (fabs (v) <= 0.25)
        q00 = q0 + 0.5 * t * t * ((((((a[6] * v + a[5]) * v + a[4]) * v + a[3]) * v + a[2]) * v + a[1]) * v + a[0]) * v;
      else
        q00 = q0 - s * t + 0.25 * t * t + (s2 + s2) * log (1.0 + v);
      /* Step 11:	 hat acceptance (h) (if q not positive go to step 8) */
      if (q00 > 0.0) {
        w = biomcmc_expm1 (q00);
        /* if t is rejected sample again at step 8 */
        if ((c * fabs (u)) <= (w * exp (e - 0.5 * t * t))) break;
      }
    }
  } /* repeat .. until  `t' is accepted */
  x = s + 0.5 * t;
  return x * x / beta;
}

double
biomcmc_rng_norm (double mu, double sigma)
{
  if ((mu == NaN) || (!biomcmc_isfinite (sigma)) || (sigma < 0.)) return NaN;
  if ((sigma == 0.) || (!biomcmc_isfinite (mu))) return mu;
  return mu + sigma * biomcmc_rng_snorm ();
}

double
biomcmc_rng_lnorm (double meanlog, double sdlog)
{
  if ((meanlog == NaN) || (!biomcmc_isfinite (sdlog)) || (sdlog < 0.)) return NaN;
  if (sdlog == 0.) return exp (meanlog);
  if (!biomcmc_isfinite (meanlog)) return meanlog;
  return exp (meanlog + sdlog * biomcmc_rng_snorm ());
}

/* Random variates from the Poisson distribution.
 * Ahrens, J.H. and Dieter, U. (1982). Computer generation of Poisson deviates from modified normal distributions.
 * ACM Trans. Math. Software 8, 163-179. */
double 
biomcmc_rng_pois (double mu)
{
  static const double a[8] = {-0.5, 0.3333333, -0.2500068, 0.2000118, -0.1661269, 0.1421878, -0.1384794, 0.1250060};
  static const double one_7 = 0.1428571428571428571, one_12 = 0.0833333333333333333, one_24 = 0.0416666666666666667;
  static const double fact[10] = {1., 1., 2., 6., 24., 120., 720., 5040., 40320., 362880.};/*Factorial Table (0:9)! */
  static int l, m;
  static double b1, b2, c, c0, c1, c2, c3;
  static double pp[36], p0, p, q, s, d, omega;
  static double big_l;/* integer "w/o overflow" */
  static double muprev = 0., muprev2 = 0.;/*, muold	 = 0.*/
  double del, difmuk= 0., E= 0., fk= 0., fx, fy, g, px, py, t, u= 0., v, x;
  double pois = -1.;
  int k; 
  bool kflag, big_mu, new_big_mu = false;

  if ((mu < 0.) || !biomcmc_isfinite (mu)) return NaN;

  big_mu = (mu >= 10.);
  if(big_mu) new_big_mu = false;

  if (!(big_mu && (mu == muprev))) {/* maybe compute new persistent par.s */
    if (big_mu) {
      new_big_mu = true;
      /* Case A. (recalculation of s,d,l	because mu has changed):
       * The poisson probabilities pk exceed the discrete normal probabilities fk whenever k >= m(mu). */
      muprev = mu;
      s = sqrt (mu);
      d = 6. * mu * mu;
      big_l = floor (mu - 1.1484); /* = an upper bound to m(mu) for all mu >= 10.*/
    }
    else { /* Small mu ( < 10) -- not using normal approx. */
      /* Case B. (start new table and calculate p0 if necessary) */
      /*muprev = 0.;-* such that next time, mu != muprev ..*/
      if (mu != muprev) {
        muprev = mu;
        m = MAX (1, (int) mu);
        l = 0; /* pp[] is already ok up to pp[l] */
        q = p0 = p = exp (-mu);
      }

      while (1) {
        /* Step U. uniform sample for inversion method */
        u = biomcmc_rng_unif ();
        if (u <= p0) return 0.;
        /* Step T. table comparison until the end pp[l] of the pp-table of cumulative poisson probabilities
           (0.458 > ~= pp[9](= 0.45792971447) for mu=10 ) */
        if (l != 0) {
          for (k = ((u <= 0.458) ? 1 : MIN (l, m));  k <= l; k++) if (u <= pp[k]) return (double) k;
          if (l == 35) continue; /* u > pp[35] */
        }
        /* Step C. creation of new poisson; probabilities p[l..] and their cumulatives q =: pp[k] */
        l++;
        for (k = l; k <= 35; k++) {
          p *= mu / k;
          q += p;
          pp[k] = q;
          if (u <= q) {
            l = k;
            return (double) k;
          }
        }
        l = 35;
      } // while (1)
    }   // mu < 10 
  } /* end {initialize persistent vars} */

  /* Only if mu >= 10 */
  /* Step N. normal sample */
  g = mu + s * biomcmc_rng_snorm ();
  if (g >= 0.) {
    pois = floor (g);
    /* Step I. immediate acceptance if pois is large enough */
    if (pois >= big_l) return pois;
    /* Step S. squeeze acceptance */
    fk = pois;
    difmuk = mu - fk;
    u = biomcmc_rng_unif ();
    if ((d * u) >= (difmuk * difmuk * difmuk)) return pois;
  }

  /* Step P. preparations for steps Q and H. (recalculations of parameters if necessary) */
  if (new_big_mu || (mu != muprev2)) {
    /* Careful! muprev2 is not always == muprev because one might have exited in step I or S */
    muprev2 = mu;
    omega = R_1_SQRT_2PI / s;
    /* The quantities b1, b2, c3, c2, c1, c0 are for the Hermite approximations to the discrete normal probs fk. */

    b1 = one_24 / mu;
    b2 = 0.3 * b1 * b1;
    c3 = one_7 * b1 * b2;
    c2 = b2 - 15. * c3;
    c1 = b1 - 6. * b2 + 45. * c3;
    c0 = 1. - b1 + 3. * b2 - 15. * c3;
    c = 0.1069 / mu; /* guarantees majorization by the 'hat'-function. */
  }
  /* 'Subroutine' F is called (kflag=0 for correct return) */
  if (g >= 0.) { kflag = false; goto Step_F; }

  while (1) {
    /* Step E. Exponential Sample */
    E =  - log (biomcmc_rng_unif_pos ());

    /*  sample t from the laplace 'hat' (if t <= -0.6744 then pk < fk for all mu >= 10.) */
    u = 2 * biomcmc_rng_unif () - 1.;
    t = (u >= 0) ? (1.8 + fabs (E)) : (1.8 - fabs (E));
    if (t > -0.6744) {
      pois = floor (mu + s * t);
      fk = pois;
      difmuk = mu - fk;
      /* 'subroutine' F is called (kflag=1 for correct return) */
      kflag = true;

Step_F: /* 'subroutine' F : calculation of px,py,fx,fy. */
      if (pois < 10) { /* use factorials from table fact[] */
        px = -mu;
        py = pow(mu, pois) / fact [(int) pois];
      }
      else {
        /* Case pois >= 10 uses polynomial approximation a0-a7 for accuracy when advisable */
        del = one_12 / fk;
        del = del * (1. - 4.8 * del * del);
        v = difmuk / fk;
        if (fabs(v) <= 0.25) px = fk * v * v * (((((((a[7] * v + a[6]) * v + a[5]) * v + a[4]) * v + a[3]) * 
                                                  v + a[2]) * v + a[1]) * v + a[0]) - del;
        else px = fk * log(1. + v) - difmuk - del; /* |v| > 1/4 */
        py = R_1_SQRT_2PI / sqrt (fk);
      }
      x = (0.5 - difmuk) / s;
      x *= x; /* x^2 */
      fx = -0.5 * x;
      fy = omega * (((c3 * x + c2) * x + c1) * x + c0);
      if (kflag) {
        /* Step H. Hat acceptance (E is repeated on rejection) */
        if (c * fabs(u) <= py * exp(px + E) - fy * exp(fx + E)) break;
      } 
      else if (fy - u * fy <= py * exp (px - fx)) break;/* Step Q. Quotient acceptance (rare case) */
    }/* t > -.67.. */
  }
  return pois;
}




/* calculates log|gamma(x)| with the sign of the gamma function to the address in the 
 * second argument if this is not NULL. */
double 
biomcmc_lgammafn (double x, int *sgn)
{
  /* For IEEE double precision DBL_EPSILON = 2^-52 = 2.220446049250313e-16 */
  static const double xmax  = 2.5327372760800758e+305; /* =DBL_MAX/log(DBL_MAX)=2^1024/(1024*log(2))=2^1014/log(2)*/
  static const double dxrel = 1.490116119384765696e-8; /* = sqrt(DBL_EPSILON) = 2^-26 = 5^26 * 1e-26 */
  double ans, y, sinpiy;

  if (sgn != NULL) *sgn = 1;

  /* fmod (x,y) returns the floating-point remainder of dividing x by y. */
  if ((sgn != NULL) && (x < 0) && (fmod (floor (-x), 2.) == 0)) *sgn = -1;

  if (x <= 0 && x == trunc(x)) return pInf; /* Negative integer arguments are forbidden */

  y = fabs(x);
  if (y <= 10) return log (fabs (biomcmc_gammafn (x)));
  if (y > xmax) return pInf; /* */ 
  if (x > 0) {
    if(x > 1e17) return(x * (log (x) - 1.));
    else if(x > 4934720.) return(R_LN_SQRT_2PI + (x - 0.5) * log (x) - x);
    else return R_LN_SQRT_2PI + (x - 0.5) * log (x) - x + lgammacor (x); /* i.e. y = x > 10 */ 
  }
  /* else: x < -10; y = -x */ 
  sinpiy = fabs (sin (R_PI * y));
  ans = R_LN_SQRT_PId2 + (x - 0.5) * log (y) - x - log (sinpiy) - lgammacor (y);

#ifdef BIOMCMC_DEBUG
  if(fabs((x - trunc(x - 0.5)) * ans / x) < dxrel) 
    biomcmc_error ("The answer is less than half precision because the argument is too near a negative integer");
#endif

  return ans;
}

/*  This function computes the value of the gamma function. */
double 
biomcmc_gammafn (double x)
{
  static const double gamcs[42] = {
    +.8571195590989331421920062399942e-2, +.4415381324841006757191315771652e-2,
    +.5685043681599363378632664588789e-1, -.4219835396418560501012500186624e-2,
    +.1326808181212460220584006796352e-2, -.1893024529798880432523947023886e-3,
    +.3606925327441245256578082217225e-4, -.6056761904460864218485548290365e-5,
    +.1055829546302283344731823509093e-5, -.1811967365542384048291855891166e-6,
    +.3117724964715322277790254593169e-7, -.5354219639019687140874081024347e-8,
    +.9193275519859588946887786825940e-9, -.1577941280288339761767423273953e-9,
    +.2707980622934954543266540433089e-10, -.4646818653825730144081661058933e-11,
    +.7973350192007419656460767175359e-12, -.1368078209830916025799499172309e-12,
    +.2347319486563800657233471771688e-13, -.4027432614949066932766570534699e-14,
    +.6910051747372100912138336975257e-15, -.1185584500221992907052387126192e-15,
    +.2034148542496373955201026051932e-16, -.3490054341717405849274012949108e-17,
    +.5987993856485305567135051066026e-18, -.1027378057872228074490069778431e-18,
    +.1762702816060529824942759660748e-19, -.3024320653735306260958772112042e-20,
    +.5188914660218397839717833550506e-21, -.8902770842456576692449251601066e-22,
    +.1527474068493342602274596891306e-22, -.2620731256187362900257328332799e-23,
    +.4496464047830538670331046570666e-24, -.7714712731336877911703901525333e-25,
    +.1323635453126044036486572714666e-25, -.2270999412942928816702313813333e-26,
    +.3896418998003991449320816639999e-27, -.6685198115125953327792127999999e-28,
    +.1146998663140024384347613866666e-28, -.1967938586345134677295103999999e-29,
    +.3376448816585338090334890666666e-30, -.5793070335782135784625493333333e-31 };

  int i, n;
  double y, sinpiy, value;
  static int ngam = 0;
  static double xmin = 0, xmax = 0., xsml = 0., dxrel = 0.;

  /* Initialize machine dependent constants, the first time gamma() is called. */
  if (!ngam) {
    double m1, m2;
    /*  Typical values (if below fails): ngam 22, xmax  = 171.61447887182298, xmin = -170.5674972726612, 
     *  xsml = 2.2474362225598545e-308, dxrel = 1.490116119384765696e-8; */
    ngam = chebyshev_init (gamcs, 42, DBL_EPSILON/20.);
    gammalims (&xmin, &xmax);
    m1 = log (DBL_MIN); m2 = -log (DBL_MAX);
    xsml = exp (MAX (m1, m2) + 0.01);/* = exp(.01)*DBL_MIN = 2.247e-308 for IEEE */
    dxrel = sqrt (DBL_EPSILON); /* = sqrt(DBL_EPSILON) = 2^{-26}  */ 
  }

  if(x == NaN) return x;
  /* If the argument is exactly zero or a negative integer then return NaN. */
  if ((x == 0) || ((x < 0) && (x == (long)x))) return NaN;
  y = fabs (x);
  if (y <= 10) {
    /* Compute gamma(x) for -10 <= x <= 10 Reduce the interval and find gamma(1 + y) for 0 <= y < 1 first of all. */
    n = x;
    if(x < 0) --n;
    y = x - n;/* n = floor(x)  ==>	y in [ 0, 1 ) */
    --n;
    value = chebyshev_eval (y * 2 - 1, gamcs, ngam) + .9375;
    if (n == 0) return value;/* x = 1.dddd = 1+y */

    if (n < 0) { /* compute gamma(x) for -10 <= x < 1 exact 0 or "-n" checked already above. */
#ifdef BIOMCMC_DEBUG
      if ((x < -0.5) && (fabs (x - (int)(x - 0.5) / x) < dxrel)) 
        biomcmc_error ("less than half precision answer in biomcmc_gammafn() since x=%g too near a neg integer", x);
#endif
      if (y < xsml) { /* The argument is so close to 0 that the result would overflow. */
        if(x > 0) return pInf;
        else return mInf;
      }

      n = -n;
      for (i = 0; i < n; i++) value /= (x + i);
      return value;
    }
    else { /* gamma(x) for 2 <= x <= 10 */
      for (i = 1; i <= n; i++) value *= (y + i);
      return value;
    }
  }
  else {
    /* gamma(x) for	 y = |x| > 10. */
    if (x > xmax) return pInf; /* Overflow */
    if (x < xmin) return 0.;   /* Underflow */

    if((y <= 50) && (y == (int)y)) { /* compute (n - 1)! */
      value = 1.;
      for (i = 2; i < y; i++) value *= i;
    }
    else { /* normal case  (0.918... = ln(sqrt(2pi)) )*/
      value = exp ((y - 0.5) * log (y) - y + R_LN_SQRT_2PI + ((2*y == (int)2*y) ? stirlerr (y) : lgammacor (y)));
    }
    if (x > 0) return value;

#ifdef BIOMCMC_DEBUG
    if (fabs((x - (int)(x - 0.5))/x) < dxrel)
      biomcmc_error ("less than half precision answer in biomcmc_gammafn() since argument too near a neg integer");
#endif

    sinpiy = sin(R_PI * y); /* 3.141592 = PI */
    if (sinpiy == 0) return pInf; /* Negative integer arg - overflow */

    return - R_PI / (y * sinpiy * value);
  }
}


double 
biomcmc_log1p(double x)
{/* Compute the relative error logarithm [ log(1 + x) ] */
  /* series for biomcmc_log1p on the interval [-.375,.375] with weighted error =  6.35e-32, log weighted error = 31.20, 
   * significant figures required = 30.93 and decimal places required = 32.01 */
  static const double alnrcs[43] = {
    +.10378693562743769800686267719098e+1, -.13364301504908918098766041553133e+0,
    +.19408249135520563357926199374750e-1, -.30107551127535777690376537776592e-2,
    +.48694614797154850090456366509137e-3, -.81054881893175356066809943008622e-4,
    +.13778847799559524782938251496059e-4, -.23802210894358970251369992914935e-5,
    +.41640416213865183476391859901989e-6, -.73595828378075994984266837031998e-7,
    +.13117611876241674949152294345011e-7, -.23546709317742425136696092330175e-8,
    +.42522773276034997775638052962567e-9, -.77190894134840796826108107493300e-10,
    +.14075746481359069909215356472191e-10, -.25769072058024680627537078627584e-11,
    +.47342406666294421849154395005938e-12, -.87249012674742641745301263292675e-13,
    +.16124614902740551465739833119115e-13, -.29875652015665773006710792416815e-14,
    +.55480701209082887983041321697279e-15, -.10324619158271569595141333961932e-15,
    +.19250239203049851177878503244868e-16, -.35955073465265150011189707844266e-17,
    +.67264542537876857892194574226773e-18, -.12602624168735219252082425637546e-18,
    +.23644884408606210044916158955519e-19, -.44419377050807936898878389179733e-20,
    +.83546594464034259016241293994666e-21, -.15731559416479562574899253521066e-21,
    +.29653128740247422686154369706666e-22, -.55949583481815947292156013226666e-23,
    +.10566354268835681048187284138666e-23, -.19972483680670204548314999466666e-24,
    +.37782977818839361421049855999999e-25, -.71531586889081740345038165333333e-26,
    +.13552488463674213646502024533333e-26, -.25694673048487567430079829333333e-27,
    +.48747756066216949076459519999999e-28, -.92542112530849715321132373333333e-29,
    +.17578597841760239233269760000000e-29, -.33410026677731010351377066666666e-30,
    +.63533936180236187354180266666666e-31, };

  if (x == 0.)  return 0.;
  if (x == -1.) return mInf;
  if (x < -1)   return NaN;

  if (fabs(x) <= .375) {
    /* Improve on speed (only); again give result accurate to IEEE double precision: */
    if(fabs(x) < .5 * DBL_EPSILON) return x;
    if( ((0 < x) && (x < 1e-8)) || (-1e-9 < x && x < 0)) return x * (1 - .5 * x);
    /* "22" for IEEE double precision where DBL_EPSILON =  2.22044604925031e-16 */
    return x * (1 - x * chebyshev_eval(x / .375, alnrcs, 22));
  }
#ifdef BIOMCMC_DEBUG
  if (x < -0.999999985) fprintf (stderr, "biomcmc DEBUG: Answer less than half precision since x=%lf too near -1 in biomcmc_log1p()\n", x);
#endif
  return log(1 + x);
}

double biomcmc_log1pmx (double x)
{/* Accurate calculation of log(1+x)-x, particularly for small x.  */
  static const double minLog1Value = -0.79149064;
  static const double two[4] = { 2./9., 2./7., 2./5., 2./3.};

  if (x > 1. || x < minLog1Value) return biomcmc_log1p(x) - x;
  else { /* expand in	 [x/(2+x)]^2 */
    double term = x / (2. + x);
    double y = term * term;
    if (fabs(x) < 1e-2) return term * ((((two[0] * y + two[1]) * y + two[2]) * y + two[3]) * y - x);
    else return term * (2 * y * logcf (y, 3, 2, 1e-14) - x);
  }
}

double biomcmc_expm1 (double x)
{/* Compute exp(x) - 1 accurately also when x is close to zero, i.e. |x| << 1 */
  double y, a = fabs(x);

  if (a < DBL_EPSILON) return x;
  if (a > 0.697) return exp(x) - 1.;  /* negligible cancellation */

  if (a > 1e-8) y = exp(x) - 1;
  else y = (x / 2. + 1.) * x; /* Taylor expansion, more accurate in this range */
  /* Newton step for solving log(1 + y) = x for y : WARNING: does not work for y ~ -1: bug in 1.5.0 */
  return (y - ((1. + y) * (biomcmc_log1p (y) - x)));
}

/* code from Gnu Scientific Library gsl-1.14/randist/discrete.c [http://www.gnu.org/software/gsl/]
 * original Copyright (C) 1996 - 2009 James Theiler, Brian Gough (GNU public license)

 Random Discrete Events: Given K discrete events with different probabilities P[k] produce a value k consistent with its probability.

 * Based on: Alastair J Walker, An efficient method for generating discrete random variables with general distributions, ACM Trans
 * Math Soft 3, 253-256 (1977).  See also: D. E. Knuth, The Art of Computer Programming, Volume 2 (Seminumerical algorithms), 3rd
 * edition, Addison-Wesley (1997), p120.

 * Walker's algorithm does some preprocessing, and provides two arrays: floating point F[k] and integer A[k].  A value k is chosen
 * from 0..K-1 with equal likelihood, and then a uniform random number u is compared to F[k].  If it is less than F[k], then k is
 * returned.  Otherwise, A[k] is returned. 

 * Walker spoke of using two random numbers (an integer 0..K-1, and a floating point u in [0,1]), but Knuth points out that one can just
 * use the integer and fractional parts of K*u where u is in [0,1]. In fact, Knuth further notes that taking F'[k]=(k+F[k])/K, one can
 * directly compare u to F'[k] without having to explicitly set u=K*u-int(K*u).
 */

dsample_stack
new_dsample_stack (size_t vector_size) 
{
  dsample_stack s;
  s = (dsample_stack) biomcmc_malloc (sizeof (struct dsample_stack_struct));
  s->N = vector_size;
  s->i = 0; /* indicates stack is empty */
  s->v = (size_t *) malloc (vector_size * sizeof (size_t));
  return s;
}

void 
del_dsample_stack (dsample_stack s)
{
  if (s) {
    if (s->v) free (s->v);
    free(s);
  }
}

discrete_sample
new_discrete_sample_from_frequencies (double *prob, size_t size)
{ /* preprocessing for Walker's Algorithm */
  size_t  nBigs, nSmalls, i, b, s;
  discrete_sample g;
  dsample_stack Bigs, Smalls;
  double *E, pTotal = 0.0, mean, d;

  if (size < 1) biomcmc_error ("must have a positive number of events for sampling the discrete distribution");

  /* Make sure elements of ProbArray[] are positive, and normalize s.t. sums up to one */
  for (i = 0; i < size; i++) {
    if (prob[i] < 0.) biomcmc_error ("found a negative probability on sampling of empirical frequencies");
    pTotal += prob[i];
  }

  g = (discrete_sample) biomcmc_malloc (sizeof (struct discrete_sample_struct));
  g->K = size;

  if (size == 1) { /* not actually a vector, only one element */
    g->F = NULL; g->A = NULL;
    return g;
  }

  g->F = (double *) biomcmc_malloc (size * sizeof (double));
  g->A = (size_t *) biomcmc_malloc (size * sizeof (size_t));

  E = (double *) biomcmc_malloc (size * sizeof (double));
  for (i = 0; i < size; i++) E[i] = prob[i] / pTotal;

  mean = 1.0/size;
  nSmalls = nBigs = 0;

  for (i = 0; i < size; i++) { /* count number of bigs and smalls */
    if (E[i] < mean) { g->A[i] = 0; nSmalls++; }
    else             { g->A[i] = 1; nBigs++; }
  }

  Bigs   = new_dsample_stack (nBigs);
  Smalls = new_dsample_stack (nSmalls);

  for (i = 0; i < size; i++) {
    if (g->A[i]) {
      if (Bigs->i >= Bigs->N) biomcmc_error ("stack overflow on discrete sampling (for big values)");
      Bigs->v[ Bigs->i++ ] = i; /* stack push */
    }
    else {
      if (Smalls->i >= Smalls->N) biomcmc_error ("stack overflow on discrete sampling (for small values)");
      Smalls->v[ Smalls->i++ ] = i; /* stack push */
    }
  }

  /* Now work through the smalls */
  while (Smalls->i > 0) {
    s = Smalls->v[ --Smalls->i ]; /* stack pop */
    if (Bigs->i == 0) {
      g->A[s] = s;
      g->F[s] = 1.0;
      continue;
    }
    b = Bigs->v[ --Bigs->i ]; /* stack pop */
    g->A[s] = b;
    g->F[s] = size * E[s];
    // fprintf (stderr, "DEBUG new_discrete_sample| s=%2zd, A=%2zd, F=%.4f\n", s, g->A[s], g->F[s]); // DEBUG
    d = mean - E[s];
    E[s] += d;              /* now E[s] == mean */
    E[b] -= d;
    if (E[b] < mean) {
      if (Smalls->i >= Smalls->N) biomcmc_error ("stack overflow on discrete sampling (for small values, when reordering)");
      Smalls->v[ Smalls->i++ ] = b; /* stack push - no longer big, move to stack of smalls */
    }
    else if (E[b] > mean) {
      if (Bigs->i >= Bigs->N) biomcmc_error ("stack overflow on discrete sampling (for big values, when reordering)");
      Bigs->v[ Bigs->i++ ] = b; /* stack push - still big, return it to stack of bigs */
    }
    else {
      /* E[b] == mean implies it is finished too */
      g->A[b] = b;
      g->F[b] = 1.0;
    }
  } // while (Smalls)

  while (Bigs->i > 0) {
    b = Bigs->v[ --Bigs->i ]; /* stack pop */
    g->A[b] = b;
    g->F[b] = 1.0;
  }

  /* Stacks have been emptied, and A and F have been filled */
  if (Smalls->i) biomcmc_error ("problem emptying Small stack of discrete frequencies' sampling"); 

  /* For convenience, set F'[k]=(k+F[k])/K (Knuth convention for saving some computation) */
  for (i = 0; i < size; i++) {
    g->F[i] += i;
    g->F[i] /= size;
  }

  del_dsample_stack (Bigs);
  del_dsample_stack (Smalls);
  if (E) free(E);

  return g;
}

void 
del_discrete_sample (discrete_sample g)
{
  if (g) {
    if (g->A) free(g->A);
    if (g->F) free(g->F);
    free(g);
  }
}

size_t
biomcmc_rng_discrete (discrete_sample g)
{
  size_t c = 0;
  double u, f;

  if (g->K == 1) return 0;

  u = biomcmc_rng_unif();
  c = u * (g->K);
  f = g->F[c];

  // fprintf (stderr, "DEBUG rng_discrete| c, f, u: %zd  %.4f  %f\n", c, f, u); // DEBUG

  if (f == 1.0) return c;
  if (u < f) return c;
  else       return g->A[c];
}

double
biomcmc_discrete_sample_pdf (discrete_sample g, size_t k)
{
  size_t i;
  double f, p = 0;
  if (k > g->K) return 0;
  if (g->K == 1) return ((k == 0)? 1. : 0.);
  for (i = 0; i < g->K; i++) {
    f = g->F[i];
    f = (g->K * f) - i;
    if (i == k) p += f; 
    else if (k == g->A[i]) p += 1.0 - f;
  }

  return p/(double)(g->K);
}

/* Compute log (exp (logx) + exp (logy)) without overflow and without loss of accuracy. */
double 
biomcmc_logspace_add (double logx, double logy)
{
  if (logx > logy) return logx + biomcmc_log1p (exp (logy - logx));
  else             return logy + biomcmc_log1p (exp (logx - logy));
}

/* Compute log (exp (logx) - exp (logy)) without overflow and without loss of accuracy. */
double 
biomcmc_logspace_sub (double logx, double logy)
{
  return logx + biomcmc_log1p (-exp (logy - logx));
}

bool biomcmc_isfinite (double x)
{
  return ((x != NaN) & (x != pInf) & (x != mInf));
}

