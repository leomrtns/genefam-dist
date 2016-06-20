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

#include "prob_distribution_aux.h"

/* log gamma correction factor for x >= 10 so that log(gamma(x)) = .5*log(2*pi) + (x-.5)*log(x) -x + lgammacor(x) */  
double 
lgammacor (double x)
{
  static const double algmcs[15] = {
    +.1666389480451863247205729650822e+0, -.1384948176067563840732986059135e-4,
    +.9810825646924729426157171547487e-8, -.1809129475572494194263306266719e-10,
    +.6221098041892605227126015543416e-13, -.3399615005417721944303330599666e-15,
    +.2683181998482698748957538846666e-17, -.2868042435334643284144622399999e-19,
    +.3962837061046434803679306666666e-21, -.6831888753985766870111999999999e-23,
    +.1429227355942498147573333333333e-24, -.3547598158101070547199999999999e-26,
    +.1025680058010470912000000000000e-27, -.3401102254316748799999999999999e-29,
    +.1276642195630062933333333333333e-30 };
  /* below: 94906265.62425156 = 2^26.5 and  3.745194030963158e306 = DBL_MAX / 48 =  2^1020 / 3 */

  if (x < 10) return NaN;
#ifdef BIOMCMC_DEBUG
  if (x >= 3.745194030963158e306) biomcmc_error ("underflow occurred in lgammacor");
#endif
  if (x < 94906265.62425156) {
    double tmp = (10./x);
    return chebyshev_eval (tmp * tmp * 2 - 1, algmcs, 5) / x;
  }
  return 1. / (x * 12.);
}

/* determines the number of terms for the double precision orthogonal series "dos" needed to insure
   the error is no larger than "eta".  Ordinarily eta will be chosen to be one-tenth machine precision. */
int 
chebyshev_init (const double *dos, int nos, double eta)
{
  int i, ii;
  double err;

  if (nos < 1) return 0;

  err = 0.0;
  i = 0;      /* just to avoid compiler warnings */
  for (ii = 1; ii <= nos; ii++) {
    i = nos - ii;
    err += fabs (dos[i]);
    if (err > eta) return i;
  }
  return i;
}

double 
chebyshev_eval (double x, const double *a, const int n)
{
  double b0, b1, b2, twox;
  int i;

  if ((n < 1) || (n > 1000) || (x < -1.1) || (x > 1.1)) return NaN;

  twox = x * 2;
  b2 = b1 = 0;
  b0 = 0;
  for (i = 1; i <= n; i++) {
    b2 = b1;
    b1 = b0;
    b0 = twox * b1 - b2 + a[n - i];
  }
  return (b0 - b2) * 0.5;
}

/* calculates the minimum and maximum legal bounds for x in biomcmc_gammafn(x). */
void 
gammalims (double *xmin, double *xmax)
{
  double alnbig, alnsml, xln, xold;
  int i;

  alnsml = log (DBL_MIN);
  *xmin = -alnsml;
  for (i=1; i<=10; ++i) {
    xold = *xmin;
    xln = log (*xmin);
    *xmin -= *xmin * ((*xmin + .5) * xln - *xmin - .2258 + alnsml) /
    (*xmin * xln + .5);
    if (fabs (*xmin - xold) < .005) { *xmin = -(*xmin) + .01; goto find_xmax; }
  }

  /* unable to find xmin: in R-mathlib they give up; here we try the IEEE_754 value and proceed */
  *xmin = -170.5674972726612; goto find_xmax;

find_xmax:

  alnbig = log (DBL_MAX);
  *xmax = alnbig;
  for (i=1; i<=10; ++i) {
    xold = *xmax;
    xln = log(*xmax);
    *xmax -= *xmax * ((*xmax - .5) * xln - *xmax + .9189 - alnbig) /
    (*xmax * xln - .5);
    if (fabs (*xmax - xold) < .005) { *xmax += -.01; goto done; }
  }

  /* unable to find xmax: try the IEEE_754 value and proceed */
  *xmax =  171.61447887182298; goto done;

done:

  *xmin = MAX (*xmin, -(*xmax) + 1);
}

/* Continued fraction for calculation of 1/i + x/(i+d) + x^2/(i+2*d) + x^3/(i+3*d) + ... = sum_{k=0}^Inf x^k/(i+k*d)
 * auxilary in biomcmc_log1pmx() and lgamma1p(); eps = relative tolerance */
double
logcf (double x, double i, double d, double eps)
{
  double c1 = 2 * d, c2 = i + d, c4 = c2 + d, a1 = c2;
  double b1 = i * (c2 - i * x), b2 = d * d * x, a2 = c4 * c2 - b2;

  b2 = c4 * b1 - i * b2;

  while (fabs(a2 * b1 - a1 * b2) > fabs(eps * b1 * b2)) {
    double c3 = c2*c2*x;
    c2 += d; c4 += d;
    a1 = c4 * a2 - c3 * a1;
    b1 = c4 * b2 - c3 * b1;

    c3 = c1 * c1 * x;
    c1 += d; c4 += d;
    a2 = c4 * a1 - c3 * a2;
    b2 = c4 * b1 - c3 * b2;

    if (fabs (b2) > scalefactor) { 
      a1 /= scalefactor; b1 /= scalefactor; a2 /= scalefactor; b2 /= scalefactor; 
    }
    else if (fabs (b2) < 1 / scalefactor) {
      a1 *= scalefactor; b1 *= scalefactor; a2 *= scalefactor; b2 *= scalefactor;
    }
  }

  return a2 / b2;
}

double 
lgamma1p (double a)
{/* Compute  log(gamma(a+1))  accurately also for small a (0 < a < 0.5). */
  const double eulers_const =	 0.5772156649015328606065120900824024;
  const int N = 40; /* coeffs[i] holds (zeta(i+2)-1)/(i+2) , i = 0:(N-1), N = 40 */
  static const double coeffs[40] = { /* (zeta(2)-1)/2, (zeta(3)-1)/3 ... (zeta(40+1)-1)/(40+1) */
    0.3224670334241132182362075833230126e-0, 0.6735230105319809513324605383715000e-1,
    0.2058080842778454787900092413529198e-1, 0.7385551028673985266273097291406834e-2,
    0.2890510330741523285752988298486755e-2, 0.1192753911703260977113935692828109e-2,
    0.5096695247430424223356548135815582e-3, 0.2231547584535793797614188036013401e-3,
    0.9945751278180853371459589003190170e-4, 0.4492623673813314170020750240635786e-4,
    0.2050721277567069155316650397830591e-4, 0.9439488275268395903987425104415055e-5,
    0.4374866789907487804181793223952411e-5, 0.2039215753801366236781900709670839e-5,
    0.9551412130407419832857179772951265e-6, 0.4492469198764566043294290331193655e-6,
    0.2120718480555466586923135901077628e-6, 0.1004322482396809960872083050053344e-6,
    0.4769810169363980565760193417246730e-7, 0.2271109460894316491031998116062124e-7,
    0.1083865921489695409107491757968159e-7, 0.5183475041970046655121248647057669e-8,
    0.2483674543802478317185008663991718e-8, 0.1192140140586091207442548202774640e-8,
    0.5731367241678862013330194857961011e-9, 0.2759522885124233145178149692816341e-9,
    0.1330476437424448948149715720858008e-9, 0.6422964563838100022082448087644648e-10,
    0.3104424774732227276239215783404066e-10, 0.1502138408075414217093301048780668e-10,
    0.7275974480239079662504549924814047e-11, 0.3527742476575915083615072228655483e-11,
    0.1711991790559617908601084114443031e-11, 0.8315385841420284819798357793954418e-12,
    0.4042200525289440065536008957032895e-12, 0.1966475631096616490411045679010286e-12,
    0.9573630387838555763782200936508615e-13, 0.4664076026428374224576492565974577e-13,
    0.2273736960065972320633279596737272e-13, 0.1109139947083452201658320007192334e-13 };
  const double c = 0.2273736845824652515226821577978691e-12;/* zeta(N+2)-1 */
  double lgam;
  int i;

  if (fabs (a) >= 0.5) return biomcmc_lgammafn (a + 1, NULL);
  /* Abramowitz & Stegun 6.1.33, also  http://functions.wolfram.com/06.11.06.0008.01 */
  lgam = c * logcf (-a / 2, N + 2, 1, 1e-14);
  for (i = N - 1; i >= 0; i--) lgam = coeffs[i] - a * lgam;
  return (a * lgam - eulers_const) * a - biomcmc_log1pmx (a);
}

/* dpois_wrap (x_P_1,  lambda, log_p) == dpois (x_P_1 - 1, lambda, log_p) */
double
dpois_wrap (double x_plus_1, double lambda, bool log_p)
{
  if (!biomcmc_isfinite (lambda)) return log_p ? mInf : 0.;
  if (x_plus_1 > 1.) return dpois_raw (x_plus_1 - 1, lambda, log_p);
  if (lambda > (fabs (x_plus_1 - 1) * M_cutoff)) 
    return log_p? (-lambda - biomcmc_lgammafn (x_plus_1, NULL)) : exp (-lambda - biomcmc_lgammafn (x_plus_1, NULL));
  else {
    double d = dpois_raw (x_plus_1, lambda, log_p);
    return log_p ? d + log (x_plus_1 / lambda) : d * (x_plus_1 / lambda);
  }
}

double 
dpois_raw (double x, double lambda, bool log_p)
{
  /* x >= 0 ; integer for dpois(), but not e.g. for pgamma()! lambda >= 0 */
  if (lambda == 0) return( (x == 0) ? log_p ? 0. : 1. : log_p ? mInf : 0. );
  if (!biomcmc_isfinite (lambda)) return log_p ? mInf : 0.;
  if (x < 0) return log_p ? mInf : 0.;
  if (x <= lambda * DBL_MIN) return log_p ? -lambda : exp (-lambda);
  if (lambda < x * DBL_MIN) 
    return log_p ? (-lambda + x*log(lambda) - biomcmc_lgammafn (x+1, NULL)) : 
    exp (-lambda + x*log(lambda) - biomcmc_lgammafn (x+1, NULL));

  return log_p ? (-0.5 * log (R_2PI * x) - stirlerr (x) - bd0 (x,lambda)) : 
  exp (-stirlerr (x) - bd0 (x,lambda)) / sqrt (R_2PI * x);
}

/* stirlerr(n) = log(n!) - log( sqrt(2*pi*n)*(n/e)^n )
 *             = log Gamma(n+1) - 1/2 * [log(2*pi) + log(n)] - n*[log(n) - 1]
 *             = log Gamma(n+1) - (n + 1/2) * log(n) + n - log(2*pi)/2
 * see also lgammacor() in which computes almost the same! */
double 
stirlerr (double n)
{
  static const double S[5] = {1./12. ,1./360., 1./1260., 1./1680., 1./1188};
  static const double sferr_halves[31] = {/* error for 0, 0.5, 1.0, 1.5, ..., 14.5, 15.0. */
    0.0, 0.1534264097200273452913848, 0.0810614667953272582196702, 0.0548141210519176538961390,
    0.0413406959554092940938221,  0.03316287351993628748511048, 0.02767792568499833914878929,
    0.02374616365629749597132920, 0.02079067210376509311152277, 0.01848845053267318523077934,
    0.01664469118982119216319487, 0.01513497322191737887351255, 0.01387612882307074799874573,
    0.01281046524292022692424986, 0.01189670994589177009505572, 0.01110455975820691732662991,
    0.010411265261972096497478567, 0.009799416126158803298389475, 0.009255462182712732917728637,
    0.008768700134139385462952823, 0.008330563433362871256469318, 0.007934114564314020547248100,
    0.007573675487951840794972024, 0.007244554301320383179543912, 0.006942840107209529865664152,
    0.006665247032707682442354394, 0.006408994188004207068439631, 0.006171712263039457647532867,
    0.005951370112758847735624416, 0.005746216513010115682023589, 0.005554733551962801371038690};
  double nn;

  /* For n > 15, uses the series 1/12n - 1/360n^3 + ...
   * For n <=15, integers or half-integers, uses stored values.
   * For other n < 15, uses lgamma directly */
  if (n <= 15.0) {
    nn = n + n;
    if (nn == (int)nn) return sferr_halves[(int)nn];
    return (biomcmc_lgammafn (n + 1., NULL) - (n + 0.5)*log (n) + n - R_LN_SQRT_2PI); 
  }
  nn = n*n;
  if (n>500.) return ((S[0]-S[1]/nn)/n);
  if (n> 80.) return ((S[0]-(S[1]-S[2]/nn)/nn)/n);
  if (n> 35.) return ((S[0]-(S[1]-(S[2]-S[3]/nn)/nn)/nn)/n);
  return ((S[0]-(S[1]-(S[2]-(S[3]-S[4]/nn)/nn)/nn)/nn)/n); /* 15 < n <= 35 : */
}

/* Evaluates the "deviance part"
 *  bd0(x,M) :=  M * D0(x/M) = M*[ x/M * log(x/M) + 1 - (x/M) ] =  x * log(x/M) + M - x
 *  where M = E[X] = n*p (or = lambda), for   x, M > 0
 *  In a manner that should be stable (with small relative error) for all x and M=np. In particular for x/np close 
 *  to 1, direct evaluation fails, and evaluation is based on the Taylor series of log((1+v)/(1-v)) with 
 *  v = (x-np)/(x+np). */
double 
bd0 (double x, double np)
{
  double ej, s, s1, v;
  int j;

  if(!biomcmc_isfinite (x) || !biomcmc_isfinite (np) || (np == 0.)) return NaN; 

  if (fabs(x-np) < 0.1*(x+np)) {
    v = (x-np)/(x+np);
    s = (x-np)*v;/* s using v -- change by MM */
    ej = 2*x*v;
    v = v*v;
    for (j=1; ; j++) { /* Taylor series */
      ej *= v;
      s1 = s+ej/((j<<1)+1);
      if (s1==s) return (s1); /* last term was effectively 0 */
      s = s1;
    }
  }
  return (x*log (x/np) + np - x); /* else:  | x - np |  is not too small */
}

double
pgamma_smallx (double x, double alph, bool log_p)
{/* Abramowitz and Stegun 6.5.29 [right] */
  double f2, sum = 0, c = alph, n = 0, term;

  /* Relative to 6.5.29 all terms have been multiplied by alph and the first, thus being 1, is omitted. */
  do {
    n++;
    c *= -x / n;
    term = c / (alph + n);
    sum += term;
  } while (fabs (term) > DBL_EPSILON * fabs (sum));

  double f1 = log_p ? biomcmc_log1p (sum) : 1. + sum;
  if (alph > 1) {
    f2 = dpois_raw (alph, x, log_p);
    f2 = log_p ? f2 + x : f2 * exp (x);
  }
  else if (log_p) f2 = alph * log (x) - lgamma1p (alph);
  else f2 = pow (x, alph) / exp (lgamma1p (alph));
  return log_p ? f1 + f2 : f1 * f2;
}

double
pd_upper_series (double x, double y, bool log_p)
{
  double term = x / y;
  double sum = term;

  /* sum =  \sum_{n=1}^ oo  x^n     / (y*(y+1)*...*(y+n-1))    =  \sum_{n=0}^ oo  x^(n+1) / (y*(y+1)*...*(y+n)) = 
   *	   =  x/y * (1 + \sum_{n=1}^oo	x^n / ((y+1)*...*(y+n))) ~  x/y + o(x/y) {which happens when alph -> Inf} */
  do { y += 1.; term *= x / y; sum += term; } while (term > sum * DBL_EPSILON);
  return log_p ? log (sum) : sum;
}

double
pd_lower_series (double lambda, double y)
{
  double term = 1, sum = 0;

  /* sum =  \sum_{n=0}^ oo  y*(y-1)*...*(y - n) / lambda^(n+1) 
   *     =  y/lambda * (1 + \sum_{n=1}^Inf  (y-1)*...*(y-n) / lambda^n   ~  y/lambda + o(y/lambda) */
  while (y >= 1 && term > sum * DBL_EPSILON) { term *= y / lambda; sum += term; y -= 1.; }

  if (y != floor (y)) {
    /* The series does not converge as the terms start getting bigger (besides flipping sign) for y < -lambda. */
    double f;
    f = pd_lower_cf (y, lambda + 1 - y);
    sum += term * f;
  }
  return sum;
}

double
pd_lower_cf (double i, double d)
{/* Continued fraction for calculation of (i/d) + o(i/d) */
  static const int max_it = 200000;
  double f = 0, of, c1, c2, c3, c4, a1, b1, a2, b2;

  if (i < d * 1e-20) return (i/d); /* includes d = Inf,  or i = 0 < d */

  a1 = 0; b1 = 1;
  a2 = i; b2 = d;
  while (b2 > scalefactor) { a1/=scalefactor; b1/=scalefactor; a2/=scalefactor; b2/=scalefactor;}

  if(a2 == 0) return 0;/* just in case, e.g. d=i=0 */

  c2 = a2; c4 = b2; c1 = 0;
  while (c1 < max_it) {
    c1++;	c2--;	c3 = c1 * c2;	c4 += 2;
    a1 = c4 * a2 + c3 * a1;
    b1 = c4 * b2 + c3 * b1;

    c1++;	c2--;	c3 = c1 * c2;	c4 += 2;
    a2 = c4 * a1 + c3 * a2;
    b2 = c4 * b1 + c3 * b2;

    if (b2 > scalefactor) { a1/=scalefactor; b1/=scalefactor; a2/=scalefactor; b2/=scalefactor;}

    if (b2 != 0) {
      of = f;
      f = a2 / b2;
      /* convergence check: relative; absolute for small f : */
      if (fabs (f - of) <= DBL_EPSILON * MAX(1., fabs(f)))
        return f;
    }
  }
#ifdef BIOMCMC_DEBUG
  biomcmc_error ("NON-convergence in pd_lower_cf() f= %g. Shouldn't happen, but not fatal. \n", f);
#endif
  return f;
}

/* Compute the following ratio with higher accuracy that would be had from doing it directly.
 *   \f$ \frac{ dnorm (x, 0, 1, FALSE) } { pnorm (x, 0, 1, FALSE) } \f$
 * Abramowitz & Stegun 26.2.12 */
double
dpnorm (double x, double lp)
{
  /* So as not to repeat a pnorm call, we expect lp == pnorm (x, 0, 1, TRUE)
   * but use it only in the non-critical case where either x is small or p==exp(lp) is close to 1. */
  bool upper_tail = false;
  if (x < 0) { x = -x; upper_tail = true; }
  if ((x > 10) && upper_tail) {
    double term = 1. / x;
    double sum = term;
    double x2 = x * x;
    double i = 1.;

    do { term *= -i / x2; sum += term; i += 2.; } while (fabs (term) > DBL_EPSILON * sum);
    return 1. / sum;
  }
  else return biomcmc_dnorm (x, 0, 1, false) / exp (lp);
}

/* Asymptotic expansion to calculate the probability that Poisson variate has value <= x.
 * Various assertions about this are made (without proof) at http://members.aol.com/iandjmsmith/PoissonApprox.htm */
double
ppois_asymp (double x, double lambda, bool log_p)
{ /* It is called by pgamma_raw with (lower_tail = false); needs change if it is possible to be called with 
     lower_tail = T */
  static const double coefs_a[8] = {
    -1e99, /* placeholder used for 1-indexing */
    2/3., -4/135., 8/2835., 16/8505., -8992/12629925., -334144/492567075., 698752/1477701225. };
  static const double coefs_b[8] = {
    -1e99, /* placeholder */
    1/12., 1/288., -139/51840., -571/2488320., 163879/209018880., 5246819/75246796800., -534703531/902961561600.};
  double elfb, elfb_term;
  double res12, res1_term, res1_ig, res2_term, res2_ig;
  double dfm, pt_, s2pt, f, np;
  int i;

  dfm = lambda - x;
  /* If lambda is large, the distribution is highly concentrated about lambda.  So representation error in x or 
   * lambda can lead to arbitrarily large values of pt_ and hence divergence of the coefficients of this 
   * approximation. */
  pt_ = - biomcmc_log1pmx (dfm / x);
  s2pt = sqrt (2 * x * pt_);
  if (dfm < 0) s2pt = -s2pt;

  res12 = 0;
  res1_ig = res1_term = sqrt (x);
  res2_ig = res2_term = s2pt;
  for (i = 1; i < 8; i++) {
    res12 += res1_ig * coefs_a[i];
    res12 += res2_ig * coefs_b[i];
    res1_term *= pt_ / i ;
    res2_term *= 2 * pt_ / (2 * i + 1);
    res1_ig = res1_ig / x + res1_term;
    res2_ig = res2_ig / x + res2_term;
  }

  elfb = x;
  elfb_term = 1;
  for (i = 1; i < 8; i++) {
    elfb += elfb_term * coefs_b[i];
    elfb_term /= x;
  }
  elfb = -elfb; /* if (upper_tail) in R-mathlib but here we only call within pgamma_raw()  with upper tail) */
  f = res12 / elfb;
  np = biomcmc_pnorm (s2pt, 0.0, 1.0, log_p);

  if (log_p) return np + biomcmc_log1p (f * dpnorm (s2pt, np));
  else return np + f * biomcmc_dnorm (s2pt, 0., 1., log_p);
} 

double 
pgamma_raw (double x, double alph, bool log_p)
{ /* Here, assume that  (x,alph) are not NA  &  alph > 0 . */
  double res, dpsum;

  if(x <= 0.)   return log_p ? mInf : 0.;
  if(x >= pInf) return log_p ?   0. : 1.;

  if (x < 1) res = pgamma_smallx (x, alph, log_p);
  else if ((x <= (alph - 1)) && (x < 0.8 * (alph + 50))) {
    /* incl. large alph compared to x */
    double sum = pd_upper_series (x, alph, log_p);/* = x/alph + o(x/alph) */
    double d = dpois_wrap (alph, x, log_p);
    res = log_p ? sum + d : sum * d;
  } 
  else if (((alph - 1) < x) && (alph < 0.8 * (x + 50))) {
    /* incl. large x compared to alph */
    double sum;
    double d = dpois_wrap (alph, x, log_p);
    if (alph < 1) {
      if (x * DBL_EPSILON > 1 - alph) sum = log_p ? 0. : 1.;
      else {
        double f = pd_lower_cf (alph, x - (alph - 1)) * x / alph;
        /* = [alph/(x - alph+1) + o(alph/(x-alph+1))] * x/alph = 1 + o(1) */
        sum = log_p ? log (f) : f;
      }
    }
    else {
      sum = pd_lower_series (x, alph - 1);/* = (alph-1)/x + o((alph-1)/x) */
      sum = log_p ? biomcmc_log1p (sum) : 1 + sum;
    }
    dpsum = d + sum;
    res = log_p ? ((dpsum > -R_LN2) ? log (-expm1 (dpsum)) : biomcmc_log1p (-exp (dpsum))) : (1. - d * sum);
  } 
  else res = ppois_asymp (alph - 1, x, log_p); /* x > 1 and x fairly near alph. */

  /* We lose a fair amount of accuracy to underflow in the cases where the final result is very close to DBL_MIN. 
   * In those cases, simply redo via log space. with(.Machine, double.xmin / double.eps) #|-> 1.002084e-292 */
  if (!log_p && res < DBL_MIN / DBL_EPSILON) return exp (pgamma_raw (x, alph, true));
  else return res;
}

/* g = log Gamma (nu/2) and tol = EPS1 */
double 
qchisq_appr (double p, double nu, double g, bool log_p, double tol)
{
  static const double C[4] = {4.67, 6.66, 6.73, 13.32};
  double alpha, a, c, ch, p1, p2, q, t, x, rdtmacro;

  if ((p==NaN) || (nu == NaN)) return NaN;
  if ( (log_p  && (p > 0)) || ((!log_p) && ((p < 0) || (p > 1))) ) return pInf;
  if (nu <= 0) return NaN;

  alpha = 0.5 * nu;/* = [pq]gamma() shape */
  c = alpha-1;

  if (nu < (-1.24)*(p1 = (log_p ? p : log (p)))) ch = exp ((log (alpha) + p1 + g)/alpha + R_LN2);	/*small chi-squared*/
  else if (nu > 0.32) {	/*  using Wilson and Hilferty estimate */
    x = biomcmc_qnorm (p, 0, 1, log_p);
    p1 = 2./(9*nu);
    ch = nu * pow (x * sqrt(p1) + 1 - p1, 3);
    /* approximation for p tending to 1: */
    rdtmacro = log_p ? (p > -R_LN2 ? log (-expm1 (p)) : biomcmc_log1p (-exp (p))) : (biomcmc_log1p(-p));
    if( ch > 2.2*nu + 6 ) ch = -2 * (rdtmacro - c * log (0.5 * ch) + g);
  } 
  else { /* "small nu" : 1.24*(-log(p)) <= nu <= 0.32 */
    ch = 0.4;
    rdtmacro = log_p ? (p > -R_LN2 ? log (-expm1 (p)) : biomcmc_log1p (-exp (p))) : (biomcmc_log1p(-p));
    a = rdtmacro + g + c*R_LN2;
    do {
      q = ch;
      p1 = 1. / (1 + ch * (C[0] + ch));
      p2 = ch * (C[2] + ch * (C[1] + ch));
      t = -0.5 + (C[0] +2 * ch) * p1 - (C[2] + ch * (C[3] + 3 * ch)) / p2;
      ch -= (1- exp (a + 0.5 * ch) * p2 * p1) / t;
    } while(fabs(q - ch) > tol * fabs(ch));
  }
  return ch;
}

void pnorm_both(double x, double *cum, double *ccum, int i_tail, bool log_p)
{
  /* i_tail in {0,1,2} means: "lower", "upper", or "both" :
     if(lower) return  *cum := P[X <= x] ; if(upper) return *ccum := P[X >  x] = 1 - P[X <= x] */
  static const double a[5] = {
    2.2352520354606839287, 161.02823106855587881, 1067.6894854603709582, 18154.981253343561249, 0.065682337918207449113
  };
  static const double b[4] = {
    47.20258190468824187, 976.09855173777669322, 10260.932208618978205, 45507.789335026729956
  };
  static const double c[9] = {
    0.39894151208813466764, 8.8831497943883759412, 93.506656132177855979, 597.27027639480026226, 2494.5375852903726711,
    6848.1904505362823326, 11602.651437647350124, 9842.7148383839780218, 1.0765576773720192317e-8
  };
  static const double d[8] = {
    22.266688044328115691, 235.38790178262499861, 1519.377599407554805, 6485.558298266760755,
    18615.571640885098091, 34900.952721145977266, 38912.003286093271411, 19685.429676859990727
  };
  static const double p[6] = {
    0.21589853405795699, 0.1274011611602473639, 0.022235277870649807, 0.001421619193227893466,
    2.9112874951168792e-5, 0.02307344176494017303
  };
  static const double q[5] = {
    1.28426009614491121, 0.468238212480865118, 0.0659881378689285515, 0.00378239633202758244, 7.29751555083966205e-5
  };
  double xden, xnum, temp, del, eps, xsq, y;
  double min = DBL_MIN;
  int i, lower, upper;

  if(x == NaN) { *cum = *ccum = x; return; }
  
  eps = DBL_EPSILON * 0.5; /* Consider changing these */

  /* i_tail in {0,1,2} =^= {lower, upper, both} */
  lower = i_tail != 1; upper = i_tail != 0;

  y = fabs(x);
  if (y <= 0.67448975) { /* qnorm(3/4) = .6744.... -- earlier had 0.66291 */
    if (y > eps) {
      xsq = x * x;
      xnum = a[4] * xsq;
      xden = xsq;
      for (i = 0; i < 3; ++i) {
        xnum = (xnum + a[i]) * xsq;
        xden = (xden + b[i]) * xsq;
      }
    } 
    else xnum = xden = 0.0;

    temp = x * (xnum + a[3]) / (xden + b[3]);
    if(lower)  *cum = 0.5 + temp;
    if(upper) *ccum = 0.5 - temp;
    if(log_p) {
      if(lower)  *cum = log(*cum);
      if(upper) *ccum = log(*ccum);
    }
  }
  else if (y <= R_SQRT_32) { /* Evaluate pnorm for 0.674.. = qnorm(3/4) < |x| <= sqrt(32) ~= 5.657 */
    xnum = c[8] * y;
    xden = y;
    for (i = 0; i < 7; ++i) {
      xnum = (xnum + c[i]) * y;
      xden = (xden + d[i]) * y;
    }
    temp = (xnum + c[7]) / (xden + d[7]);

    xsq = trunc (y * 16) / 16; 
    del = (y - xsq) * (y + xsq);
    if(log_p) {	
      *cum = (-xsq * xsq * 0.5) + (-del * 0.5) + log (temp);	
      if((lower && (x > 0.)) || (upper && (x <= 0.))) 
        *ccum = biomcmc_log1p (-exp (-xsq * xsq * 0.5) * exp(-del * 0.5) * temp);
    }
    else { 
      *cum = exp (-xsq * xsq * 0.5) * exp (-del * 0.5) * temp; 
      *ccum = 1.0 - *cum; 
    }
    /* swap  ccum <--> cum */
    if (x > 0.) { temp = *cum; if(lower) *cum = *ccum; *ccum = temp; }
  }

  /* else	|x| > sqrt(32) = 5.657 : the next two case differentiations were really for lower=T, log=F 
   * Particularly *not* for log_p ! Note that we do want symmetry(0), lower/upper -> hence use y */
  else if(log_p || (lower && -37.5193 < x  &&  x < 8.2924) || (upper && -8.2924  < x  &&  x < 37.5193) ) {
    /* Evaluate pnorm for x in (-37.5, -5.657) union (5.657, 37.5) */
    xsq = 1.0 / (x * x);
    xnum = p[5] * xsq;
    xden = xsq;
    for (i = 0; i < 4; ++i) {
      xnum = (xnum + p[i]) * xsq;
      xden = (xden + q[i]) * xsq;
    }
    temp = xsq * (xnum + p[4]) / (xden + q[4]);
    temp = (R_1_SQRT_2PI - temp) / y; 

    xsq = trunc (x * 16) / 16; 
    del = (x - xsq) * (x + xsq);
    if(log_p) {	
      *cum = (-xsq * xsq * 0.5) + (-del * 0.5) + log (temp);	
      if((lower && (x > 0.)) || (upper && (x <= 0.))) 
        *ccum = biomcmc_log1p (-exp (-xsq * xsq * 0.5) * exp(-del * 0.5) * temp);
    }
    else { 
      *cum = exp (-xsq * xsq * 0.5) * exp (-del * 0.5) * temp; 
      *ccum = 1.0 - *cum; 
    }
    /* swap  ccum <--> cum */
    if (x > 0.) { temp = *cum; if(lower) *cum = *ccum; *ccum = temp; }
  }
  else { /* no log_p , large x such that probs are 0 or 1 */
    if(x > 0) { *cum = 1.; *ccum = 0.;	}
    else {      *cum = 0.; *ccum = 1.;	}
  }
  if(log_p) {
    if(*cum > -min)	  *cum = -0.;
    if(*ccum > -min) *ccum = -0.;
  }
  else {
    if(*cum < min)	 *cum = 0.;
    if(*ccum < min)	*ccum = 0.;
  }
  return;
}

double
do_poisson_search (double y, double *z, double p, double lambda, double incr)
{
  double new_y;
  if(*z >= p) for(;;) { /* search to the left */
    new_y = y - incr;
    if((y == 0) || (*z = biomcmc_ppois (new_y, lambda, false)) < p) return y;
    y =  MAX (0, new_y);
  }
  else for(;;) { /* search to the right */
    y = y + incr;
    if ((*z = biomcmc_ppois (y, lambda, false)) >= p) return y;
  }
}

