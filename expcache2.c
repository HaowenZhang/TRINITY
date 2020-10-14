#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include "expcache2.h"

double __expcache_exp10cache[EXP10_CACHE_STEPS+1];
double __expcache_log10cache[LOG10_INV_STEP+1];
double __expcache_logexpcache[LOGEXP_CACHE_STEPS+1];
double __expcache_invexpcache[INVEXP_CACHE_STEPS+1];
double __expcache_invexpcache_d[INVEXP_CACHE_STEPS+1];
double __expcache_invexpcache_d2[INVEXP_CACHE_STEPS+1];

#define exp10cache __expcache_exp10cache
#define log10cache __expcache_log10cache
#define logexpcache __expcache_logexpcache
#define invexpcache __expcache_invexpcache
#define invexpcache_d __expcache_invexpcache_d
#define invexpcache_d2 __expcache_invexpcache_d2
#ifdef __GNUC__
#define likely(x)       __builtin_expect((x),1)
#define unlikely(x)     __builtin_expect((x),0)
#else
#define likely(x)   (x)
#define unlikely(x) (x)
#endif

//Returns 1.0/(1+exp(10^-x))
double inv_1p_exp_10_mx(double x) {
  double f;
  int64_t b;
  if (likely(!(x < INVEXP_CACHE_MAX))) return 0.5*(1-exp10(-x));
  if (likely(!(x > INVEXP_CACHE_MIN))) return 0;
  f = (x*INVEXP_INV_STEP) - (INVEXP_CACHE_MIN*INVEXP_INV_STEP);
  b = f;
  f -= b;
  return (invexpcache[b]+f*invexpcache_d[b]+f*f*invexpcache_d2[b]);
}

//Returns Log10(exp10(x)+1)
double log10_1p_exp10fc(double x) {
  double f;
  int64_t b;
  if (unlikely(!(x > LOGEXP_CACHE_MIN))) return (exp10(x)/M_LN10);
  if (unlikely(!(x < LOGEXP_CACHE_MAX))) return x;
  f = (x*LOGEXP_INV_STEP) - (LOGEXP_CACHE_MIN*LOGEXP_INV_STEP);
  b = f;
  f -= b;
  return (logexpcache[b]+f*(logexpcache[b+1]-logexpcache[b]));
}


double exp10fc(double x) {
  double f;
  int64_t b;
  if (unlikely(!(x > EXP10_CACHE_MIN))) return 0;
  if (unlikely(!(x < EXP10_CACHE_MAX))) return HUGE_VAL;
  f = (x*EXP10_INV_STEP) - (EXP10_CACHE_MIN*EXP10_INV_STEP);
  b = f;
  f -= b;
  double x2 = f*(M_LN10/EXP10_INV_STEP);
  double x22 = x2*x2;
  return (exp10cache[b]*(1.0+x2+0.5*x22+(1.0/6.0)*x22*x2));
}


double log10fc(double x) {
  double f;
  int64_t b,c;
  union {
    double tf;
    uint64_t ti;
  } conv;
  uint64_t d;
  if (unlikely(x<=0)) return -1000;
  conv.tf = x;
  d = conv.ti; //Copy memory from float to int

  //Mask out exponent from float, namely, bits 31-24; we right-shift
  // by 23 places to move bit 24 to the place of bit 1.  The exponent of
  // single-precision floats is biased upwards by 127 from its true value
  // (as opposed to using 2's-complement for negative exponents---biasing
  //  facilitates comparison of two exponents).
  c = ((d>>52)&((1<<11) - 1)) - 1023;

  //Mask out the next 10 (or EXP10_LOG2_INV_STEP) bits---these correspond
  // directly to the index in the log cache lookup table.
  b = (d>>(52-LOG10_LOG2_INV_STEP)) & ((1<<LOG10_LOG2_INV_STEP)-1);

  //The remaining bits are the fractional part of the index;
  f = (d&(((uint64_t)1<<(52-LOG10_LOG2_INV_STEP))-1)) *
    ((double)(1.0/(double)((uint64_t)1<<(52-LOG10_LOG2_INV_STEP))));

  double res = ((double)(M_LN2/M_LN10)*(double)c + log10cache[b] + (log10cache[b+1]-log10cache[b])*f);
  return res;
}

void gen_exp10cache(void) {
  int i;
  double f;
  for (i=0; i<EXP10_CACHE_STEPS+1; i++) {
    f = EXP10_CACHE_STEP*i + EXP10_CACHE_MIN;
    exp10cache[i] = pow(10, f);
  }
  for (i=0; i<LOG10_INV_STEP+1; i++)
    log10cache[i] = log10(1.0 + LOG10_CACHE_STEP*i);

  for (i=0; i<LOGEXP_CACHE_STEPS+1; i++) {
    f = LOGEXP_CACHE_STEP*i + LOGEXP_CACHE_MIN;
    logexpcache[i] = log1p(exp10(f))/M_LN10;
  }

  for (i=0; i<INVEXP_CACHE_STEPS+1; i++) {
    f = INVEXP_CACHE_STEP*i + INVEXP_CACHE_MIN;
    invexpcache[i] = 1.0/(1.0+exp(pow(10,-1.0*f)));
    invexpcache_d[i] = INVEXP_CACHE_STEP*exp10(-f)*log(10)*exp(exp10(-f))/pow(1.0+exp(exp10(-f)),2);
    invexpcache_d2[i] = -0.5*INVEXP_CACHE_STEP*INVEXP_CACHE_STEP*
      exp(exp10(-f))*(exp(exp10(-f))*(exp10(-f)-exp10(-2.0*f))+exp10(-f)+exp10(-2.0*f))*pow(log(10),2)/pow(1.0+exp(exp10(-f)),3);
  }

  //Test the cache:
  for (i=-200; i<200; i++) {
    double g = i+0.1;
    f = exp10fc(g);
    if (fabs(f/pow(10,i+0.1)-1.0) > 2e-12) {
      printf("Error in expf calculation: exp10(%f) gave %.12e instead of %.12e\n", g, f, exp10(g));
      exit(0);
    }
    if (fabs(log10fc(pow(10, i+0.1))/(i+0.1)-1.0) > 2e-10) {
      printf("Error in logf calculation: log10(%e) gave %.12f\n", f, log10fc(f));
      exit(0);
    }
    f = log10_1p_exp10fc(g);
    if (fabs(f/(log1p(exp10(g))/M_LN10)-1.0) > 2e-7) {
      printf("Error in log10_1p_exp10fc calculation: gave %.12e instead of %.12e at %f\n", f, log1p(exp10(g))/M_LN10, g);
      exit(0);
    }
    /*f = inv_1p_exp_10_mx(g);
    if (fabs(f-1.0/(1.0+exp(pow(10,-1.0*g)))) > 2e-12) {
      printf("Error in inv_1p_exp_10_mx calculation: gave %.12e instead of %.12e at %f\n", f, 1.0/(1.0+exp(pow(10,-1.0*g))), g);
      exit(0);
      }*/
  }
}
