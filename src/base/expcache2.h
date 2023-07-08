#ifndef EXPCACHE_H
#define EXPCACHE_H

#include <math.h>
#include <stdint.h>

#define EXP10_CACHE_MIN -250
#define EXP10_CACHE_MAX 250
#define EXP10_LOG2_INV_STEP 10
#define EXP10_INV_STEP  (1<<EXP10_LOG2_INV_STEP)
#define LOG10_LOG2_INV_STEP 16
#define LOG10_INV_STEP  (1<<LOG10_LOG2_INV_STEP)
#define EXP10_CACHE_STEPS ((EXP10_CACHE_MAX-EXP10_CACHE_MIN)*EXP10_INV_STEP)
#define EXP10_CACHE_STEP (1.0/(double)(EXP10_INV_STEP))
#define LOG10_CACHE_STEP (1.0/(double)(LOG10_INV_STEP))

#define LOGEXP_CACHE_MIN -12
#define LOGEXP_CACHE_MAX  12
#define LOGEXP_LOG2_INV_STEP 12
#define LOGEXP_CACHE_STEPS ((LOGEXP_CACHE_MAX-LOGEXP_CACHE_MIN)*LOGEXP_INV_STEP)
#define LOGEXP_INV_STEP (1<<LOGEXP_LOG2_INV_STEP)
#define LOGEXP_CACHE_STEP (1.0/(double)(LOGEXP_INV_STEP))

#define INVEXP_CACHE_MIN -2
#define INVEXP_CACHE_MAX  10
#define INVEXP_LOG2_INV_STEP 12
#define INVEXP_CACHE_STEPS ((INVEXP_CACHE_MAX-INVEXP_CACHE_MIN)*INVEXP_INV_STEP)
#define INVEXP_INV_STEP (1<<INVEXP_LOG2_INV_STEP)
#define INVEXP_CACHE_STEP (1.0/(double)(INVEXP_INV_STEP))

/* Necessary constants from math.h */
#ifndef M_LN10
#define M_LN10      2.30258509299404568401799145468436421   /* log_e 10 */
#endif
#ifndef M_LN2  
#define M_LN2       0.693147180559945309417232121458176568  /* log_e 2 */
#endif

#ifndef exp10f
#define exp10f(x) expf((double)M_LN10*(x))
#endif
#ifndef exp10
#define exp10(x) exp(M_LN10*(x))
#endif

extern double __expcache_exp10cache[EXP10_CACHE_STEPS+1];
extern double __expcache_log10cache[LOG10_INV_STEP+1];
extern double __expcache_logexpcache[LOGEXP_CACHE_STEPS+1];
extern double __expcache_invexpcache[INVEXP_CACHE_STEPS+1];

#define exp10cache __expcache_exp10cache
#define log10cache __expcache_log10cache
#define logexpcache __expcache_logexpcache
#define invexpcache __expcache_invexpcache
#ifdef __GNUC__
#define likely(x)       __builtin_expect((x),1)
#define unlikely(x)     __builtin_expect((x),0)
#else
#define likely(x)   (x)
#define unlikely(x) (x)
#endif


// extern inline //Returns 1.0/(1+exp(10^-x))

double inv_1p_exp_10_mx(double x); 
/*
{
  double f;
  int64_t b;
  if (likely(!(x < INVEXP_CACHE_MAX))) return 1;
  if (likely(!(x > INVEXP_CACHE_MIN))) return 0;
  f = (x*INVEXP_INV_STEP) - (INVEXP_CACHE_MIN*INVEXP_INV_STEP);
  b = f;
  f -= b;
  return (invexpcache[b]+f*(invexpcache[b+1]-invexpcache[b]));
}
*/


//Returns Log10(exp10(x)+1)
//extern inline 
double log10_1p_exp10fc(double x); 
/*
{
  double f;
  int64_t b;
  if (unlikely(!(x > LOGEXP_CACHE_MIN))) return (exp10(x)/M_LN10);
  if (unlikely(!(x < LOGEXP_CACHE_MAX))) return x;
  f = (x*LOGEXP_INV_STEP) - (LOGEXP_CACHE_MIN*LOGEXP_INV_STEP);
  b = f;
  f -= b;
  return (logexpcache[b]+f*(logexpcache[b+1]-logexpcache[b]));
}
*/

//extern inline 
double exp10fc(double x);
/*
{
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
*/

//extern inline 
double log10fc(double x);
/*
{
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
*/


#undef exp10cache
#undef log10cache
#undef likely
#undef unlikely

void gen_exp10cache(void);

#endif
