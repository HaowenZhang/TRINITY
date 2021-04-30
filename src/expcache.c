#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include "expcache.h"

float __expcache_exp10cache[EXP10_CACHE_STEPS+1];
float __expcache_log10cache[EXP10_INV_STEP+1];

#define exp10cache __expcache_exp10cache
#define log10cache __expcache_log10cache
#ifdef __GNUC__
#define likely(x)       __builtin_expect((x),1)
#define unlikely(x)     __builtin_expect((x),0)
#else
#define likely(x)   (x)
#define unlikely(x) (x)
#endif


float exp10fc(float x) {
  float f;
  int b;
  if (unlikely(!(x > EXP10_CACHE_MIN))) return 0;
  if (unlikely(!(x < EXP10_CACHE_MAX))) return HUGE_VAL;
  f = (x*EXP10_INV_STEP) - (EXP10_CACHE_MIN*EXP10_INV_STEP);
  b = f;
  f -= b;
  return (__expcache_exp10cache[b] + (exp10cache[b+1]-exp10cache[b])*f);
}


float log10fc(float x) {
  float f;
  int b,c;
  union {
    float tf;
    uint32_t ti;
  } conv;
  uint32_t d;
  if (unlikely(x<=0)) return -1000;
  conv.tf = x;
  d = conv.ti; //Copy memory from float to int

  //Mask out exponent from float, namely, bits 31-24; we right-shift
  // by 23 places to move bit 24 to the place of bit 1.  The exponent of
  // single-precision floats is biased upwards by 127 from its true value
  // (as opposed to using 2's-complement for negative exponents---biasing
  //  facilitates comparison of two exponents).
  c = ((d>>23)&((1<<8) - 1)) - 127;

  //Mask out the next 10 (or EXP10_LOG2_INV_STEP) bits---these correspond
  // directly to the index in the log cache lookup table.
  b = (d>>(23-EXP10_LOG2_INV_STEP)) & ((1<<EXP10_LOG2_INV_STEP)-1);

  //The remaining bits are the fractional part of the index;
  f = (d&((1<<(23-EXP10_LOG2_INV_STEP))-1)) *
    ((float)(1.0/(double)(1<<(23-EXP10_LOG2_INV_STEP))));

  return ((float)(M_LN2/M_LN10)*c + log10cache[b] + (log10cache[b+1]-log10cache[b])*f);
}

void gen_exp10cache(void) {
  int i;
  double f;
  for (i=0; i<EXP10_CACHE_STEPS+1; i++) {
    f = EXP10_CACHE_STEP*i + EXP10_CACHE_MIN;
    exp10cache[i] = pow(10, f);
  }
  for (i=0; i<EXP10_INV_STEP+1; i++)
    log10cache[i] = log10(1.0 + EXP10_CACHE_STEP*i);

  //Test the cache:
  for (i=-20; i<20; i++) {
    f = exp10fc(i+0.1);
    if (fabs(f/pow(10,i+0.1))-1 > 2e-6) {
      printf("Error in expf calculation: exp10(%f) gave %e\n", i+0.1, f);
      exit(0);
    }
    if (fabs(log10fc(pow(10, i+0.1))/(i+0.1))-1 > 2e-6) {
      printf("Error in logf calculation: log10(%e) gave %f\n", f, log10fc(f));
      exit(0);
    }
  }
}
