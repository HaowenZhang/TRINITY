#ifndef EXPCACHE_H
#define EXPCACHE_H

#include <math.h>
#include <stdint.h>

#define EXP10_CACHE_MIN -30
#define EXP10_CACHE_MAX 30
#define EXP10_LOG2_INV_STEP 10
#define EXP10_INV_STEP  (1<<EXP10_LOG2_INV_STEP)
#define EXP10_CACHE_STEPS ((EXP10_CACHE_MAX-EXP10_CACHE_MIN)*EXP10_INV_STEP)
#define EXP10_CACHE_STEP (1.0/(double)(EXP10_INV_STEP))

/* Necessary constants from math.h */
#ifndef M_LN10
#define M_LN10      2.30258509299404568401799145468436421   /* log_e 10 */
#endif
#ifndef M_LN2  
#define M_LN2       0.693147180559945309417232121458176568  /* log_e 2 */
#endif

#ifndef exp10f
#define exp10f(x) expf((float)M_LN10*(x))
#endif
#ifndef exp10
#define exp10(x) exp(M_LN10*(x))
#endif

extern float __expcache_exp10cache[EXP10_CACHE_STEPS+1];
extern float __expcache_log10cache[EXP10_INV_STEP+1];

#define exp10cache __expcache_exp10cache
#define log10cache __expcache_log10cache
#ifdef __GNUC__
#define likely(x)       __builtin_expect((x),1)
#define unlikely(x)     __builtin_expect((x),0)
#else
#define likely(x)   (x)
#define unlikely(x) (x)
#endif


extern inline float exp10fc(float x) {
  float f;
  int b;
  if (unlikely(!(x > EXP10_CACHE_MIN))) return 0;
  if (unlikely(!(x < EXP10_CACHE_MAX))) return HUGE_VAL;
  f = (x*EXP10_INV_STEP) - (EXP10_CACHE_MIN*EXP10_INV_STEP);
  b = f;
  f -= b;
  return (__expcache_exp10cache[b] + (exp10cache[b+1]-exp10cache[b])*f);
}

extern inline float log10fc(float x) {
  float f;
  int b,c;
  uint32_t d;
  union {
    float tf;
    uint32_t ti;
  } conv;
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

#undef exp10cache
#undef log10cache
#undef likely
#undef unlikely

void gen_exp10cache(void);

#endif
