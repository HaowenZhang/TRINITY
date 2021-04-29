#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "universe_time.h"

// The ***remaining*** fraction of stellar mass due
// to stellar evolutions between two scale factors,
// a1 and a2. Reference: Section 2.3 of Behroozi et
// al. (2013)
float _sm_loss(float a1, float a2) 
{
  float t1, t2, t;
  if (a1>a2) return 1.0;
  t1 = scale_to_years(a1);
  t2 = scale_to_years(a2);
  t = t2-t1;
  return (1.0-0.05*log(1 + t/1.4e6));
}

// The ***remaining*** fraction of stellar mass due
// to stellar evolutions loss at a scale factor a3,
// for stars that have formed between scale factors
// a1 and a2.
float calc_sm_loss_int(float a1, float a2, float a3) 
{
#define INT_STEPS 11
  int i;
  float a, mul;
  float step = (a2-a1)/((double)(INT_STEPS-1));
  double sum = _sm_loss(a1,a3)+_sm_loss(a2,a3);

  // Divide the (a1, a2) period into INT_STEPS sub-bins.
  for (i=1; i<INT_STEPS-1; i++) 
  {
    a = a1+i*step;
    mul = (i&1) ? 4 : 2;
    sum += mul*_sm_loss(a, a3);
  }
  sum *= 1.0/((INT_STEPS-1)*3.0);
  return sum;
#undef INT_STEPS
}
