#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mah.h"
#include "universe_time.h"

#define NORMAL_BINS (1024*8+1)
double p_to_normal[NORMAL_BINS];

double perc_to_normal(double p) {
  double x1 = -14;
  double x2 = 14;
  while (x2-x1 > 1e-7) {
    double half = 0.5*(x1+x2);
    double perc = 0.5*(1.0+erf(half));
    if (perc > p) { x2 = half; }
    else { x1 = half; }    
  }
  return ((x1+x2)*M_SQRT1_2);
}

void populate_pcache(void) {
  int64_t i=0;
  for (i=0; i<NORMAL_BINS; i++) {
    double p = (double)i / ((double)NORMAL_BINS-1);
    p_to_normal[i] = perc_to_normal(p);
  }
}

double pcache(double p) {
  double f = p*(NORMAL_BINS-1.0);
  if (f<0) return p_to_normal[0];
  if (f>=(NORMAL_BINS-1)) return p_to_normal[NORMAL_BINS-1];
  int64_t b = f;
  f -= b;
  return (p_to_normal[b] + f*(p_to_normal[b+1]-p_to_normal[b]));
}

int main(void) {
  int64_t i = 0;
  double z = 0;
  double z_max = 8.0;
  init_time_table(0.27,0.7);
  populate_pcache();
  printf("%f %f %f %f\n", pcache(0.84134474606), pcache(0.15865525393), pcache(0.02275013194), pcache(0.97724986805));
  return 0;
  //  printf("%f\n", m_evolution_avg(9, 1.0/16.0, 1.0/13.5)); 
  return 0;
}
