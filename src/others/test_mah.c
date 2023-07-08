#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mah.h"
#include "universe_time.h"

int main(int argc, char **argv) {
  double m0 = atof(argv[1]);
  init_time_table(0.27, 0.7);
  for (double z = 0; z<15; z+=0.1) {
    double m = m_evolution_avg(m0, 1, 1.0/(1.0+z));
    printf("%f %f\n", z, m);
  }
  return 0;
}

