#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mah.h"
#include "universe_time.h"

int main(int argc, char **argv) {
  init_time_table(0.27, 0.7);
  double m = 10.177-log10(0.7);
  for (double z = 10.0; z<13.0; z+=0.1) {
    double mar = ma_rate_avg_mnow(m, 1.0/(1.0+z));
    printf("%f %e\n", z, mar/pow(m, 10));
  }
  return 0;
}

