#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "distance.h"

int main(int argc, char **argv) {
  float app_mag, z, d;
  init_cosmology(0.27, 0.73, 0.7);
  if (argc < 3) {
    printf("Usage: %s app_mag z\n", argv[0]);
    return 1;
  }

  app_mag = atof(argv[1]);
  z = atof(argv[2]);
  d = luminosity_distance(z);
  d/=1e-5; //10 pc
  printf("%f", app_mag - log10(d)*5.0);
  return 0;
}
