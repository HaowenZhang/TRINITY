#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <inttypes.h>
#include "distance.h"

int main(int argc, char **argv) {
  float v, d, z1, z2, z, zr;
  float true_l = 0, obs_l = 0;
  int64_t i,j;
  if (argc < 4) {
    printf("Usage: %s z1 z2 delta_z\n", argv[0]);
    exit(1);
  }
  init_cosmology(0.27, 1-0.27, 0.70);
  z1 = atof(argv[1]);
  z2 = atof(argv[2]);
  d = atof(argv[3]);
  if (d<=0) { printf("1\n"); return 0; }

  for (i=0; i<50; i++) {
    z = z1 + (z2-z1)*(i+0.5)/50.0;
    float s = d*(1.0+z);
    for (j=-10; j<11; j++) {
      zr = z+j*s/4.0;
      if (zr<=0) continue;
      float min_z = zr-s/8.0;
      if (min_z < 0) min_z = 0;
      float dv = Vc(zr+s/8.0)-Vc(min_z);
      float w = exp(-pow(j/4.0, 2)*0.5);
      true_l += dv*w;
      obs_l +=  pow(luminosity_distance(zr)/luminosity_distance(z), -2)*dv*w;
    }
  }
  printf("%f\n", obs_l/true_l);
  return 0;
}
