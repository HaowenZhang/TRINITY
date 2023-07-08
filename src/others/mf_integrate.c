#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "distance.h"
#include "check_syscalls.h"
#include "observations.h"

#define BPZ 100

int main(int argc, char **argv) {
  float d, z1, z2;
  if (argc < 4) {
    printf("Usage: %s z1 z2 mf_bolshoi.dat\n", argv[0]);
    exit(1);
  }
  load_mf_cache(argv[3]);
  
  z1 = atof(argv[1]);
  z2 = atof(argv[2]);
  if (z2 < z1) { d = z2; z2 = z1; z1 = d; }
  double dz, tw=0, m, *zs = NULL, *weights = NULL;
  int64_t i, num_zs;
  
  num_zs = fabs(z2-z1)*BPZ;
  if (num_zs < 2) num_zs = 2;

  check_realloc_s(zs, sizeof(double), num_zs);
  check_realloc_s(weights, sizeof(double), num_zs);
  dz = fabs(z2-z1)/((double)num_zs+1);
  float zstart = z1;
  if (z2 < zstart) zstart = z2;
  zstart += dz/2.0;
  for (i=0; i<num_zs; i++) {
    double z = zstart + dz*i;
    zs[i] = z;
    weights[i] = Vc(z+dz/2.0)-Vc(z-dz/2.0);
    tw += weights[i];
  }
  for (i=0; i<num_zs; i++) weights[i]/=tw;

  printf("#ND from z=%f to z=%f\n", z1, z2);
  printf("#M ND@z=%f Av.ND(%f,%f)\n", 0.5*(z1+z2), z1, z2);
  for (m=7.0; m<16.0; m+=0.01) {
    double nd = 0;
    for (i=0; i<num_zs; i++) {
      nd += pow(10, mf_cache(1.0/(1.0+zs[i]), m))*weights[i];
    }
    printf("%f %e %e\n", m, pow(10, mf_cache(1.0/(1.0+0.5*(z1+z2)), m)), nd);
  }
  return 0;
}
