#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "distance.h"

int main(int argc, char **argv) {
  float v, d, z1, z2, a, s;
  if (argc < 4) {
    printf("Usage: %s z1 z2 survey area (arcmin^2) [Om h]\n", argv[0]);
    exit(1);
  }
  double Om = 0.27;
  double h = 0.7;

  if (argc>4) Om = atof(argv[4]);
  if (argc>5) h = atof(argv[5]);
  init_cosmology(Om, 1.0-Om, h);
  
  z1 = atof(argv[1]);
  z2 = atof(argv[2]);
  if (z2 < z1) { d = z2; z2 = z1; z1 = d; }
  a = atof(argv[3])/(4*180*180/M_PI*60.0*60.0);
  s = sqrt(fabs(a)*4*M_PI);

  v = Vc(z2) - Vc(z1);
  printf("%.2e Mpc^3 #Volume\n", v*a);
  d = Dc(z2) - Dc(z1);
  printf("%.2f Mpc #Pencil-Beam Length\n", d);
  printf("%.2f Mpc #Average side length\n", sqrt(v*a/d));
  printf("%.2f Mpc #Beginning side length\n", s*Dc(z1));
  printf("%.2f Mpc #Ending side length\n", s*Dc(z2));
  printf("%.2f' #Opening Angle\n", s*180/M_PI*60);
  return 0;
}
