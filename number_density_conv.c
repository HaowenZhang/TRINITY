#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "observations.h"
#include "mah.h"


float nd_to_mass(float scale, float nd) {
  float mass = 17;
  float tot = 0;
  for (; tot < nd; mass -= 0.01)
    tot += pow(10, mf_cache(scale, mass))/100.0;
  return mass;
}

float mass_to_nd(float scale, float mass) {
  float tot = 0;
  for (; mass<17; mass+=0.01)
    tot += pow(10, mf_cache(scale, mass))/100.0;
  return tot;
}

int main(void) {
  float z1, z2, nd1, nd2;
  load_mf_cache("mf_bolshoi.dat");
  printf("Enter redshift 1: ");
  scanf("%f", &z1);
  printf("Enter log10(number density) [Mpc^-3]: ");
  scanf("%f", &nd1);
  printf("Enter redshift 2: ");
  scanf("%f", &z2);
  
  float m = nd_to_mass(1.0/(z1+1.0), pow(10, nd1));
  printf("Halo virial mass at z = %f: 10^%.2f Msun\n", z1, m);
  m = m_now(1.0/(z1+1.0), m);
  if (z1 && z2)
    printf("Halo virial mass at z = %f: 10^%.2f Msun\n", 0.0, m);
  m = m_at_a(m, 1.0/(z2+1.0));
  printf("Halo virial mass at z = %f: 10^%.2f Msun\n", z2, m);
  nd2 = log10(mass_to_nd(1.0/(z2+1.0), m));
  printf("log10(number density) at z = %f: %f [Mpc^-3]\n", z2, nd2);
  return 0;
}

