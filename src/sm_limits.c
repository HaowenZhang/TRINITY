#include "sm_limits.h"

#define ALL_MIN_MASS 8.3

float sm_min_mass(float z)
{
  float a = 1.0/(z+1.0);
  float sm = 10.819 - 3.1919*a;
  if (z < 0.2) return 7.2;
  if (z < 2) return 8.2;
  if (z < 2.1) return 8.7;
  if (z < 6.1) return 8.9;
  if (z < 7.1) return 7.1;
  if (z < 100) return 7.6;
  return sm;
}

float _sm_max_mass(float z) {
  if (z<0.2) return 12.03;
  if (z>7.5) return 9.5;
  if (z>6.5) return 10.2;
  if (z>5.9) return 10.8;
  if (z>4.9) return 11.2;
  if (z>3.9) return 11.25;
  return 11.8;
}

float sm_max_mass(float z) {
  return (_sm_max_mass(z)-(0.07+0.04*z));
}

float hm_max_mass(float z) {
  if (z<0.2) return 15.2;
  if (z<1.1) return 14.7;
  if (z<2.2) return 14.1;
  if (z<3.2) return 13.6;
  if (z<4.2) return 13.3;
  if (z<5.2) return 12.8;
  if (z<6.2) return 12.5;
  if (z<7.2) return 12.2;
  return 12;
}
