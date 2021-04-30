#include <math.h>
#include <assert.h>

#ifndef __APPLE__
#ifndef isfinite
#define isfinite(x) finitef(x)
#endif
#endif

#define MAX_Z 100.0
#define Z_BINS 1000.0
#define TOTAL_BINS (((int)MAX_Z)*((int)Z_BINS))
double Omega_M, Omega_L, h;
double Dh;
double _Dc[TOTAL_BINS];


double _E (double z) {
  double z1 = 1.0+z;
  return sqrt(Omega_M * (z1*z1*z1) + Omega_L);
}

void init_cosmology(double omega_m, double omega_l, double h0)
{
  int i;
  double z;
  double Dc_int = 0;
  Omega_M = omega_m;
  Omega_L = omega_l;
  h = h0;
  Dh = 2997.92458 / h; //In Mpc  (Speed of light / H0)
  for (i=0; i<TOTAL_BINS; i++) {
    z = (float)(i+0.5) / Z_BINS;
    _Dc[i] = Dc_int * Dh;
    Dc_int += 1.0 / (_E(z)*Z_BINS);
  }
}

double redshift(double a) {
  return (1.0/a-1.0);
}

double scale_factor(double z) {
  return (1.0/(1.0+z));
}

double comoving_distance(double z) {
  float f = z*Z_BINS;
  int bin = f;
  if (z<0) return 0;
  if (bin>(TOTAL_BINS-2)) return (_Dc[TOTAL_BINS-1]);
  f -= bin;
  return (_Dc[bin]*(1.0-f) + _Dc[bin+1]*f);
}

double comoving_distance_h(double z) {
  return (comoving_distance(z)*h);
}

double transverse_distance(double z) {
  return (comoving_distance(z));
}

double angular_diameter_distance(double z) {
  return (transverse_distance(z) / (1.0 + z));
}

double luminosity_distance(double z) {
  return ((1.0+z)*transverse_distance(z));
}

double comoving_volume_element(double z) {
  double z1da = (1.0+z)*angular_diameter_distance(z);
  return (Dh*(z1da*z1da)/_E(z));
}

double comoving_volume(double z) {
  double r = transverse_distance(z);
  return (4.0*M_PI*r*r*r/3.0);
}

double comoving_volume_to_redshift(double Vc) {
  double r = cbrt(Vc*(3.0/(4.0*M_PI)));
  if (r<=0) return 0;
  double z = 1;
  double dz = 0.1;
  double rt;
  while (dz > 1e-7) {
    rt = transverse_distance(z);
    dz = (r - rt) * dz/(transverse_distance(z+dz) - transverse_distance(z));
    if (!isfinite(dz)) return z;
    if (z+dz < 0) z /= 3.0;
    else z+=dz;
    dz = fabs(dz);
    if (dz>0.1) dz = 0.1;
  }
  return z;
}
