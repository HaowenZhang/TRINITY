#ifndef DISTANCE_H
#define DISTANCE_H

#define _DISTANCE_MAX_Z 100.0
#define _DISTANCE_Z_BINS 1000.0
#define _DISTANCE_TOTAL_BINS ((int)(_DISTANCE_MAX_Z*_DISTANCE_Z_BINS))

extern double Dh;

void init_cosmology(double omega_m, double omega_l, double h0);
double redshift(double a);
double scale_factor(double z);

double comoving_distance(double z);
double transverse_distance(double z);
double angular_diameter_distance(double z);
double luminosity_distance(double z);
double comoving_volume_element(double z);
double comoving_volume(double z);
double comoving_volume_to_redshift(double Vc);

#define Dc(z) comoving_distance(z)
#define Dch(z) comoving_distance_h(z)
#define Dm(z) transverse_distance(z)
#define Da(z) angular_diameter_distance(z)
#define Dl(z) luminosity_distance(z)
#define dVc(z) comoving_volume_element(z)
#define Vc(z) comoving_volume(z)

#endif
