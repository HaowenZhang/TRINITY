#include <math.h>

#define STEPS 1024

#ifndef __APPLE__
#ifndef isfinite
#define isfinite(x) finitef(x)
#endif
#endif

// Cosmological parameters. But these particular values
// (i.e., (omega_0, h0) = (0.27, 0.7)) are not used
// if the cached halo mass functions from N-body simulations
// contain their own parameter values.
double omega_0 = 0.27; // dimensionless matter density at z=0
double t_0 = 0;       // Time now (in Hubble time units).
// cached "time table" that can be used to quickly calculate
// the time (i.e., the age of the Universe) given the scale
// factor.
double times[STEPS+1]={0};
double H_CONV = 9.77813952e9/0.7; //Reciprocal of the Hubble constant. 1/(1 km/s/Mpc) in years

// Initialize the time table, given the local dimensionless
// mass density Om, and Hubble constant h0. Of course, we assume
// a flat cosmology.
void init_time_table(double Om, double h0) 
{
  double a = 1;
  double t = t_0;
  double dadt, dtda;
  int i;
  omega_0 = Om;
  H_CONV = 9.77813952e9/h0; // 1/(1 km/s/Mpc) in years

  times[STEPS] = t_0;
  for (i=1; i<=STEPS; i++) 
  {
    // Calculate the da/dt according to the cosmology.
    dadt = sqrt(omega_0 * (1.0/a - (a*a)) + (a*a));
    dtda = 1.0/dadt;
    a -= 1.0/((double)STEPS);
    t -= dtda*(1.0/((double)STEPS));
    times[STEPS-i] = t;
  }
}


// Linearly interpolate between calculated values.
double scale_to_time(double scale) 
{
  double s = scale;
  int l = (int)(s*STEPS);
  double f = s*STEPS - l;
  if (scale > 1) return ((scale-1.0)*STEPS*(times[STEPS]-times[STEPS-1]));
  if (scale < 0) return times[0];
  return (times[l]+f*(times[l+1]-times[l]));
}

// Calculate the age of the Universe in years given
// the scale factor.
double scale_to_years(double scale) 
{
  return (scale_to_time(scale)*H_CONV);
}

// Calculate the scale factor given the age of the 
// Universe in years.
double years_to_scale(double years) 
{
  double htime = years / H_CONV;
  if (htime<=times[0]) return 0;
  double a = 0.5;
  double da = 0.1;
  double stime;
  while (da > 1e-7) 
  {
    stime = scale_to_time(a);
    da = (htime - stime) * da/(scale_to_time(a+da) - stime);
    if (!isfinite(da)) return a;
    if (a+da < 0) a /= 3.0;
    else a+=da;
    da = fabs(da);
    if (da>0.1) da = 0.1;
  }
  return a;
}

