#include <math.h>
#include "merger_rates.h"

//Log units of mass, ratio
//Returns merger rate / halo / unit redshift / unit log ratio
//Mass in units of Msun (no h).
double merger_rate_macc(double mass, double ratio, double scale) {
  double amp = 0.01329;
  amp *= 1.0+0.249724*log(scale)/(1.0+scale);
  double mscale = pow(10, ((mass+log10(0.7))-12.0)*(0.0515+0.072*scale));
  double beta = -1;
  double gamma = 0.196535;
  double xi = pow(10,gamma*(ratio+3.09109));
  if (ratio > 0) return 0;
  return(amp*mscale*pow(10, beta*ratio)*exp(xi));
}


// double merger_rate_mpeak(double mass, double ratio, double scale) {
//   double amp = 0.0316164;
//   amp *= 1.0+0.20*log(scale)/(1.0+scale);
//   double mscale = pow(10, ((mass-12.0)*(0.03+0.05*scale)));
//   double beta = -1;
//   double gamma = 0.258579;
//   double xi = pow(10,gamma*(ratio+1.92947+0.110*(mass-12.0)));
//   if (ratio > 0) return 0;
//   return(amp*mscale*pow(10, beta*ratio)*exp(xi));
// }

#define LOG10E 0.43429448190325176
double merger_rate_mpeak(double mass, double ratio, double scale) {
  double A0 = 0.14815518 * (mass - 12) - 0.29135703;
  double b0 = 0; 
  double A_scale = -1.60947649 + 3.81591421 * scale + (-2.15254163) * scale * scale ;
  //double A_scale = -1.74183975 + 1.59525058 * scale;
  double b_scale = -1.11364173 + 1.4979694 * scale + (-0.75690147) * scale * scale;
  //double b_scale = -0.91253239 + 0.57824761 * scale;
  if (ratio > 0) return 0;
  double result = A0 + A_scale + (b0 + b_scale) * ratio - LOG10E * pow(10, ratio + 0.5);
  return (pow(10, result));
}

