#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#define HUGE 3.40282347e+38F
gsl_integration_workspace *w[40] = {0};

void init_integration(int i) {
  if (!w[i]) {
#pragma omp critical
    if (!w[i]) {
      for (i=0; i<40; i++)
	w[i] = gsl_integration_workspace_alloc(4000);
    }
  }
}

double adaptiveGauss(double (*f)(double, void*),   // ptr to function
		     void *extra_data,
		     double a, double b,  // interval [a,b]
		     double rel_err, // error tolerance
		     int integral_level)  
{ 
  gsl_function F;
  double result, error;
  F.function = f;
  F.params = extra_data;
  if (!w[integral_level]) init_integration(integral_level);
 
  int err = gsl_integration_qag(&F, a, b, HUGE, rel_err, 4000, 1, w[integral_level],
		      &result, &error);
  if (err) return 0;
  return result;
}
