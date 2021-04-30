#ifndef INTEGRATE_2_H
#define INTEGRATE_2_H

void init_integration(int i);
double adaptiveGauss(double (*f)(double, void*),   // ptr to function
		     void *extra_data,
		     double a, double b,  // interval [a,b]
		     double rel_err,  // error tolerance
		     int integral_level); 

#endif /* INTEGRATE_2_H */
