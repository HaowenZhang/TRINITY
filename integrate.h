#ifndef INTEGRATE_H
#define INTEGRATE_H

/* From wikipedia */
double adaptiveSimpsons(double (*f)(double, void*),   // ptr to function
		       void *extra_data,
                           double a, double b,  // interval [a,b]
                           double epsilon,  // error tolerance
		       int maxRecursionDepth);   // recursion cap        
 
#endif /* INTEGRATE_H */
