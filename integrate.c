#include <math.h>

/* From wikipedia */
//
// Recursive auxiliary function for adaptiveSimpsons() function below
//
double adaptiveSimpsonsAux(double (*f)(double, void *), void *extra_data,
			  double a, double b, double epsilon,
			  double S, double fa, double fb, double fc, int bottom) {
  double c = (a + b)/2, h = b - a;
  double d = (a + c)/2, e = (c + b)/2;
  double fd = f(d, extra_data), fe = f(e, extra_data);
  double Sleft = (h/12)*(fa + 4*fd + fc);
  double Sright = (h/12)*(fc + 4*fe + fb);
  double S2 = Sleft + Sright;
  if (bottom <= 0 || fabs(S2 - S) <= 15*epsilon)
    return S2 + (S2 - S)/15;
  return adaptiveSimpsonsAux(f, extra_data, a, c, epsilon/2, Sleft,  fa, fc, fd, bottom-1) +
    adaptiveSimpsonsAux(f, extra_data, c, b, epsilon/2, Sright, fc, fb, fe, bottom-1);
}

//
// Adaptive Simpson's Rule
//
double adaptiveSimpsons(double (*f)(double, void*),   // ptr to function
		       void *extra_data,
                           double a, double b,  // interval [a,b]
                           double epsilon,  // error tolerance
                           int maxRecursionDepth) {   // recursion cap
  double c = (a + b)/2, h = b - a;
  double fa = f(a, extra_data), fb = f(b, extra_data), fc = f(c, extra_data);
  double S = (h/6)*(fa + 4*fc + fb);
  return adaptiveSimpsonsAux(f, extra_data, a, b, epsilon, S, fa, fb, fc, maxRecursionDepth);
}


