#include <stdio.h>
#include "universe_time.h"
#include <math.h>
#include <inttypes.h>
//Mass accretion history fits

#define M13 13.2763173435874 //Mass of calibrated MAH, in Msun (no h).

double m13_at_a(double a) {
  double z = 1.0/a - 1.0;
  return (M13 + 3.0*log10(1.0+z) - 6.11478*log10(1.0+0.5*z) - M_LOG10E*0.503151*z);
}

//In log units
double m_at_a(double m, double a) {
  double a_0 = 0.205268 - log10(exp(M_LN10*0.18*(9.64921-m))+1.0);
  double f = (m - M13)*(1.0+exp(-(1.0-a_0)/0.215))/(1.0+exp(-(a-a_0)/0.215));
  return (m13_at_a(a) + f);
}

double m_scatter_at_a(double m, double a) {
  double scatter_slope = 0.96 + 0.13*(m-13.125+log10(0.7));
  return ((1.0-a)*scatter_slope*0.5);
}

double ma_rate(double m, double a) {
  double m1 = pow(10, m_at_a(m, a-0.005));
  double m2 = pow(10, m_at_a(m, a+0.005));
  double t1 = scale_to_years(a-0.005);
  double t2 = scale_to_years(a+0.005);
  double mar = (m2-m1) / (t2-t1);
  return (mar);
}

double ma_rate_avg(double m, double a) {
  double mar = ma_rate(m,a);
  double da = a - 0.4;
  //return mar*0.85*(0.25*da*da+0.96)*pow(10,a*0.05*(m_at_a(m,a)-11.655))/pow(1.0+pow(10, (11.5-m_at_a(m,a))), 0.06);
  return mar*(0.22*da*da+0.85)*pow(pow(10,m_at_a(m,a)-12.0), a*0.05)/pow(1.0+pow(10, (11-0.5*a-m_at_a(m,a))), 0.04+0.1*a);
}

double m_now(double a, double m0) {
  double m_at_z0 = m0;
  int k;
  for (k=0; k<20; k++)
    m_at_z0 += m0-m_at_a(m_at_z0, a);
  return m_at_z0;
}

double ma_rate_avg_mnow(double m, double a) {
  double mnow = m_now(a, m);
  return ma_rate_avg(mnow, a);
}

double m_evolution_avg(double m, double a1, double a2) {
  int64_t i;
  double m2 = m;
  double steps = 200.0;
  double mar = ma_rate_avg_mnow(m, a1);
  double a_mid = a1 + (a2-a1)/(steps*2.0);
  double m_mid = m + log10(1+mar*(scale_to_years(a_mid)-scale_to_years(a1))/pow(10,m));
  double a = a1;
  double da = (a2-a1)/steps;
  for (i=0; i<steps; i++) {
    double mar = ma_rate_avg_mnow(m_mid, a_mid);
    double dm = mar * (scale_to_years(a+da) - scale_to_years(a));
    if (fabs(dm/pow(10,m2)) > 0.001 && i > 10)
      return m_evolution_avg(m2, a, a2);
    m2 += log10(1.0+dm/pow(10,m2));
    a += da;
    double mar2 = ma_rate_avg_mnow(m2, a);
    m_mid += log10(1.0+ mar2 * (scale_to_years(a_mid+da) - scale_to_years(a_mid)) / pow(10, m_mid));
    a_mid += da;
  }
  return m2;
}
