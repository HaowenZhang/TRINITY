#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "observations.h"
#include "smf.h"
#include "all_smf.h"
#include "distance.h"
#include "integrate.h"
#include "mlist.h"
#include "calc_sfh.h"
#include "expcache2.h"
#include "universe_time.h"
//#include "fitter.h"

#define INTERP(y,x) { double s1, s2; s1 = (1.0-f)*steps[i].x[j]+f*steps[i+1].x[j];     \
    s2 = (1.0-f)*steps[i].x[j+1]+f*steps[i+1].x[j+1];			               \
    y = s1+mf*(s2-s1); }

#define LINTERP(y,x) { double s1, s2; s1 = (1.0-f)*log10(steps[i].x[j])+f*log10(steps[i+1].x[j]); \
    s2 = (1.0-f)*log10(steps[i].x[j+1])+f*log10(steps[i+1].x[j+1]);	\
    y = s1+mf*(s2-s1); }

#define ERG_PER_MSUN_KM2_S2 1.988e43
#define F_BARYON 0.157

float fitter(float *params) {
  struct smf_fit test;
  int i;
  for (i=0; i<NUM_PARAMS; i++)
    test.params[i] = params[i];

  assert_model(&test);

  for (i=0; i<NUM_PARAMS; i++)
    params[i] = test.params[i];
  //iterations++;
  float err = all_smf_chi2_err(test);
  if (!isfinite(err) || err<0) return 1e30;
  return err;
}

float calc_chi2(float *params) {
  return fitter(params);
}

float concentration(float hm, float z) {
  return (pow(10, (2.2358 - 0.10*hm))/(1.0+z));
}

float Delta_vir(float z) {
  float om = 0.27;
  float omega = om*(pow(1.0+z,3)/(om*pow(1.0+z, 3) + 1.0-om));
  float x = omega - 1.0;
  return ((1.0/(1.0+x)) * (18*M_PI*M_PI + 82*x - 39*x*x));
}

double x_to_d(double x) {
  double xp1 = 1.0+x;
  return ((log1p(x) - x/xp1)/(x*x*x));
}

double d_to_x(double d) {
  double x = d;
  double td = x_to_d(x);
  double new_x;
  while (fabs((d-td) / d) > 1e-7) {
    double td2 = x_to_d(x+0.1);
    double slope = (td2 - td) / 0.1;
    new_x = x + (d-td)/slope;
    if (new_x < 0) x/=2;
    else x = new_x;
    td = x_to_d(x);
  }
  return x;
}

float ratio_r200c_rvir(float cvir, float vir_dens, float dens_200c) {
  float d = x_to_d(cvir);
  return (d_to_x(d*dens_200c/vir_dens)/cvir);
}

float escape_velocity(float m, float z, float re_thresh) {
  float dvir = Delta_vir(z);
  float scale = 1.0/(1.0+z);
  float uvol = scale*scale*scale;
  float b_dens =  0.27*2.77519737e11/uvol;
  float g = 4.30117902e-9;
  float c = concentration(m, z);
  float rvir = pow((pow(10,m)/(4.0*M_PI*dvir*b_dens/3.0)),1.0/3.0);
  float c_dens = b_dens / (0.27/uvol / (0.27/uvol + 0.73));
  float dens_200c = 200*c_dens;
  float vir_dens = b_dens*dvir;
  float fr_start = (0.015)*ratio_r200c_rvir(c, vir_dens, dens_200c);
  float fr_end = re_thresh * fr_start;
  //float rs = rvir/c;
  //float rho = pow(10, m) / (4.0*M_PI*(pow(rs, 3))*(log(1.0+c)-1.0/(1.0+1.0/c)));
  float d = 1.0/(log(1.0+c)-c/(1.0+c));
  float phi = -d*g*pow(10, m)*(log(1.0+c*fr_start)/(rvir*fr_start)-log(1.0+c*fr_end)/(rvir*fr_end)); // - g*pow(10,m)/(rvir);
  float vesc = sqrt(2.0*fabs(phi));
  return vesc;
}

int main(int argc, char **argv)
{
  int64_t i;
  struct smf_fit smf;
  if (argc<2+NUM_PARAMS) {
    fprintf(stderr, "Usage: %s mass_cache (mcmc output)\n", argv[0]);
    exit(1);
  }
  for (i=0; i<NUM_PARAMS; i++)
    smf.params[i] = atof(argv[i+2]);
  assert_model(&smf);
  gsl_set_error_handler_off();
  nonlinear_luminosity = 1;
  setup_psf(1);
  load_mf_cache(argv[1]);
  init_timesteps();
  INVALID(smf) = 0;
  //double chi2 = calc_chi2(smf.params);
  //printf("Actual chi2=%e\n", chi2);
  calc_sfh(&smf);
  printf("#Is the model invalid? %e\n", INVALID(smf));
  double t,m;
  // int t;
  printf("#1+z M_h M_bh E_SMBH E_vesc vesc\n");
 
 // for (t=0; t<num_outputs-1; t++)
 //  {
 //    for (i = 0; i < M_BINS; ++i)
 //    {
 //      printf("%f %f %f\n", 1 / steps[t].scale, M_MIN + (i + 0.5) / BPDEX, log10(steps[t].bh_mass_avg[i]));
 //    }
 //  }

  for (t=0; t<num_outputs-1; t+=1.0/3.0) {
    i = t;
    double f = t-i;
    double zp1 = (1.0-f)/steps[i].scale + f/steps[i+1].scale;
    double eff = exp10((1.0-f)*steps[i].smhm.bh_efficiency_rad + f*steps[i+1].smhm.bh_efficiency_rad);
    double z = zp1 - 1;
    double mass_real = 13.5351-0.23712*z+2.0187*exp(-z/4.48394);
    for (m=8; m<15; m+=0.05) {
      //double mass_real = 13.5351-0.23712*z+2.0187*exp(-z/4.48394); 
      if (m >= mass_real) continue;
      double mf = (m-M_MIN)*BPDEX+0.5;
      int64_t j = mf;
      mf -= j;
      double log_bh_mass, log_bh_acc_rate;
      double l_kin;
      LINTERP(log_bh_mass,bh_mass_avg);

      LINTERP(log_bh_acc_rate,bh_acc_rate);

      
      //if (!isfinite(log_bh_acc_rate)) continue;
      double E_smbh = eff * exp10(log_bh_mass) * 9e10 * ERG_PER_MSUN_KM2_S2; //9e10 is the c^2 in km^2/s^2.
      // double E_smbh = eff * exp10(log_bh_acc_rate) * 1e9 * 9e10 * ERG_PER_MSUN_KM2_S2; //9e10 is the c^2 in km^2/s^2.
      double vesc = escape_velocity(m, z, 7);
      double E_vesc = 0.5 * F_BARYON * exp10(m) * vesc*vesc * ERG_PER_MSUN_KM2_S2;
      printf("%f %f %f %e %e %e %e\n", zp1, m, log_bh_mass, E_smbh, E_vesc, vesc, eff);
    }
  }
  return 0;
}
