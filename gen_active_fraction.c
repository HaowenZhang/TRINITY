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

int main(int argc, char **argv)
{
  int64_t i;
  struct smf_fit smf;
  if (argc<4+NUM_PARAMS) {
    fprintf(stderr, "Usage: %s ledd_or_lum? log(lum_limit/ledd_limit) mass_cache (mcmc output)\n", argv[0]);
    exit(1);
  }
  for (i=0; i<NUM_PARAMS; i++)
    smf.params[i] = atof(argv[i+4]);
  int ledd_or_lum = atoi(argv[1]); //if zero, then luminosity, else ledd.
  double log_lim = atof(argv[2]);
  
  assert_model(&smf);
  gsl_set_error_handler_off();
  nonlinear_luminosity = 1;
  setup_psf(1);
  load_mf_cache(argv[3]);
  init_timesteps();
  INVALID(smf) = 0;
  //double chi2 = calc_chi2(smf.params);
  //printf("Actual chi2=%e\n", chi2);
  calc_sfh(&smf);
  
  for (i=1; i<num_outputs; i++) 
  {

    calc_active_bh_fraction_lim(i, &smf, ledd_or_lum, log_lim);

  }
  
  printf("#Is the model invalid? %e\n", INVALID(smf));
  double t,m;
  // int t;
  printf("#1+z M_h SM M_bh f_active\n");
 
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
    double mu = (1.0 - f) * steps[i].smhm.mu + f * steps[i].smhm.mu;
    double z = zp1 - 1;
    double mass_real = 13.5351-0.23712*z+2.0187*exp(-z/4.48394);
    for (m=8; m<15; m+=0.05) {
      //double mass_real = 13.5351-0.23712*z+2.0187*exp(-z/4.48394); 
      if (m >= mass_real) continue;
      double mf = (m-M_MIN)*BPDEX+0.5;
      int64_t j = mf;
      mf -= j;
      double log_bh_mass, frac_active, log_sm;
      double l_kin;
      INTERP(log_bh_mass,log_bh_mass);

      INTERP(log_sm, log_sm);

      LINTERP(frac_active, f_active);
      if (!isfinite(frac_active)) frac_active = -10;
      printf("%f %f %f %f %f\n", zp1, m, log_sm, log_bh_mass, frac_active);
    }
  }
  return 0;
}
