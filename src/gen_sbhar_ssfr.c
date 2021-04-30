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
  if (argc<2+NUM_PARAMS) {
    fprintf(stderr, "Usage: %s mass_cache (mcmc output)\n", argv[0]);
    exit(1);
  }
  for (i=0; i<NUM_PARAMS; i++)
    smf.params[i] = atof(argv[i+2]);

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
  printf("#1+z M_h SM M_bh SFR SFR_smloss new_sm dM_bh/dt new_bhm old_bhm old_bhm_sum\n");
 
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

    for (m=8; m<15; m+=0.05) {
      double mf = (m-M_MIN)*BPDEX+0.5;
      int64_t j = mf;
      mf -= j;
      double log_sm, log_bh_mass, 
              log_bh_acc_rate,
              log_sfr, log_sfr_smloss,
              log_new_sm, log_new_bhm,
              log_old_bhm, log_old_bhm_sum;

      // LINTERP(log_sm, sm_avg);
      // LINTERP(log_bh_mass,bh_mass_avg);
      // LINTERP(log_bh_acc_rate,bh_acc_rate);
      // LINTERP(log_sfr, sfr);
      // LINTERP(log_sfr_smloss, sfr_smloss);
      // LINTERP(log_new_sm, new_sm);
      // LINTERP(log_new_bhm, new_bhm);
      // LINTERP(log_old_bhm, old_bh_mass);
      // LINTERP(log_old_bhm_sum, old_bh_mass_sum);

      INTERP(log_sm, sm_avg);
      INTERP(log_bh_mass,bh_mass_avg);
      INTERP(log_bh_acc_rate,bh_acc_rate);
      INTERP(log_sfr, sfr);
      INTERP(log_sfr_smloss, sfr_smloss);
      INTERP(log_new_sm, new_sm);
      INTERP(log_new_bhm, new_bhm);
      INTERP(log_old_bhm, old_bh_mass);
      INTERP(log_old_bhm_sum, old_bh_mass_sum);

      log_sm = log10(log_sm);
      log_bh_mass = log10(log_bh_mass);
      log_bh_acc_rate = log10(log_bh_acc_rate);
      log_sfr = log10(log_sfr);
      log_sfr_smloss = log10(log_sfr_smloss);
      log_new_sm = log10(log_new_sm);
      log_new_bhm = log10(log_new_bhm);
      log_old_bhm = log10(log_old_bhm);
      log_old_bhm_sum = log10(log_old_bhm_sum);

      
      // log_sm += mu;

      if (!isfinite(log_bh_acc_rate)) continue;
      printf("%f %f %f %f %f %f %f %f %f %f\n", zp1, m, log_sm, log_bh_mass, 
                                                    log_sfr, log_sfr_smloss, log_new_sm,
                                                    log_new_bhm, log_old_bhm, log_old_bhm_sum);
    }
  }
  return 0;
}