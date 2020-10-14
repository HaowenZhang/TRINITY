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

#define REAL_ND_CUTOFF 1e-9

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

  // nonlinear_luminosity = 1;
  setup_psf(1);
  load_mf_cache(argv[1]);
  init_timesteps();
  INVALID(smf) = 0;
  gsl_set_error_handler_off();
  //double chi2 = calc_chi2(smf.params);
  //printf("Actual chi2=%e\n", chi2);
  calc_sfh(&smf);
  printf("#Is the model invalid? %e\n", INVALID(smf));
  double t,m;
  // int t;
  printf("#1+z M_h SM SFR ND smf new_sm_sf new_sm_merger sfrac lv\n");
 

  for (t=0; t<num_outputs-1; t+=1.0) {
    i = t;
    double f = t-i;
    double zp1 = (1.0-f)/steps[i].scale + f/steps[i+1].scale;
    double mu = (1.0 - f) * steps[i].smhm.mu + f * steps[i].smhm.mu;
    double csfr_completeness = (1.0 - f) * steps[i].smhm.csfr_completeness + f * steps[i].smhm.csfr_completeness;
    double smloss = steps[i].smloss[i];
    double icl_frac = steps[i].smhm.icl_frac;
    for (m=7.1; m<=16.1; m+=0.05) {
      double mf = (m-M_MIN)*BPDEX+0.5;
      int64_t j = mf;
      mf -= j;
      double log_sm, log_sfr, nd, sfrac, sfr, new_sm_sf, new_sm_merger, sm_icl;

      LINTERP(log_sfr, sfr);
      INTERP(sfr, sfr);
      INTERP(log_sm, log_sm);
      INTERP(nd, t);
      INTERP(sfrac, sfrac);
      //sfr = steps[i].sfr[j];
      //log_sm = log10(steps[i].sm_avg[j]);
      nd = steps[i].t[j];
      sm_icl = steps[i].sm_icl[j];
      new_sm_sf = steps[i].new_sm[j] * smloss;
      new_sm_merger = steps[i].sm_icl[j] * icl_frac;
      sfrac = steps[i].sfrac[j];
      double lv = steps[i].lv[j];
      //log_sm += mu;
      if (nd < REAL_ND_CUTOFF) continue;

      printf("%f %f %f %e %e %e %e %e %e %f\n", zp1, m, log_sm, sfr, sfrac, sm_icl, new_sm_sf, new_sm_merger, sfrac, lv);
    }
  }
  return 0;
}
