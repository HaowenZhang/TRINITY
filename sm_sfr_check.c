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
  int64_t i, j;
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
  printf("#z+1 M_h SM SM_obs SFR SFR_obs\n");

  for (i=0; i<num_outputs; i++) {

    double zp1 = (1.0)/steps[i].scale;
    double mu = steps[i].smhm.mu;
    double kappa = steps[i].smhm.kappa;
    double z = zp1 - 1;
    double sfr_corr = mu + kappa * exp(-(z - 2) * (z - 2) * 0.5);
    double mass_real = 13.5351-0.23712*z+2.0187*exp(-z/4.48394);
    for (j=0; j<M_BINS; j++) {
      //double mass_real = 13.5351-0.23712*z+2.0187*exp(-z/4.48394); 
      m = M_MIN + (j + 0.5) * INV_BPDEX;
      if (m >= mass_real) continue;
   
      double sm = steps[i].log_sm[j];
      double sfr = log10(steps[i].sfr[j]);

      printf("%f %f %f %f %f %f\n", zp1, m, sm, sm + mu, sfr, sfr + sfr_corr); 
    }
  }
  return 0;
}
