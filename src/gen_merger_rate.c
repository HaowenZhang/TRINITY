// Calculate the numbers of SMBH mergers above mass ratios 1:10 and 1:100
// per halo, as functions of halo mass and redshift. This is for now 
// ***DEPRECATED*** because we do not store the SMBH mass distributions
// that are available for mergers.
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

// This file generates interpolated BH merger rates with mass ratio larger than 1:10 and 1:100.

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
  
  if (argc < 3) 
  {
    fprintf(stderr, "Usage: %s mass_cache parameter_file (> output_file)\n", argv[0]);
    exit(1);
  }

  // Read in model parameters
  FILE *param_input = check_fopen(argv[2], "r");
  char buffer[2048];
  fgets(buffer, 2048, param_input);
  read_params(buffer, smf.params, NUM_PARAMS);

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
  printf("#1+z M_h SM M_bh ND n_merge10 n_merge100\n");
 

  for (t=0; t<num_outputs-1; t+=1.0/3.0) {
    i = t;
    double f = t-i;
    double zp1 = (1.0-f)/steps[i].scale + f/steps[i+1].scale;
    double mu = (1.0 - f) * steps[i].smhm.mu + f * steps[i].smhm.mu;

    for (m=8; m<15; m+=0.05) {
      double mf = (m-M_MIN)*BPDEX+0.5;
      int64_t j = mf;
      mf -= j;
      double log_bh_mass, log_bh_acc_rate, log_sm, log_n_merge10, log_n_merge100, log_nd;

      LINTERP(log_bh_mass,bh_mass_avg);
      LINTERP(log_bh_acc_rate,bh_acc_rate);
      LINTERP(log_sm, sm_avg);
      LINTERP(log_n_merge10, n_merge10);
      LINTERP(log_n_merge100, n_merge100);
      LINTERP(log_nd, t);

      //log_n_merge10 = log10(log_n_merge10);
      //log_n_merge100 = log10(log_n_merge100);
      //log_nd = log10(log_nd);
      log_sm += mu;

      if (!isfinite(log_bh_acc_rate)) continue;
      
	float years = scale_to_years(1.0 / zp1);
      printf("%f %f %f %f %f %f %f\n", zp1, m, log_sm, log_bh_mass, log_nd, log_n_merge10, log_n_merge100);
    }
  }
  return 0;
}
