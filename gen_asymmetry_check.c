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
  printf("#1+z M_h SM obs_uv sfr  M_bh bhar bh_eta BH_Merge_Rate bh_unmerged dt\n");
 
 // for (t=0; t<num_outputs-1; t++)
 //  {
 //    for (i = 0; i < M_BINS; ++i)
 //    {
 //      printf("%f %f %f\n", 1 / steps[t].scale, M_MIN + (i + 0.5) / BPDEX, log10(steps[t].bh_mass_avg[i]));
 //    }
 //  }

  for (i=0; i<num_outputs; i++) {
    double dt = steps[i].dt; 
    double zp1 = 1.0/steps[i].scale;
    double mu = steps[i].smhm.mu;
    double z = zp1 - 1;
    double mass_real = 13.5351-0.23712*z+2.0187*exp(-z/4.48394);
    for (j=0; j<M_BINS; j++) {
      m = M_MIN + (j + 0.5) * INV_BPDEX;  
      
     
    
   
      double log_bh_mass, bh_acc_rate, bh_merge_rate, bh_unmerged, sfr, bh_eta, obs_uv, log_sm;

     log_bh_mass = steps[i].log_bh_mass[j];
     bh_acc_rate = steps[i].bh_acc_rate[j];
     bh_merge_rate = steps[i].bh_merge_rate[j];
     bh_unmerged = steps[i].bh_unmerged[j];
     sfr = steps[i].sfr[j];
     bh_eta = steps[i].bh_eta[j];
     obs_uv = steps[i].obs_uv[j];
     log_sm = steps[i].log_sm[j];
     double nd = steps[i].t[j];
      printf("%f %f %f %f %e %f %e %e %e %e %e %f %f %f %e\n", zp1, m, log_sm, obs_uv, sfr, log_bh_mass, bh_acc_rate, bh_eta, bh_merge_rate, bh_unmerged, dt, steps[i].smhm.scatter, steps[i].smhm.bh_gamma, steps[i].smhm.bh_scatter, nd);
    }
  }
  return 0;
}
