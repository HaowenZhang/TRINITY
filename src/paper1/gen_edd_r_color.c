#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include "../base/observations.h"
#include "../base/smf.h"
#include "../base/all_smf.h"
#include "../base/distance.h"
#include "../base/integrate.h"
#include "../base/mlist.h"
#include "../base/calc_sfh.h"
#include "../base/expcache2.h"
#include "../base/universe_time.h"
#include "../base/param_default.h"

extern double param_default[];

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
  // if (argc<2+NUM_PARAMS) {
  //   fprintf(stderr, "Usage: %s mass_cache (mcmc output)\n", argv[0]);
  //   exit(1);
  // }
  if (argc<2) {
    fprintf(stderr, "Usage: %s mass_cache (mcmc output)\n", argv[0]);
    exit(1);
  }

  // for (i=0; i<NUM_PARAMS; i++)
  //   smf.params[i] = atof(argv[i+2]);

  // Read in the model parameter values if provided by the user
  if (argc >= 2+NUM_PARAMS)
    for (i=0; i<NUM_PARAMS; i++)
      smf.params[i] = atof(argv[i+2]);
  // Otherwise use our default values
  else
    for (i=0; i<NUM_PARAMS; i++)
      smf.params[i] = param_default[i];

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
  printf("#1+z yrs M_h M_b M_bh SFR dM_bh/dt Edd_r. bh_eta Obs.ER L_typ BH_Merge_Rate f_active acc_rate_obs SM eta_avg L_kin obs_uv gal_mr\n");
 
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
    // double mu = (1.0 - f) * steps[i].smhm.mu + f * steps[i].smhm.mu;
    double z = zp1 - 1;
    double mass_real = 13.5351-0.23712*z+2.0187*exp(-z/4.48394);
    for (m=8; m<15; m+=0.05) {
      //double mass_real = 13.5351-0.23712*z+2.0187*exp(-z/4.48394); 
      if (m >= mass_real) continue;
      double mf = (m-M_MIN)*BPDEX+0.5;
      int64_t j = mf;
      mf -= j;
      double log_bh_mass, log_bh_acc_rate, log_edd, log_edd_obs, avg_l, efficiency, bh_merge_rate, bh_unmerged, log_bm, log_sfr, gal_mr, bh_eta, obs_uv, frac_active, log_bh_acc_rate_obs, log_sm, eta_rad_avg, eta_kin_avg;
      double l_kin;
      //LINTERP(log_bh_mass,bh_mass_avg);
      INTERP(obs_uv, obs_uv);
      INTERP(log_bh_mass,log_bh_mass);
      LINTERP(log_bh_acc_rate,bh_acc_rate);
      LINTERP(log_bh_acc_rate_obs, bh_acc_rate_obs);
      LINTERP(log_sfr, sfr);
      //INTERP(log_bh_acc_rate_obs, bh_acc_rate_obs);
      //log_bh_acc_rate = log10(log_bh_acc_rate);
      //INTERP(log_sfr,sfr);
      //log_sfr = log10(log_sfr);
      LINTERP(gal_mr, mr);
      INTERP(log_bm, log_bm);
      INTERP(log_sm, log_sm);
      //log_sm += mu;

      INTERP(bh_eta, bh_eta);
      LINTERP(eta_rad_avg, bh_eta_rad_avg);
      INTERP(eta_kin_avg, bh_eta_kin_avg);
      INTERP(bh_merge_rate,bh_merge_rate);
      INTERP(bh_unmerged,bh_unmerged);
      INTERP(frac_active, f_active);
      
      //if (!isfinite(log_bh_acc_rate)) continue;
      log_edd = log_bh_acc_rate-log_bh_mass + log10(4.5e7);
      efficiency = (1.0-f)*steps[i].smhm.bh_efficiency + 
	f*steps[i+1].smhm.bh_efficiency;
      log_edd_obs = log_edd + log10(efficiency) + 1.0;
      avg_l = -5.26 -2.5*(log_edd_obs+log_bh_mass);
      avg_l = -1.0*avg_l;
      bh_merge_rate = log10(bh_merge_rate);
      l_kin = 38.1 + log_bh_mass + eta_kin_avg;
      // bh_unmerged = (bh_unmerged > 0.1 ) ? log10(bh_unmerged) : -1;
      // bh_unmerged = log10(bh_unmerged);
	float years = scale_to_years(1.0 / zp1);
      printf("%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", zp1, years, m, log_bm, log_bh_mass, log_sfr, log_bh_acc_rate, log_edd, bh_eta, log_edd_obs, avg_l, bh_merge_rate, bh_unmerged, frac_active, log_bh_acc_rate_obs, log_sm, eta_rad_avg, l_kin, obs_uv, gal_mr);
    }
  }
  return 0;
}
