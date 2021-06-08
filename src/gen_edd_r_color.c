// Calculate the median stellar mass, SMBH mass, average SFR,
// SMBH accretion rates, merger rates, Eddington ratios, as
// functions of halo mass and redshift.
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

// Helper functions to interpolate between halo mass bins and snapshots.
#define INTERP(y,x) { double s1, s2; s1 = (1.0-f)*steps[i].x[j]+f*steps[i+1].x[j];     \
    s2 = (1.0-f)*steps[i].x[j+1]+f*steps[i+1].x[j+1];			               \
    y = s1+mf*(s2-s1); }

#define LINTERP(y,x) { double s1, s2; s1 = (1.0-f)*log10(steps[i].x[j])+f*log10(steps[i+1].x[j]); \
    s2 = (1.0-f)*log10(steps[i].x[j+1])+f*log10(steps[i+1].x[j+1]);	\
    y = s1+mf*(s2-s1); }

int main(int argc, char **argv)
{
  int64_t i;
  struct smf_fit smf;
  if (argc<2) 
  {
    fprintf(stderr, "Usage: %s mass_cache parameter_file\n", argv[0]);
    exit(1);
  }

  // Read in model parameters

  FILE *param_input = check_fopen(argv[2], "r");
  char buffer[2048];
  fgets(buffer, 2048, param_input);
  read_params(buffer, smf.params, NUM_PARAMS);

  // for (i=0; i<NUM_PARAMS; i++)
  //   smf.params[i] = atof(argv[i+2]);
  // Fix some model parameters.
  assert_model(&smf);
  // Turn off the built-in GSL error handler that kills the program
  // when an error occurs. We handle the errors manually.
  gsl_set_error_handler_off();
  // We use non-linear scaling relation between the radiative and total Eddington ratios.
  nonlinear_luminosity = 1;
  // Set up the PSF for stellar mass functions. See observations.c.
  setup_psf(1);
  // Load cached halo mass functions.
  load_mf_cache(argv[1]);
  // Initialize all the timesteps/snapshots.
  init_timesteps();
  INVALID(smf) = 0;
  // Calculate the star-formation histories and black hole histories. See calc_sfh.c.
  calc_sfh(&smf);
  printf("#Is the model invalid? %e\n", INVALID(smf));
  double t,m;

  printf("#1+z M_h M_star M_bulge M_bh SFR BHAR bh_eta_avg bh_eta_typical BH_Merge_Rate (All in log units except for 1+z)\n");

  for (t=0; t<num_outputs-1; t+=1.0/3.0) 
  {
    i = t;
    double f = t-i;
    double zp1 = (1.0-f)/steps[i].scale + f/steps[i+1].scale;
    double mu = (1.0 - f) * steps[i].smhm.mu + f * steps[i].smhm.mu;
    double z = zp1 - 1;
    // mass_real is the halo mass threshold above which we deem the halos
    // as too rare to care about.
    double mass_real = 13.5351-0.23712*z+2.0187*exp(-z/4.48394);
    for (m=8; m<15; m+=0.05) 
    {
      
      if (m >= mass_real) continue;
      double mf = (m-M_MIN)*BPDEX+0.5;
      int64_t j = mf;
      mf -= j;
      double log_bh_mass, log_bh_acc_rate, log_edd, log_edd_obs, avg_l, efficiency, bh_merge_rate, bh_unmerged, log_bm, log_sfr, gal_mr, bh_eta, obs_uv, frac_active, log_bh_acc_rate_obs, log_sm, eta_rad_avg, eta_kin_avg;
      double l_kin;
      LINTERP(log_bh_mass,bh_mass_avg);
      INTERP(obs_uv, obs_uv);
      LINTERP(log_bh_acc_rate,bh_acc_rate);
      LINTERP(log_bh_acc_rate_obs, bh_acc_rate_obs);
      LINTERP(log_sfr, sfr);
      LINTERP(gal_mr, mr);
      INTERP(log_bm, log_bm);
      INTERP(log_sm, log_sm);

      // bh_eta is the ***typical*** SMBH Eddington ratio, i.e., the
      // break point of the double power-law distribution.
      INTERP(bh_eta, bh_eta);
      INTERP(bh_merge_rate,bh_merge_rate);
      INTERP(bh_unmerged,bh_unmerged);
      INTERP(frac_active, f_active);
      
      //log_edd is the ***average*** SMBH Eddington ratio.
      log_edd = log_bh_acc_rate-log_bh_mass + log10(4.5e7);

      // age of the Universe at the given redshift.
	    float years = scale_to_years(1.0 / zp1);
      printf("%f %f %f %f %f %f %f %f %f %f\n", zp1, m, log_sm, log_bm, log_bh_mass, log_sfr, log_bh_acc_rate, log_edd, bh_eta, bh_merge_rate);
    }
  }
  return 0;
}
