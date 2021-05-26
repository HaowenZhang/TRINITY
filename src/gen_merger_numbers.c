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
  if (argc<2+NUM_PARAMS) 
  {
    fprintf(stderr, "Usage: %s mass_cache (mcmc output)\n", argv[0]);
    exit(1);
  }
  for (i=0; i<NUM_PARAMS; i++)
    smf.params[i] = atof(argv[i+2]);
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
  printf("#1+z M_h SM M_bh ND n_merge10 n_merge100\n");
 

  for (t=0; t<num_outputs-1; t+=1.0/3.0) 
  {
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

      log_sm += mu;

      if (!isfinite(log_bh_acc_rate)) continue;
      
	    float years = scale_to_years(1.0 / zp1);
      printf("%f %f %f %f %f %f %f\n", zp1, m, log_sm, log_bh_mass, log_nd, log_n_merge10, log_n_merge100);
    }
  }
  return 0;
}
