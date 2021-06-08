// Calculate the AGN duty cycle (or active fractions) as a function
// of halo mass and redshift. Here, duty cycles are defined as the
// fraction of AGNs with bolometric luminosities or radiative Eddington 
// ratios above an input threshold
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

// Helper functions to do interpolations.
#define INTERP(y,x) { double s1, s2; s1 = (1.0-f)*steps[i].x[j]+f*steps[i+1].x[j];     \
    s2 = (1.0-f)*steps[i].x[j+1]+f*steps[i+1].x[j+1];                    \
    y = s1+mf*(s2-s1); }

#define LINTERP(y,x) { double s1, s2; s1 = (1.0-f)*log10(steps[i].x[j])+f*log10(steps[i+1].x[j]); \
    s2 = (1.0-f)*log10(steps[i].x[j+1])+f*log10(steps[i+1].x[j+1]); \
    y = s1+mf*(s2-s1); }

int main(int argc, char **argv)
{
  int64_t i;
  struct smf_fit smf;
  if (argc < 5) 
  {
    fprintf(stderr, "Usage: %s ledd_or_lum? log(lum_limit/ledd_limit) mass_cache param_file (> output_file)\n", argv[0]);
    exit(1);
  }
  
  // Read in model parameters
  FILE *param_input = check_fopen(argv[4], "r");
  char buffer[2048];
  fgets(buffer, 2048, param_input);
  read_params(buffer, smf.params, NUM_PARAMS);

  int ledd_or_lum = atoi(argv[1]); //if zero, then the threshold is for bolometric luminosity, and
                                   //for Eddington ratio otherwise.
  double log_lim = atof(argv[2]); //The threshold value.

  // Turn off the built-in GSL error handler that kills the program
  // when an error occurs. We handle the errors manually.
  gsl_set_error_handler_off();
  // We use non-linear scaling relation between the radiative and total Eddington ratios.
  nonlinear_luminosity = 1;
  // Set up the PSF for stellar mass functions. See observations.c.
  setup_psf(1);
  // Load cached halo mass functions.
  load_mf_cache(argv[3]);
  // Initialize all the timesteps/snapshots.
  init_timesteps();
  INVALID(smf) = 0;
  // Calculate the star-formation histories and black hole histories. See calc_sfh.c.
  calc_sfh(&smf);
  
  // Calculate the active fractions for halos in each snapshot. Here we ignore the
  // very first snapshot.
  for (i=1; i<num_outputs; i++) 
  {
    calc_active_bh_fraction_lim(i, &smf, ledd_or_lum, log_lim);
  }
  
  printf("#Is the model invalid? %e\n", INVALID(smf));
  double t,m;
  printf("#1+z M_h SM M_bh f_active\n");
 
  // Interpolate between snapshots and redshift bins.
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
