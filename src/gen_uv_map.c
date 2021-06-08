// Calculate the median stellar mass, star formation rate, 
// median observed UV magnitude and its scatter, as
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

#define REAL_ND_CUTOFF 1e-9

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
  printf("#1+z M_h SM SFR M_UV scatter_UV\n");
 

  for (t=0; t<num_outputs-1; t+=1.0) 
  {
    i = t;
    double f = t-i;
    double zp1 = (1.0-f)/steps[i].scale + f/steps[i+1].scale;
    double mu = (1.0 - f) * steps[i].smhm.mu + f * steps[i].smhm.mu;
    double csfr_completeness = (1.0 - f) * steps[i].smhm.csfr_completeness + f * steps[i].smhm.csfr_completeness;
    double smloss = steps[i].smloss[i];
    double icl_frac = steps[i].smhm.icl_frac;
    for (m=7.1; m<=16.1; m+=0.20) 
    {
      double mf = (m-M_MIN)*BPDEX+0.5;
      int64_t j = mf;
      mf -= j;
      double log_sm, log_sfr, nd, sfrac, sfr, new_sm_sf, new_sm_merger, sm_icl, obs_uv, std_uv;

      sfr = steps[i].sfr[j];
      obs_uv = steps[i].obs_uv[j];
      std_uv = steps[i].std_uv[j];
      log_sm = log10(steps[i].sm_avg[j]);

      printf("%f %f %f %e %f %f\n", zp1, m, log_sm, sfr, obs_uv, std_uv);
    }
  }
  return 0;
}
