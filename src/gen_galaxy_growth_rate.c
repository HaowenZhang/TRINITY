// Calculate the galaxy total galaxy growth rates
// (star formation + mergers) as
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
  printf("#1+z Mh GalaxyGrowthRate [Msun/yr]\n");


  for (t=0; t<num_outputs-1; t+=1.0/3.0) 
  {
    i = t;
    double f = t-i;
    double zp1 = (1.0-f)/steps[i].scale + f/steps[i+1].scale;
    double mu = (1.0 - f) * steps[i].smhm.mu + f * steps[i].smhm.mu;
    double dt = (1.0 - f) * steps[i].dt + f * steps[i].dt;

    for (m=8; m<15; m+=0.05) {
      double mf = (m-M_MIN)*BPDEX+0.5;
      int64_t j = mf;
      mf -= j;
      double log_bh_acc_rate, sm, sm_old;

      LINTERP(log_bh_acc_rate,bh_acc_rate);
      INTERP(sm, sm_avg);
      INTERP(sm_old, old_sm);
      sm += mu;
      sm_old += mu;
      // The growth rate is calculated by comparing new and old stellar masses.
      double galaxy_growth_rate = log10((sm - sm_old) / dt);

      if (!isfinite(log_bh_acc_rate)) continue;

      printf("%f %f %f\n", zp1, m, galaxy_growth_rate);
    }
  }
  return 0;
}
