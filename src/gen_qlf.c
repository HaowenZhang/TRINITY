// Calculate the quasar bolometric luminosity functions 
// at redshift z, as a function of bolometric luminosity.
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

int main(int argc, char **argv)
{
  int64_t i;
  struct smf_fit smf;
  double l;
  if (argc<3+NUM_PARAMS) 
  {
    fprintf(stderr, "Usage: %s z mass_cache (mcmc output)\n", argv[0]);
    exit(1);
  }
  // Read in the redshift and model parameters. 
  double z = atof(argv[1]);
  for (i=0; i<NUM_PARAMS; i++)
    smf.params[i] = atof(argv[i+3]);

  // Turn off the built-in GSL error handler that kills the program
  // when an error occurs. We handle the errors manually.
  gsl_set_error_handler_off();
  // We use non-linear scaling relation between the radiative and total Eddington ratios.
  nonlinear_luminosity = 1;
  // Set up the PSF for stellar mass functions. See observations.c.
  setup_psf(1);
  // Load cached halo mass functions.
  load_mf_cache(argv[2]);
  // Initialize all the timesteps/snapshots.
  init_timesteps();
  INVALID(smf) = 0;
  // Calculate star formation histories.
  calc_sfh(&smf);
  // Find the closest snapshot to the input redshift.
  int64_t step;
  double f;
  calc_step_at_z(z, &step, &f);

  printf("#BH_alpha: %f\n", steps[step].smhm.bh_alpha);
  printf("#BH_delpha: %f\n", steps[step].smhm.bh_delta);
  // Calculate and output quasar luminosity functions.
  for (l=-10; l>-34.1; l-=0.1) 
  {
    printf("%f %f\n", l, log10(calc_quasar_lf_new(l, z))); //, log10(calc_quasar_lf(l+6.25, z))-2.5);
  }


  return 0;
}
