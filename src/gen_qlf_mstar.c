// Generate quasar luminosity functions broken up into
// the contributions from different galaxy mass bins,
// as functions of redshift.
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
  double l, mstar;
  if (argc < 4) 
  {
    fprintf(stderr, "Usage: %s z mass_cache param_file (> output_file)\n", argv[0]);
    exit(1);
  }
  // Read in the redshift and model parameters. 
  double z = atof(argv[1]);

  // Read in model parameters
  FILE *param_input = check_fopen(argv[3], "r");
  char buffer[2048];
  fgets(buffer, 2048, param_input);
  read_params(buffer, smf.params, NUM_PARAMS);

  nonlinear_luminosity = 1;
  // Turn off the built-in GSL error handler that kills the program
  // when an error occurs. We handle the errors manually.
  gsl_set_error_handler_off();
  // Set up the PSF for stellar mass functions. See observations.c.
  setup_psf(1);
  // Load cached halo mass functions.
  load_mf_cache(argv[2]);
  // Initialize all the timesteps/snapshots.
  init_timesteps();
  INVALID(smf) = 0;
  // Calculate the star-formation histories and black hole histories. See calc_sfh.c.
  calc_sfh(&smf);
  int64_t step;
  double f;
  // Calculate the # of snapshot that is the closest to z.
  calc_step_at_z(z, &step, &f);

  // Calculate quasar luminosity functions contributed by different galaxy mass bins.
  for (l=-10; l>-34.1; l-=0.1) 
  {
    for (mstar=8; mstar<=11; mstar+=1)
      printf("%f %f %f %f\n", l, mstar, mstar+1, log10(calc_quasar_lf_mstar(l, z, mstar, mstar+1))); //, log10(calc_quasar_lf(l+6.25, z))-2.5);
  }

  return 0;
}
