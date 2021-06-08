// Output SMBH Eddington ratio distributions at a given redshift, z,
// for each halo mass bin.
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



int main(int argc, char **argv)
{
  int64_t i, j;
  struct smf_fit smf;

  if (argc < 4) 
  {
    fprintf(stderr, "Usage: %s z mass_cache param_file (> output_file)\n", argv[0]);
    exit(1);
  }
  
  // Read in model parameters and redshift.
  double z = atof(argv[1]);
  FILE *param_input = check_fopen(argv[3], "r");
  char buffer[2048];
  fgets(buffer, 2048, param_input);
  read_params(buffer, smf.params, NUM_PARAMS);

  // Fix model parameters
  assert_model(&smf);
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
  // Calculate the star-formation histories and black hole histories. See calc_sfh.c.
  calc_sfh(&smf);
  // Find the snapshot that is the closest to the given redshift.
  int64_t step;
  double f;
  calc_step_at_z(z, &step, &f);

  printf("#z=%f, snapshot=%d\n", z, step);

  // Output SMBH Eddington ratio distributions.
  for (i=0; i<M_BINS; i++)
  {
    for (j=0; j<BHER_BINS; j++)
    {
      printf("%e ", steps[step].bher_dist[i*BHER_BINS+j]);
    }
    printf("\n");
  } 

  return 0;
}
