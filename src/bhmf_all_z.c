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
  double m;
  
  if (argc < 3) 
  {
    fprintf(stderr, "Usage: %s mass_cache param_file (> output_file)\n", argv[0]);
    exit(1);
  }
  // Read in model parameters, redshift, BH mass, and lower and upper limits of Eddington ratio.
  FILE *param_input = check_fopen(argv[2], "r");
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
  load_mf_cache(argv[1]);
  // Initialize all the timesteps/snapshots.
  init_timesteps();
  INVALID(smf) = 0;
  // Calculate the star-formation histories and black hole histories. See calc_sfh.c.
  calc_sfh(&smf);

  double z;
  // The redshifts at which quasar luminosities are calculated.
  double zs[12] = {0};
  for (i=2; i<12; i++) zs[i] = i - 1;
  zs[0] = 0.1; zs[1] = 0.5;

  // Calculate quasar luminosity functions at a series of redshifts.
  fprintf(stdout, "#z BH_Mass ND ND(Active)\n");
  for (i=0; i<12; i++)
  {
    int64_t step;
    double f;
    // Find the snapshot that is the closest to the given redshift.
    calc_step_at_z(zs[i], &step, &f);

    for (m=4; m<10.5; m+=0.1) 
    {
      fprintf(stdout, "%f %f %e %e\n", zs[i], m, calc_bhmf(m, zs[i]), calc_active_bhmf(m,zs[i]));
    }
  }

  //// Redo the calculation, but now without merger contributions.
  remove_bh_merger();
  fprintf(stderr, "#z BH_Mass ND ND(Active)\n");
  for (i=0; i<12; i++)
  {
    int64_t step;
    double f;
    calc_step_at_z(zs[i], &step, &f);

    for (m=4; m<10.5; m+=0.1) 
    {
      fprintf(stderr, "%f %f %e %e\n", zs[i], m, calc_bhmf(m, zs[i]), calc_active_bhmf(m,zs[i]));
    }
  }

  return 0;
}
