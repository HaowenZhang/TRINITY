// Print out the SMBH bolometric luminosity distributions
// at a given redshift z, as functions of ***intrinsic***
// SMBH masses.
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

  if (argc < 5) 
  {
    fprintf(stderr, "Usage: %s z scatter_obs mass_cache param_file (> output_file)\n", argv[0]);
    exit(1);
  }
  double z = atof(argv[1]);
  double scatter_obs = atof(argv[2]);

  // Read in model parameters
  FILE *param_input = check_fopen(argv[4], "r");
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
  load_mf_cache(argv[3]);
  // Initialize all the timesteps/snapshots.
  init_timesteps();
  INVALID(smf) = 0;
  // Calculate the star-formation histories and black hole histories. See calc_sfh.c.
  calc_sfh(&smf);
  int64_t step;
  double f;
  // Calculate the star-formation histories and black hole histories. See calc_sfh.c.
  calc_step_at_z(z, &step, &f);
  
  
  printf("#BH_Alpha: %f\n", steps[step].smhm.bh_alpha);
  printf("#BH_delpha: %f\n", steps[step].smhm.bh_delta);
  double s1 = steps[step].smhm.scatter;
  double s2 = steps[step].smhm.bh_scatter;
  double s = sqrt(s1*s1+s2*s2);

  double mbh_min = steps[step].bh_mass_min;
  double mbh_max = steps[step].bh_mass_max;
  double mbh_bpdex = MBH_BINS / (steps[step].bh_mass_max - steps[step].bh_mass_min);
  double mbh_inv_bpdex = 1.0 / mbh_bpdex;

  // Print out the pre-calculated distributions.
  for (i=0; i<MBH_BINS; i++)
  {
    if (mbh_min + (i + 0.5) * mbh_inv_bpdex < 8) continue;
    fprintf(stdout, "%.6f ", mbh_min + (i + 0.5) * mbh_inv_bpdex);
    for (int j=0; j<LBOL_BINS; j++)
      fprintf(stdout, "%.6e ", steps[step].lum_dist_full[i*LBOL_BINS+j]);
    fprintf(stdout, "\n");
  }
  return 0;
}
