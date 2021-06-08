// Calculate quasar probability distribution functions, i.e.,
// the conditional probabilities of galaxies' hosting AGN at
// an input redshift, as a function of galaxy stellar mass,
// and AGN's bolometric luminosity (i.e., the same form as
// in Bongiorno et al. 2012).
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

  // Calculate and output quasar probability distribution functions.
  printf("#SM LB Prob(LB|SM,z)\n");
  float sm, lb;
  for (sm=8; sm<13; sm+=0.1) 
  {
    for (lb=40; lb<49; lb += 0.25) 
    {
      double prob = calc_qpdf_at_l_m_z(lb, sm, z);
      if (prob <= 0) prob = -1000;
      else prob = log10(prob);
      printf("%f %f %f\n", sm, lb, prob);
    }
  }
  return 0;
}
