// Calculate the chi2's and prior probabilities contributed by each 
// different type of observational data (e.g., stellar mass functions 
// and quasar luminosity functions).
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
  int64_t i;
  struct smf_fit smf;

  if (argc < 2) 
  {
    fprintf(stderr, "Usage: %s mass_cache < ./input_mcmc_chain\n", argv[0]);
    exit(1);
  }

  char buffer[2048];

  // Turn off the built-in GSL error handler that kills the program
  // when an error occurs. We handle the errors manually.
  gsl_set_error_handler_off();
  // We use non-linear scaling relation between the radiative and total Eddington ratios.
  nonlinear_luminosity = 1;
  // Set up the PSF for stellar mass functions. See observations.c.
  setup_psf(1);
  // Initialize the MCMC configurations.
  init_mcmc_from_args(argc, argv);
  // for chi2_smhm, chi2_bhar, chi2_sfr, chi2_icl, chi2_rad, bhz0 and bhz1, see all_smf.c.
  fprintf(stderr, "#chi2_type0...9 chi2_smhm chi2_bhar chi2_sfr chi2_icl chi2_rad bhz0 bhz1\n");
  // Read in the model parameters from the standard input stream.
  while (fgets(buffer, 2048, stdin))
  {
    read_params(buffer, smf.params, NUM_PARAMS);
    INVALID(smf) = 0;
    // Calculate the chi2's and prior probabilities and print them out.
    all_smf_chi2_err_write(smf);
  }
  
  return 0;
}
