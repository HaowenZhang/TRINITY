// Calculate the average SMBH merger rate as a function of SMBH mass
// at a given redshift.
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
  // Calculate the # of the snapshot that is the closest to
  // the input redshift.
  int64_t step;
  double f;
  calc_step_at_z(z, &step, &f);
  
  
  printf("#BH_Mass BHMR_avg\n");


  for (m=5; m<10.5; m+=0.1) 
  {
    printf("%f %e\n", m, calc_bhmr_mbh(m, z)); //See observations.c for the documentation for calc_bhmr_mbh().
  }
  // fclose(pfile);
  return 0;
}
