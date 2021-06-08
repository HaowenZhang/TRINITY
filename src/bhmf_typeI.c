// Generate the ***type I quasar*** black hole 
// mass functions at a given redshift. Such predictions
// are compared with SDSS observations from Kelly & Shen (2013).
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
  // Set up the PSF for stellar mass functions. See observations.c.
  setup_psf(1);
  // Load cached halo mass functions.
  load_mf_cache(argv[2]);
  // Initialize all the timesteps/snapshots.
  init_timesteps();
  INVALID(smf) = 0;
  // Calculate the star-formation histories and black hole histories. See calc_sfh.c.
  calc_sfh(&smf);
  // Calculate the # of snapshot that is the closest to the input redshift.
  int64_t step;
  double f;
  calc_step_at_z(z, &step, &f);

  printf("#BH_Alpha: %f\n", steps[step].smhm.bh_alpha);
  printf("#BH_delpha: %f\n", steps[step].smhm.bh_delta);
  // calculate the total scatter in BH mass at fixed ***halo mass***,
  // which is a quadratic sum of the scatter around the median black
  // hole mass--bulge mass relation, and that of the median stellar
  // mass--halo mass relation, enhanced by the slope of the black
  // hole mass--bulge mass relation.
  double s1 = steps[step].smhm.scatter*steps[step].smhm.bh_gamma;
  double s2 = steps[step].smhm.bh_scatter;
  double s = sqrt(s1*s1+s2*s2);
  printf("#BH_Scatter at fixed halo mass: %f\n", s);
  printf("#BH_Mass ND\n");

  for (m=4; m<11; m+=0.03) 
  {
    printf("%f %e\n", m, calc_bhmf_typeI(m, z));
  }
  return 0;
}
