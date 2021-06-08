// Generate ***total*** (active+inactive) black hole mass functions, broken up
// into the contributions from different halo mass bins, at a given redshift.
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
  double mbh, mh;

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
  int64_t step;
  double f;
  // Calculate the # of snapshot that is the closest to the input redshift.
  calc_step_at_z(z, &step, &f);
  

  double s1 = steps[step].smhm.scatter;
  double s2 = steps[step].smhm.bh_scatter;
  double s = sqrt(s1*s1+s2*s2);
  printf("#BH_Scatter at fixed halo mass: %f\n", s);
  printf("Mbh ND ND(11<Mh<12) ND(12<Mh<13) ND(13<Mh<14) ND(14<Mh<15)\n");

  for (mbh=4; mbh<10.5; mbh+=0.1) 
  {
    // Calculate the total BHMF
    printf("%f %e", mbh, calc_bhmf(mbh, z));
    for (mh=11; mh<15; mh+=1)
    {
      // Calculate the contributions from each halo mass bin.
      printf(" %e", calc_bhmf_mh(mbh, z, mh, mh+1));
    }
    printf("\n");
  }
  return 0;
}
