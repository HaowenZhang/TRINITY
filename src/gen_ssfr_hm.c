// Calculate the galaxy specific star formation rate
// as a function of halo mass and redshift.
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

extern struct timestep *steps;
extern int64_t num_outputs;

int main(int argc, char **argv)
{
  float m;
  struct smf_fit the_smf;
  int i, j;
  if (argc<3+NUM_PARAMS) 
  {
    fprintf(stderr, "Usage: %s mass_cache (mcmc output)\n", argv[0]);
    exit(1);
  }
  // Read in the model parameters
  for (i=0; i<NUM_PARAMS; i++)
    the_smf.params[i] = atof(argv[i+2]);
  the_smf.params[NUM_PARAMS] = 0;

  // Turn off the built-in GSL error handler that kills the program
  // when an error occurs. We handle the errors manually.
  gsl_set_error_handler_off();
  // Set up the PSF for stellar mass functions. See observations.c.
  setup_psf(1);
  // Load cached halo mass functions.
  load_mf_cache(argv[1]);
  // Initialize all the timesteps/snapshots.
  init_timesteps();
  // Calculate the star-formation histories and black hole histories. See calc_sfh.c.
  calc_sfh(&the_smf);

  printf("#z log10(M_halo)[Msun] SSFR[yr^-1]\n");
  for (j=0; j<num_outputs; j++) 
  {
    if (steps[j].scale < 1.0/9.0) continue;
    float z = 1.0/steps[j].scale - 1.0;
    for (i=0; i<M_BINS; i++) 
    {
      m = M_MIN + (i+0.5)*INV_BPDEX;
      if (m<9 || steps[j].t[i] < 1e-9) continue;
      // Since we already calculated and stored the average stellar mass and star formation
      // rates for every halo mass bin, we calcualte the specific star formation rate 
      // using its definition: SSFR = SFR/Mstar.
      printf("%f %f %g\n", z, m, steps[j].sfr[i]/steps[j].sm_avg[i]);
    }
  }
  return 0;
}
