// Print out the mass distributions of unmerged (or wandering) SMBHs for 
// each snapshot and halo mass bin. This is deprecated for now because
// we currently do not store such distributions.
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
  int64_t i, j, k;
  struct smf_fit smf;
  if (argc<2+NUM_PARAMS) 
  {
    fprintf(stderr, "Usage: %s mass_cache (mcmc output)\n", argv[0]);
    exit(1);
  }
  for (i=0; i<NUM_PARAMS; i++)
    smf.params[i] = atof(argv[i+2]);
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
  printf("#Is the model invalid? %e\n", INVALID(smf));
  double t,m;s
  // Printing.
  printf("#1+z M_h SM M_bh bh_unmerged_dist\n");
  for (i=0; i<num_outputs-1; i++) 
  {
    double zp1 = (1.0)/steps[i].scale;
    double mu = steps[i].smhm.mu;

    for (j=0; j<M_BINS; j++) 
    {
      double m = M_MIN + (j + 0.5) * INV_BPDEX;
      printf("%f %f %f %f ", zp1, m, log10(steps[i].sm_avg[j]), log10(steps[i].bh_mass_avg[j]));
      for (k=0; k<MBH_BINS; k++) printf(" %e", steps[i].bh_unmerged_dist[j*MBH_BINS+k]);
      printf("\n");
    }
  }
  return 0;
}
