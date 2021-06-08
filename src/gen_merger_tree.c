// This code outputs the dark matter halo merger rates to text files.
// The merger rates are calculated using the fitting functions from
// the Appendix of Zhang et al. (2021) (Trinity Paper I).
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
  if (argc < 2) 
  {
    fprintf(stderr, "Usage: %s mass_cache >mass_function_file 2>merger_rate_file\n", argv[0]);
    exit(1);
  }

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
  load_mf_cache(argv[1]);
  // Initialize all the timesteps/snapshots.
  init_timesteps();

  // The output files for halo mass functions (f_mf) and merger rates (f_mr)
  FILE *f_mf = fopen("/media/hzhang/DataBackup/UM_DR1/Trinity/HaloMF.dat", "w");
  FILE *f_mr = fopen("/media/hzhang/DataBackup/UM_DR1/Trinity/HaloMR.dat", "w");
  fprintf(f_mf, "#num_snapshot scale nd[cMpc^(-3)/bin], BPDEX=%d\n", BPDEX);
  fprintf(f_mr, "#num_snapshot scale mr[cMpc^(-3)/bin^2], BPDEX=%d\n", BPDEX);
  // Output the mass functions and merger rates.
  for (i=0; i<num_outputs; i++) 
  {
    fprintf(f_mf, "%d %.6f", i, steps[i].scale);
    fprintf(f_mr, "%d %.6f", i, steps[i].scale);
    for (j=0; j<M_BINS; j++)
    {
      fprintf(f_mf, " %.6e", steps[i].t[j]);
      // The merger rates are M_BINS*M_BINS matrices, with the (i, j)
      // entry being the merger rates between halos in the i-th and
      // j-th bins.
      for (k=0; k<M_BINS; k++)
        fprintf(f_mr, " %.6e", steps[i].merged[k*M_BINS+j]);
    }
    fprintf(f_mf, "\n");
    fprintf(f_mr, "\n");
  }
  fclose(f_mf);
  fclose(f_mr);
  return 0;
}
