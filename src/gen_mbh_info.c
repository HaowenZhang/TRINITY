// Print out the median SMBH masses, the slope and scatter of the SMBH mass--bulge mass
// relation, and the scatter around the stellar mass--halo mass relation.
// These information will be used to calculate the number density of very massive SMBHs
// at (or up to) a given redshift.
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
  int64_t i, j;
  struct smf_fit smf;
  if (argc<2+NUM_PARAMS) 
  {
    fprintf(stderr, "Usage: %s mass_cache (mcmc output)\n", argv[0]);
    exit(1);
  }
  for (i=0; i<NUM_PARAMS; i++)
    smf.params[i] = atof(argv[i+2]);

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

  printf("#Is the model invalid? %e\n", INVALID(smf));
  double t,m;

  printf("#1+z M_h SM obs_uv sfr  M_bh bhar bh_eta BH_Merge_Rate bh_unmerged dt\n");


  for (i=0; i<num_outputs; i++) 
  {
    double dt = steps[i].dt; 
    double zp1 = 1.0/steps[i].scale;
    double mu = steps[i].smhm.mu;
    double z = zp1 - 1;
    for (j=0; j<M_BINS; j++) 
    {
      m = M_MIN + (j + 0.5) * INV_BPDEX;  
      double log_bh_mass, bh_acc_rate, bh_merge_rate, bh_unmerged, sfr, bh_eta, obs_uv, log_sm;
      log_bh_mass = steps[i].log_bh_mass[j];
      bh_acc_rate = steps[i].bh_acc_rate[j];
      bh_merge_rate = steps[i].bh_merge_rate[j];
      bh_unmerged = steps[i].bh_unmerged[j];
      sfr = steps[i].sfr[j];
      bh_eta = steps[i].bh_eta[j];
      obs_uv = steps[i].obs_uv[j];
      log_sm = steps[i].log_sm[j];
      double nd = steps[i].t[j];
      printf("%f %f %f %f %e %f %e %e %e %e %e %f %f %f %e\n", zp1, m, log_sm, obs_uv, sfr, log_bh_mass, bh_acc_rate, bh_eta, bh_merge_rate, bh_unmerged, dt, steps[i].smhm.scatter, steps[i].smhm.bh_gamma, steps[i].smhm.bh_scatter, nd);
    }
  }
  return 0;
}
