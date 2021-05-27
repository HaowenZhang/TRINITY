// Print out the SMBH bolometric luminosity distributions
// at a given redshift z, as functions of ***observed***
// SMBH masses. By "observed", we mean the single-epoch
// spectroscopic SMBH mass measurements.
// The user needs to provide a random scatter in the 
// observed SMBH masses around intrinsic ones in this calculation. 
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
  if (argc<3+NUM_PARAMS) 
  {
    fprintf(stderr, "Usage: %s z scatter_obs mass_cache (mcmc output)\n", argv[0]);
    exit(1);
  }
  // Read in the redshift, observed scatter, and the model parameters.
  double z = atof(argv[1]);
  double scatter_obs = atof(argv[2]);
  for (i=0; i<NUM_PARAMS; i++)
    smf.params[i] = atof(argv[i+4]);
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

  // Calculate the SMBH ***intrinsic*** mass function
  double bhmf_z_int[MBH_BINS] = {0};
  for (i=0; i<MBH_BINS; i++)
  {
    m = mbh_min + (i + 0.5) * mbh_inv_bpdex;
    bhmf_z_int[i] = calc_bhmf(m, z);
  }

  // Calculate the distribution of SMBH ***observed*** masses.
  double P_Mo[MBH_BINS] = {0};//P(Mbh_obs) = Integral(P(Mbh_obs|Mbh_int)*P(Mbh_int));
  for (i=0; i<MBH_BINS; i++)
  {
    // double P_Mo = 0; 
    double Mo = mbh_min + (i + 0.5) * mbh_inv_bpdex;
    for (int j=0; j<MBH_BINS; j++)
    {
      double Mi = mbh_min + (j + 0.5) * mbh_inv_bpdex;
      double delta = (Mo - Mi) / scatter_obs;
      P_Mo[i] += exp(-0.5 * delta * delta) * bhmf_z_int[j];
    }
  }

  // Calculate the probability that the SMBH in j-th intrinsic mass bin
  // falls in the i-th observed mass bin. 
  double P_Mi_Mo[MBH_BINS*MBH_BINS] = {0};
  for (i=0; i<MBH_BINS; i++) //observed Mbh
  {
    double Mo = mbh_min + (i + 0.5) * mbh_inv_bpdex;
    for (int j=0; j<MBH_BINS; j++) //intrinsic Mbh
    {
      double Mi = mbh_min + (j + 0.5) * mbh_inv_bpdex;
      double delta = (Mo - Mi) / scatter_obs;
      P_Mi_Mo[i*MBH_BINS+j] = exp(-0.5 * delta * delta) * bhmf_z_int[j] / P_Mo[i];
    }
  }

  // Then the luminosity distribution for each observed mass bin
  // is the weighted sum of that for each intrinsic mass bin, weighted
  // by their contribution to the observed mass bin.
  double P_L_Mo[LBOL_BINS*MBH_BINS] = {0};
  for (i=0; i<LBOL_BINS; i++) //luminosity
  {
    for (int j=0; j<MBH_BINS; j++) //observed Mbh
    {
      for (int k=0; k<MBH_BINS; k++) //intrinsic Mbh
        P_L_Mo[i*MBH_BINS+j] += steps[step].lum_dist_full[k*LBOL_BINS+i] * P_Mi_Mo[j*MBH_BINS+k];
    }
  }

  // Print them out.
  for (i=0; i<MBH_BINS; i++)
  {
    if (mbh_min + (i + 0.5) * mbh_inv_bpdex < 8) continue;
    fprintf(stdout, "%.6f ", mbh_min + (i + 0.5) * mbh_inv_bpdex);
    for (int j=0; j<LBOL_BINS; j++)
      fprintf(stdout, "%.6e ", P_L_Mo[j*MBH_BINS+i]);
    fprintf(stdout, "\n");
  }
  return 0;
}
