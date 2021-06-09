// Calculate the halo mass distribution for a certain quasar population between (z_low, z_high), and
// with BH mass and bolometric luminosity above bh_mass_low and bh_Lbol_low, respectively.
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
  int64_t i, j, k, l;
  struct smf_fit smf;
  
  if (argc < 7) 
  {
    fprintf(stderr, "Usage: %s mass_cache param_file z_low z_high bh_mass_low(in Msun) bh_Lbol_low(in erg/s) (> output_file)\n", argv[0]);
    exit(1);
  }

  // Read in model parameters
  FILE *param_input = check_fopen(argv[2], "r");
  char buffer[2048];
  fgets(buffer, 2048, param_input);
  read_params(buffer, smf.params, NUM_PARAMS);

  double z_low = atof(argv[3]);
  double z_high = atof(argv[4]);
  double Mbh_low = atof(argv[5]);
  double Lbol_low = atof(argv[6]);

  // Fix some model parameters.
  assert_model(&smf); 
  // We use non-linear scaling relation between the radiative and total Eddington ratios.
  nonlinear_luminosity = 1;
  // Turn off the built-in GSL error handler that kills the program
  // when an error occurs. We handle the errors manually.
  gsl_set_error_handler_off();
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
  printf("#z_low=%.2f, z_high=%.2f, Mbh_low=%.6f, Lbol_low=%.6f\n", z_low, z_high, Mbh_low, Lbol_low);
  printf("#Note: prob = Const.*prob_Mbh*prob_eta*nd_halo, where the constant is chosen so that the summation of prob equals unity.\n");
  printf("#Mh prob prob_Mbh prob_eta nd_halo\n");
  double t,m;

  // Calculate the #'s of snapshots that are the closest to the input lower and upper redshifts.
  int64_t step_low, step_high;
  double f;
  calc_step_at_z(z_low, &step_high, &f);
  calc_step_at_z(z_high, &step_low, &f);

  double prob_Mh[M_BINS] = {0};
  
  // Calculate the probabilities integrated over
  for (i=step_low; i<=step_high; i++)
  {
    double sm_scatter = steps[i].smhm.scatter * steps[i].smhm.bh_gamma;
    double scatter = sqrt(sm_scatter*sm_scatter 
                        + steps[i].smhm.bh_scatter*steps[i].smhm.bh_scatter);
    double mbh_min = steps[i].bh_mass_min, mbh_max = steps[i].bh_mass_max;
    double mbh_inv_bpdex = (mbh_max - mbh_min) / MBH_BINS;

    double prob_lbol[MBH_BINS] = {0}; //The cumulative probabilities of being more luminous than Lbol_low,
                                      //as a function of black hole mass.
    // Count the probabilities for each SMBH bin.
    for (k=0; k<MBH_BINS; k++)
    {
      double Mbh = mbh_min + (k + 0.5) * mbh_inv_bpdex;
      //fprintf(stderr, "i=%d, z=%f, mbh=%f\n", i, 1/steps[i].scale-1, Mbh);
      if (Mbh < Mbh_low) continue;
      double lbol_f = (Lbol_low - LBOL_MIN) * LBOL_BPDEX;
      int64_t lbol_b = lbol_f;
      lbol_f -= lbol_b;
      if (lbol_b >= LBOL_BINS - 1) {lbol_b = LBOL_BINS - 2; lbol_f = 1;}
      prob_lbol[k] = (1 - lbol_f) * steps[i].lum_dist_full[k*LBOL_BINS + lbol_b];
      for (l=lbol_b+1; l<LBOL_BINS; l++) prob_lbol[k] += steps[i].lum_dist_full[k*LBOL_BINS + l];
    }


    // Allocate each SMBH bin's contribution to halo mass bins, to account for halo number densities.
    for (j=0; j<M_BINS; j++)
    {
      for (k=0; k<MBH_BINS; k++)
      {
        double Mbh = mbh_min + (k + 0.5) * mbh_inv_bpdex;
        if (Mbh < Mbh_low) continue;
        double dMbh = Mbh - steps[i].log_bh_mass[j];
        double prob_Mbh = 1 / (sqrt(2*M_PI) * scatter) * exp(-dMbh*dMbh / (2*scatter*scatter)) * mbh_inv_bpdex; //dimension: dlogMbh^{-1}
        prob_Mh[j] += prob_Mbh * prob_lbol[k];
      }
      prob_Mh[j] *= steps[i].t[j];
    }
  }

  // Output
  for (i = 0; i < M_BINS; i++) 
  {
    printf("%.1f %.6e\n", (M_MIN + (i + 0.5) * INV_BPDEX), prob_Mh[i]);
  }
  
  return 0;
}
