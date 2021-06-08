// Calculate the accretion-rate-weighted average AGN efficiency at each snapshot (redshift).
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
  
  if (argc < 3) 
  {
    fprintf(stderr, "Usage: %s mass_cache parameter_file (> output_file)\n", argv[0]);
    exit(1);
  }

  // Read in model parameters
  FILE *param_input = check_fopen(argv[2], "r");
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
  load_mf_cache(argv[1]);
  // Initialize all the timesteps/snapshots.
  init_timesteps();
  INVALID(smf) = 0;
  // Calculate the star-formation histories and black hole histories. See calc_sfh.c.
  calc_sfh(&smf);

  printf("#sn_id z total_accretion eff_rad\n");

  double weight = 0;
  double tot = 0;
  double bh_density = 0;
  double bh_unmerged = 0;

  for (i=0; i<num_outputs-1; i++) 
  {

    double eff_rad = steps[i].smhm.bh_efficiency_rad;
    double dt = steps[i].dt;
    double acc_tot = 0;


    for (int64_t j = 0; j < M_BINS; j++) 
    {
      double m;

      double bh_acc_rate = steps[i].bh_acc_rate[j] == 1e-8 ? 0 : steps[i].bh_acc_rate[j];
      acc_tot += bh_acc_rate * steps[i].t[j] * dt;
    }


    // Weigh the efficiency values by cosmic SMBH accretion rates.
    weight += acc_tot;
    tot += acc_tot * eff_rad;

    printf("%d %f %f %f\n", i, 1 / steps[i].scale - 1, acc_tot, eff_rad);

    if (i == num_outputs-2)
    {
      for (int64_t j = 0; j < M_BINS; j++) 
      {
        bh_density += steps[i].bh_mass_avg[j] * steps[i].t[j];
        bh_unmerged += steps[i].bh_unmerged[j] * steps[i].t[j];
      }
    }

  }

  printf("The mass averaged radiative efficiency is: %.6f\n", tot / weight);
  printf("The cosmic central BH density: : %.3e\n", bh_density);
  printf("The cosmic unmerged BH density: : %.3e\n", bh_unmerged);
  return 0;
}
