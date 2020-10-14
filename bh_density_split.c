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
  if (argc<2+NUM_PARAMS) {
    fprintf(stderr, "Usage: %s mass_cache (mcmc output)\n", argv[0]);
    exit(1);
  }
  for (i=0; i<NUM_PARAMS; i++)
    smf.params[i] = atof(argv[i+2]);
  gsl_set_error_handler_off();
  setup_psf(1);
  load_mf_cache(argv[1]);
  init_timesteps();
  INVALID(smf) = 0;
  calc_sfh(&smf);

  double sqrt12 = sqrt(0.5);
  //int64_t mbh_bins = 1000;

  printf("#z rho_5m6 rho_6m7 rho_7m8 rho_8m9 rho_9m10\n");
  // for (i=0; i<200; i++) {
  //   double z = i*0.05;
  //   printf("%f %g %g %g %g %g\n", z, log10(cosmic_bh_density(z,0,0,NULL)), log10(cosmic_old_bh_density(z,0,0,NULL)), 
	 //   log10(cosmic_bh_density(z,-2,0,NULL)), log10(cosmic_bh_density(z,0,43,NULL)), log10(cosmic_bh_density(z,0,46,NULL)));
  // }
  double mbh_grid[] = {0, 8, 18};
  //double mbh_grid[] = {5,6,7,8,9,10};
  int64_t mbh_bpdex = 20;
  //int64_t mbh_bins = (mbh_grid[2] - mbh_grid[0]) * mbh_bpdex;
  double mbh_inv_bpdex = 1.0 / mbh_bpdex;
  for (i=0; i < num_outputs; i++)
  {
    double z = 1 / steps[i].scale - 1;
    double bh_unmerged = 0;
    double scatter = sqrt(steps[i].smhm.bh_scatter*steps[i].smhm.bh_scatter+
                          steps[i].smhm.scatter*steps[i].smhm.scatter
                          *steps[i].smhm.bh_gamma*steps[i].smhm.bh_gamma);
    double inv_scatter = 1.0 / scatter;
    double norm_gauss = sqrt12 / scatter / sqrt(M_PI);
    //double rho_bh_split[5] = {0};
    double rho_bh_split[2] = {0};
    for (int j = 0; j < M_BINS; j++) 
    {
      for (int k=0; k<2; k++)
      //for (int k=0; k<5; k++)
      {
        // double log_bh_mass_avg = log10(steps[i].bh_mass_avg[j]);
        double mbh_low = mbh_grid[k];
        double mbh_high = mbh_grid[k+1];
        // double sigma_low = (mbh_low - steps[i].log_bh_mass[j]) / scatter;
        // double sigma_high = (mbh_high - steps[i].log_bh_mass[j]) / scatter;
        // double frac_between = 0.5 * (erf(sigma_high*sqrt12) - erf(sigma_low*sqrt12));
        int64_t mbh_bins = (mbh_high - mbh_low) * mbh_bpdex;
        for (int l=0; l<mbh_bins; l++)
        {
          double mbh = mbh_low + (l + 0.5) * mbh_inv_bpdex;
          double prob = norm_gauss * exp(-0.5 * (mbh - steps[i].log_bh_mass[j]) *
                                                (mbh - steps[i].log_bh_mass[j]) *
                                                inv_scatter * inv_scatter) * mbh_inv_bpdex;
          rho_bh_split[k] += exp10(mbh) * prob * steps[i].t[j];

        }
      }
    }
    printf("%f", z);
    for (int k=0; k<2; k++) printf(" %e", rho_bh_split[k]);
    //for (int k=0; k<5; k++) printf(" %e", rho_bh_split[k]);
    printf("\n");
  }
  return 0;
}
