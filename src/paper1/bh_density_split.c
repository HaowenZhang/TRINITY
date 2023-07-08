// Generate cosmic black hole mass densities belw and above 10^8
// Msuns as functions of redshift. Ref: Fig. 25 of Zhang et al. (2021).
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

  double sqrt12 = sqrt(0.5);

  printf("#z rho_BH(Mbh < 10^8 Msun) rho_BH(Mbh > 10^8 Msun)\n");
  // The (log10 of) SMBH masses defining the bins.
  double mbh_grid[] = {0, 8, 18};
  // bins in SMBH mass per dex, when calculating the contribution from
  // each halo mass bin to the SMBH mass density.
  int64_t mbh_bpdex = 20;
  double mbh_inv_bpdex = 1.0 / mbh_bpdex;
  for (i=0; i < num_outputs; i++)
  {
    double z = 1 / steps[i].scale - 1;
    double bh_unmerged = 0;
    // The scatter in SMBH mass at fixed halo mass. This is a quadratic summation
    // of the scatter around the black hole mass--bulge mass relation and the
    // scatter around the stellar mass--halo mass relation, scaled by the slope
    // of the black hole mass--bulge mass relation.
    double scatter = sqrt(steps[i].smhm.bh_scatter*steps[i].smhm.bh_scatter+
                          steps[i].smhm.scatter*steps[i].smhm.scatter
                          *steps[i].smhm.bh_gamma*steps[i].smhm.bh_gamma);
    double inv_scatter = 1.0 / scatter;
    // Pre-calculate the normalization of Gaussian distributions. 
    double norm_gauss = sqrt12 / scatter / sqrt(M_PI);
    // SMBH mass densities in different bins.
    double rho_bh_split[2] = {0};
    // Iterate over each halo mass bin
    for (int j = 0; j < M_BINS; j++) 
    {
      // Two broad SMBH mass bins
      for (int k=0; k<2; k++)
      {
        double mbh_low = mbh_grid[k];
        double mbh_high = mbh_grid[k+1];
        int64_t mbh_bins = (mbh_high - mbh_low) * mbh_bpdex;
        // Divide each broad SMBH mass bin into finer sub-bins
        for (int l=0; l<mbh_bins; l++)
        {
          // Calculate the corresponding SMBH mass, and the contribution
          // from the j-th halo mass bin to this SMBH mass bin.
          double mbh = mbh_low + (l + 0.5) * mbh_inv_bpdex;
          double prob = norm_gauss * exp(-0.5 * (mbh - steps[i].log_bh_mass[j]) *
                                                (mbh - steps[i].log_bh_mass[j]) *
                                                inv_scatter * inv_scatter) * mbh_inv_bpdex;
          rho_bh_split[k] += exp10(mbh) * prob * steps[i].t[j];

        }
      }
    }
    // Print the results.
    printf("%f", z);
    for (int k=0; k<2; k++) printf(" %e", rho_bh_split[k]);
    printf("\n");
  }
  return 0;
}
