// Generate cosmic black hole mass densities (including the total density and
// densities of SMBHs above certain luminosity and Eddington ratio limits as 
// functions of redshift.
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include "../base/observations.h"
#include "../base/smf.h"
#include "../base/all_smf.h"
#include "../base/distance.h"
#include "../base/integrate.h"
#include "../base/mlist.h"
#include "../base/calc_sfh.h"
#include "../base/expcache2.h"

int main(int argc, char **argv)
{
  int64_t i;
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
  // Set up the PSF for stellar mass functions. See observations.c.
  setup_psf(1);
  // Load cached halo mass functions.
  load_mf_cache(argv[1]);
  // Initialize all the timesteps/snapshots.
  init_timesteps();
  INVALID(smf) = 0;
  // Calculate the star-formation histories and black hole histories. See calc_sfh.c.
  calc_sfh(&smf);

  printf("#z Total Active(Ledd>0.01) Active(Lbol>43) Active(Lbol>46) unmerged dt\n");


  for (i=0; i < num_outputs; i++)
  {
    double z = 1 / steps[i].scale - 1;
    double bh_unmerged = 0;
    // Count the cosmic unmerged/wandering black hole mass density.
    for (int j = 0; j < M_BINS; j++) 
      if (steps[i].bh_unmerged[j] > 0)
      {
        bh_unmerged += steps[i].bh_unmerged[j] * steps[i].t[j];
      }
    // Calculate the cosmic BH mass densities, including: the total BH mass density,
      // the BH mass densities of BHs with Eddington ratios above 10^-2,
      // with bolometric luminosities above 10^43 erg/s, and above 10^46 erg/s.
    fprintf(stdout, "%f %g %g %g %g %g %g", z, log10(cosmic_bh_density(z,0,0,NULL)), 
     log10(cosmic_bh_density(z,-2,0,NULL)), log10(cosmic_bh_density(z,0,43,NULL)), log10(cosmic_bh_density(z,0,46,NULL)), log10(bh_unmerged), steps[i].dt);
    // Calculate the BH mass density contributed by different halo mass bins.
    for (double m=7.0; m < 16.0; m += 1.0) fprintf(stdout, " %g", log10(cosmic_bh_density_split(z, m, m + 1.0, NULL)));
    fprintf(stdout, "\n");
  }

  // Remove the contribution from SMBH mergers and do the same calculation over again.
  remove_bh_merger();
  for (i=0; i < num_outputs; i++)
  {
    double z = 1 / steps[i].scale - 1;
    double bh_unmerged = 0;
    for (int j = 0; j < M_BINS; j++) 
      if (steps[i].bh_unmerged[j] > 0)
      {
        bh_unmerged += steps[i].bh_unmerged[j] * steps[i].t[j];
      }
      
    fprintf(stderr, "%f %g %g %g %g %g %g", z, log10(cosmic_bh_density(z,0,0,NULL)), 
     log10(cosmic_bh_density(z,-2,0,NULL)), log10(cosmic_bh_density(z,0,43,NULL)), log10(cosmic_bh_density(z,0,46,NULL)), log10(bh_unmerged), steps[i].dt);
    for (double m=7.0; m < 16.0; m += 1.0) fprintf(stderr, " %g", log10(cosmic_bh_density_split(z, m, m + 1.0, NULL)));
    fprintf(stderr, "\n");
  }

  return 0;
}
