// Calculate cosmic radiative and kinetic AGN energy densities
// as a function of redshift.
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

#define MSUN_YR_TO_ERG_S 5.6631e46

int main(int argc, char **argv)
{
  int64_t i, j, k;
  struct smf_fit smf;
  if (argc<2+NUM_PARAMS) 
  {
    fprintf(stderr, "Usage: %s mass_cache (mcmc output)\n", argv[0]);
    exit(1);
  }
  // Read in model parameters.
  for (i=0; i<NUM_PARAMS; i++)
    smf.params[i] = atof(argv[i+2]);
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
  INVALID(smf) = 0;
  // Calculate the star-formation histories and black hole histories. See calc_sfh.c.
  calc_sfh(&smf);

  printf("#Is the model invalid? %e\n", INVALID(smf));
  double t,m;

  printf("#z rho_tot rho_rad rho_rad_6m7, rho_rad_7m8, rho_rad_8m9, rho_rad_9m10, rho_kin\n");
  // SMBH mass bins where energy densities are counted.
  double mbh_lims[] = {6,7,8,9,10};

  for (i=0; i<num_outputs; i++) 
  {
    double z = 1.0 / steps[i].scale - 1;
    double rho_rad = 0, rho_tot = 0, rho_kin = 0;
    double rho_rad_split[4] = {0};
    double eff_rad = steps[i].smhm.bh_efficiency_rad;
    double mbh_min = steps[i].bh_mass_min, mbh_max = steps[i].bh_mass_max;
    double mbh_inv_bpdex = (mbh_max - mbh_min) / MBH_BINS;
    // calculate the total scatter in BH mass at fixed ***halo mass***,
    // which is a quadratic sum of the scatter around the median black
    // hole mass--bulge mass relation, and that of the median stellar
    // mass--halo mass relation, enhanced by the slope of the black
    // hole mass--bulge mass relation.
    double scatter_tot = sqrt(steps[i].smhm.bh_scatter*steps[i].smhm.bh_scatter+
			steps[i].smhm.scatter*steps[i].smhm.scatter*steps[i].smhm.bh_gamma*steps[i].smhm.bh_gamma);
    // The correction from median to average BH masses for log-normal distributions.
    double corr = exp(0.5 * pow((scatter_tot * M_LN10), 2));
    // The critical Eddington ratio below which radiative Eddington ratio is proportional
    // to (total Eddington ratio)^2.
    double eta_crit = steps[i].smhm.bh_eta_crit;

    // Count radiative and kinetic energy densities from all SMBHs, 
    // and radiative energy densities from SMBHs in each SMBH mass bin.
    for (j=0; j<MBH_BINS; j++) 
    {
      double rho_rad_mbh = 0;
      double mbh = mbh_min + (j + 0.5) * mbh_inv_bpdex;
      for (k=0; k<LBOL_BINS; k++)
      {
        double lbol = LBOL_MIN + (k + 0.5) * LBOL_INV_BPDEX;
        rho_rad_mbh += exp10(lbol) * steps[i].lum_func_full[j*LBOL_BINS + k] * LBOL_INV_BPDEX;

        double eta_rad = lbol - 38.1 - mbh;
        if (eta_rad < eta_crit)
        {
          double eta_kin = log10(exp10(0.5*(eta_crit + eta_rad)) - exp10(eta_rad));
          double lkin = eta_kin + mbh + 38.1;
          if (eta_kin < -6 || mbh < 8) continue;

          rho_kin += exp10(lkin) * steps[i].lum_func_full[j*LBOL_BINS + k] * LBOL_INV_BPDEX;
          fprintf(stderr, "scale=%f, mbh=%f, lbol=%f, eta_rad=%f, eta_kin=%f, lkin=%f, contribution=%e\n", steps[i].scale, mbh, lbol, eta_rad, eta_kin, lkin, exp10(lkin) * steps[i].lum_func_full[j*LBOL_BINS + k] * LBOL_INV_BPDEX);
	      }
      }
      rho_rad += rho_rad_mbh;
      for (k=0; k<4; k++)
      {
        if ((mbh >= mbh_lims[k]) && (mbh < mbh_lims[k+1])) rho_rad_split[k] += rho_rad_mbh;
      }

    }
    rho_tot = rho_kin + rho_rad;

    fprintf(stdout, "%f %e %e ", z, rho_tot, rho_rad);
    for (j=0; j<4; j++)
    {
      fprintf(stdout, "%e ", rho_rad_split[j]);
    }
    fprintf(stdout, "%e\n", rho_kin);
  }
  return 0;
}
