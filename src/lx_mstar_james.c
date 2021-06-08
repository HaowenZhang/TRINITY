// Calculate average AGN X-ray luminosities as functions of galaxy
// stellar mass at a given redshift, z.
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

#define MASS_START 7
#define MASS_STOP 12.5
#define MASS_BPDEX 10
#define MASS_BINS (int)((MASS_STOP-MASS_START)*(MASS_BPDEX))
#define MASS_STEP (1.0/(double)MASS_BPDEX)


int main(int argc, char **argv)
{
  int64_t i, j, k, l;
  struct smf_fit smf;

  if (argc<2+NUM_PARAMS) 
  {
    fprintf(stderr, "Usage: %s mass_cache (mcmc output) z\n", argv[0]);
    exit(1);
  }
  // Read in model parameters and redshift.
  for (i=0; i<NUM_PARAMS; i++)
    smf.params[i] = atof(argv[i+2]);
  double z = atof(argv[i+4]);

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
  // Fix some model parameters.
  assert_model(&smf); 

  printf("#Is the model invalid? %e\n", INVALID(smf));
  printf("#z=%.2f\n", z);
  printf("#Mbh Phi(Mbh)\n");
  double t,m;
  int64_t step;
  double f;
  // Calculate the # of snapshot that is the closest to the input redshift.
  calc_step_at_z(z, &step, &f);

  double lx_avg[MBH_BINS] = {0}; //The probabilities of having luminosities between (Lbol_low, Lbol_high), as a function of BH mass.
  double lx_at_lbol[LBOL_BINS] = {0}; //The corresponding Lx at each Lbol bin center.
  for (i=0; i<LBOL_BINS; i++)
  {
    double lbol_tmp = LBOL_MIN + (i + 0.5) * LBOL_INV_BPDEX;
    // use a constant bolometric correction kbol = 25 to be consistent with Aird et al. (2018).
    lx_at_lbol[i] = lbol_tmp - log10(25);
  }

  // calculate the total scatter in BH mass at fixed ***halo mass***,
  // which is a quadratic sum of the scatter around the median black
  // hole mass--bulge mass relation, and that of the median stellar
  // mass--halo mass relation, enhanced by the slope of the black
  // hole mass--bulge mass relation.
  double sm_scatter = steps[step].smhm.scatter * steps[step].smhm.bh_gamma;
  double scatter = sqrt(sm_scatter*sm_scatter 
                      + steps[step].smhm.bh_scatter*steps[step].smhm.bh_scatter);
  double mbh_min = steps[step].bh_mass_min, mbh_max = steps[step].bh_mass_max;
  double mbh_inv_bpdex = (mbh_max - mbh_min) / MBH_BINS;
  double mu = steps[step].smhm.mu;

  // Calculate the average Lx as a function of BH mass using
  // luminosity distributions.
  for (i=0; i<MBH_BINS; i++)
  {
    double mbh = mbh_min + (i + 0.5) * mbh_inv_bpdex;
    double norm = 0;

    for (j=0; j<LBOL_BINS; j++) 
    {
      lx_avg[i] += exp10(lx_at_lbol[j]) * steps[step].lum_dist_full[i*LBOL_BINS+j];
      norm += steps[step].lum_dist_full[i*LBOL_BINS+j];
    }
    lx_avg[i] /= norm;

    double dc = steps[step].smhm.bh_duty;
    double f_mass = exp((mbh - steps[step].smhm.dc_mbh) / steps[step].smhm.dc_mbh_w);
    f_mass = f_mass / (1 + f_mass);
    dc *= f_mass;
    if (dc < 1e-4) dc = 1e-4;
    lx_avg[i] *= dc;
    // fprintf(stderr, "%f %e\n", mbh, prob_lbol[i]);
  }

  // Count the contribution to different stellar mass bins from different SMBH masses.
  for (i=0; i<MASS_BINS; i++) 
  {
    m = MASS_START + i*MASS_STEP; //observed stellar mass.
    // double sfr = calc_sfr_mstar(m, z);
    double mb = bulge_mass(m + mu, 1.0 / (1 + z));
    double mbh_med = calc_bh_at_bm(mb, steps[step].smhm);
    double lx = 0;

    for (j=0; j<MBH_BINS; j++)
    {
      double mbh = mbh_min + (j + 0.5) * mbh_inv_bpdex;
      double dmbh = (mbh - mbh_med) / steps[step].smhm.bh_scatter;
      // The probability of having a BH mass.
      double prob_tmp = 1 / sqrt(2*M_PI) / steps[step].smhm.bh_scatter * exp(-0.5*dmbh*dmbh) * mbh_inv_bpdex;
      lx += prob_tmp * lx_avg[j];
    }

    printf("%f %e\n", m, lx);

  }

  
  return 0;
}
