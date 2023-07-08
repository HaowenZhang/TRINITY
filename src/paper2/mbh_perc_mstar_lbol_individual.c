#include <stdlib.h>
#include <stdio.h>
#include <string.h>
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
#include "../base/param_default.h"

extern double param_default[];

// Calculate the 16th, 50th, and 84th percentiles of log10(SMBH mass)[Msun]
// as a function of galaxy stellar mass and redshift,
// for quasars ***around*** a bolometric luminosity threshold.
// Input parameters:
// lbol: The ***log10*** of the lower limit in bolometric AGN luminosity [erg/s]
// z: redshift
// sm: The ***log10*** of galaxy stellar mass [Msun]
// sigma_mbh: The scatter in the log10 of observed SMBH mass around
// the intrinsic mass, in dex (e.g., ~0.4-0.6 dex from virial estimates).
// perc: the array to store the calculated percentiles of SMBH mass.
void mbh_perc_mstar_lbol_individual(double lbol, double z, double sm, double sigma_mbh, double *perc)

{
  // Find out which snapshot the input redshift corresponds to.
  int64_t step; double sf = 0;
  calc_step_at_z(z, &(step), &sf);

  int64_t i, j, k;

  // Define the minimum and maximum SMBH masses, and thus the bin width in SMBH mass
  double mbh_min = steps[step].bh_mass_min; double mbh_max = steps[step].bh_mass_max;
  double mbh_inv_bpdex = (mbh_max - mbh_min) / MBH_BINS;

  // Intrinsic scatter of the BH mass--bulge mass relation
  double bh_scatter = steps[step].smhm.bh_scatter;
  // Intrinsic scatter of the stellar mass--halo mass relation
  double scatter = steps[step].smhm.scatter;
  // Systematic offset between the observed and the intrinsic galaxy mass
  double mu = steps[step].smhm.mu;

  // Normalization of the Gaussian distributions of SMBH mass at fixed galaxy mass
  double gauss_norm_bh = 1.0 / sqrt(2*M_PI) / bh_scatter;
  // Normalization of the Gaussian distributions of galaxy mass at fixed halo mass
  double gauss_norm_gal = 1.0 / sqrt(2*M_PI) / scatter;

  // Find out the luminosity bin number corresponding to the input lower limit
  // of bolometric luminosity
  double lbol_f = (lbol - LBOL_MIN) * LBOL_BPDEX;
  int64_t lbol_b = lbol_f;
  lbol_f -= lbol_b;
  
  // Calculate the median SMBH mass at the input galaxy mass, ***regardless***
  // of SMBH luminosity
  double bm = bulge_mass(sm + mu, steps[step].scale);
  double mbh_med = calc_bh_at_bm(bm, steps[step].smhm);

  // distribution of intrinsic SMBH masses of quasars
  double mbh_dist[MBH_BINS] = {0};
  // distribution of ***observed*** (i.e., with the Eddington bias induced by 
  // the scatter in the observed SMBH mass SMBH masses of quasars
  double mbh_dist_edd[MBH_BINS] = {0};
  
  // traverse each halo mass bin
  for (j=0; j<M_BINS; j++)
  {
    // Ignore the halo mass bins without meaningful number densities
    if (!steps[step].t[j]) continue;    

    // Calculate the probability of the galaxy with the input stellar mass
    // to be hosted by the halos in this halo mass bin
    double dsm = (sm - steps[step].log_sm[j]) / scatter;
    double prob_mstar = gauss_norm_gal * exp(-0.5*dsm*dsm);

    // We only include the active SMBHs, so the AGN duty cycle as a function
    // of halo mass is also included.
    double dc = steps[step].bh_duty[j];
    if (dc < 1e-4) dc = 1e-4;

    // traverse each SMBH mass bin.
    for (k=0; k<MBH_BINS; k++)
    {
      // For each SMBH mass bin, calculate the temporary SMBH mass
      double mbh = mbh_min + (k + 0.5) * mbh_inv_bpdex;
      double dmbh = (mbh - mbh_med) / bh_scatter;
      // Calculate the probability of the SMBH with the temporary SMBH mass
      // to be hosted by the halos in this halo mass bin
      double prob_mbh = gauss_norm_bh * exp(-0.5*dmbh*dmbh) * mbh_inv_bpdex;

      // Calculate the Eddington ratio corresponding to the luminosity limit,
      // ***scaled*** by the typical Eddington ratio of SMBHs in this halo mass
      // bin.
      double eta_frac = lbol - 38.1 - mbh - (steps[step].bh_eta[j] + log10(steps[step].bh_f_occ[j]) 
        - (steps[step].smhm.rho_bh - 1) * (mbh - steps[step].log_bh_mass[j] - steps[step].smhm.log_bh_scatter_corr));
      
      // Calculate the probability of SMBHs to be at the input luminosity (+/- 0.1 dex)
      // by integrating over the Eddington ratio distribution.
      double bher_f = (eta_frac - 0.1 - steps[step].ledd_min[j])*steps[step].ledd_bpdex[j];
      int64_t bher_b = bher_f;
      bher_f -= bher_b;
 

      double bher_f2 = (eta_frac + 0.1 - steps[step].ledd_min[j])*steps[step].ledd_bpdex[j];
      int64_t bher_b2 = bher_f2;
      bher_f2 -= bher_b2;

      if (bher_b < 0) {bher_b = 0; bher_f = 0;}
      else if (bher_b >= BHER_BINS - 1) {bher_b = BHER_BINS - 2; bher_f = 1;}
      
      if (bher_b2 < 0) {bher_b2 = 0; bher_f2 = 0;}
      else if (bher_b2 >= BHER_BINS - 1) {bher_b2 = BHER_BINS - 2; bher_f2 = 1;}

      double prob = (1 - bher_f) * steps[step].bher_dist[j*BHER_BINS+bher_b];
      for (i=bher_b+1; i<bher_b2; i++)  
        prob += steps[step].bher_dist[j*BHER_BINS+i];
      prob += bher_f2 * steps[step].bher_dist[j*BHER_BINS+bher_b2];

      // And normalize the probability
      prob /= steps[step].bher_dist_norm[j];

      // accumulate the contribution to each SMBH mass bin from every halo
      // mass bin.
      mbh_dist[k] += prob_mstar * prob_mbh * prob * steps[step].t[j] * dc;
    }

  }

  // If there is a scatter in observed SMBH mass, we convolve the intrinsic SMBH
  // mass distribution with a gaussian kernel. 
  if (sigma_mbh > 0)
  {
    double dmbh_max = 8 * sigma_mbh;
    int nbins_stencils = 2 * (dmbh_max * MBH_BPDEX) + 1;
    double *weights_edd = NULL;
    weights_edd = (double *)malloc(sizeof(double) * nbins_stencils);
    memset(weights_edd, 0, sizeof(double) * nbins_stencils);
    double gauss_norm_edd = 1.0 / sqrt(2*M_PI) / sigma_mbh;
    for (k=0; k<nbins_stencils; k++)
    {
      double dmbh = (k - nbins_stencils / 2) * MBH_INV_BPDEX / sigma_mbh;
      weights_edd[k] = gauss_norm_edd * exp(-0.5 * dmbh * dmbh) *  MBH_INV_BPDEX;
    }
    for (k=0; k<MBH_BINS; k++)
    {
      int emin = k - (dmbh_max * MBH_BPDEX) >= 0? k - (dmbh_max * MBH_BPDEX) : 0;
      int emax = k + (dmbh_max * MBH_BPDEX) < MBH_BINS ? k + (dmbh_max * MBH_BPDEX) : MBH_BINS - 1;
      for (i=emin; i<=emax; i++)
        mbh_dist_edd[i] += weights_edd[nbins_stencils / 2 - (k - i)] * mbh_dist[k];
    }
    for (k=0; k<MBH_BINS; k++) mbh_dist[k] = mbh_dist_edd[k];
    
  }
  
  // Calculate the 16th, 50th, and 84th percentiles of SMBH mass
  // with the calculated quasar SMBH mass distribution at the input
  // galaxy mass.
  double mbh_dist_norm = 0;
  double mbh_cdf[MBH_BINS] = {0};
  for (k=0; k<MBH_BINS; k++) mbh_dist_norm += mbh_dist[k];
  for (k=0; k<MBH_BINS; k++) mbh_dist[k] /= mbh_dist_norm;
  mbh_cdf[0] = mbh_dist[0];
  for (k=0; k<MBH_BINS-1; k++) mbh_cdf[k+1] = mbh_cdf[k] + mbh_dist[k+1];
  //fprintf(stderr, "#%f\n", sm);
  //for (k=0; k<MBH_BINS; k++) fprintf(stderr, "%f %f\n", MBH_MIN + (k + 0.5) * MBH_INV_BPDEX, mbh_cdf[k]);
  fprintf(stderr, "%f ", sm);
  for (k=0; k<MBH_BINS; k++) fprintf(stderr, "%e ", mbh_dist[k]);
  fprintf(stderr, "\n");
  //for (k=0; k<MBH_BINS; k++) mbh_dist[k] /= mbh_dist_norm;
  double thres[3] = {0.16, 0.50, 0.84};
  for (i=0; i<3; i++)
  {
    for (k=0; k<MBH_BINS-1; k++)
    {
      if (mbh_cdf[k] <= thres[i] && thres[i] <= mbh_cdf[k+1])
        break;
    }
    double mbh_f = (thres[i] - mbh_cdf[k]) / (mbh_cdf[k+1] - mbh_cdf[k]);
    perc[i] = MBH_MIN + (k + 1 + mbh_f) * MBH_INV_BPDEX;
  }

}

int main(int argc, char **argv)
{
  int64_t i;
  struct smf_fit smf;

  // Abort with error if at least one of the following files is not provided:
  // the cached halo mass function files (provided along with the code: mf_bolshoi_planck.dat)
  // and the quasar catalog in a text file (see the code below for the required format).
  // mcmc_output is the model parameter values based on which the calculations are gonna
  // be done. If not provided, the code will use the best-fitting parameter values from
  // TRINITY Paper I and its erratum: https://ui.adsabs.harvard.edu/abs/2023MNRAS.518.2123Z/abstract,
  // which are stored in ../base/param_default.h.
  if (argc < 3) {
    fprintf(stderr, "Usage: %s mass_cache quasar_catalog (mcmc output)\n", argv[0]);
    exit(1);
  }



  // Read in the model parameter values if provided by the user
  if (argc >= 3+NUM_PARAMS)
    for (i=0; i<NUM_PARAMS; i++)
      smf.params[i] = atof(argv[i+3]);
  // Otherwise use our default values
  else
    for (i=0; i<NUM_PARAMS; i++)
      smf.params[i] = param_default[i];

  // Assume non-linear scaling relation between the total (radiative+kinetic) Eddington ratio
  // and its radiative component (see Section 2.7 of Paper I)
  nonlinear_luminosity = 1;
  // We will handle the gsl errors manually, so we shut off the automatic handler.
  gsl_set_error_handler_off();
  // make a cache for the standard gaussian distribution
  setup_psf(1);

  // Read in cached halo mass functions
  load_mf_cache(argv[1]);
  FILE * input_file = fopen(argv[2], "r");

  // Initialize timesteps
  init_timesteps();
  INVALID(smf) = 0;

  // Calculate star formation histories (including SMBH histories)
  calc_sfh(&smf);


  char buffer[1024];
  float z, sm, lbol, mbh, sigma_mbh;
  fprintf(stdout, "#z logMstar[Msun] logMbh[Msun] logLbol[erg/s] dmbh[sigma]\n");
  
  // Read in the quasar catalog line by line.
  while (fgets(buffer, 1023, input_file))
  {
    // The required format for the quasar catalog:
    // z: redshift
    // mbh: log10 of SMBH mass
    // sm: log10 of galaxy mass
    // sigma_mbh: uncertainty in SMBH mass
    int64_t n = sscanf(buffer, "%f %f %f %f %f", &z, &mbh, &sm, &lbol, &sigma_mbh);

    // Skip the lines that do not contain quasars
    if (buffer[0] == '#') continue;
    if (n != 5) {fprintf(stderr, "Invalid line, skip.\n"); continue;}

    // Calculate the 16th, 50th, and 84th percentiles in BH mass as a function of galaxy mass
    // and lower limit of AGN bolometric luminosity
    double mbh_perc[3] = {0};
    mbh_perc_mstar_lbol_individual(lbol, z, sm, sigma_mbh, mbh_perc);

    // It turns out that the distance from the median to 16th and 84th percentiles are
    // nearly the same. So here we use the average distance as the 1-sigma scatter.
    double dmbh_in_sigma = 2 * (mbh - mbh_perc[1]) / (mbh_perc[2] - mbh_perc[0]);
    fprintf(stdout, "%f %f %f %f %f\n", z, mbh, sm, lbol, dmbh_in_sigma); 
   }
  fclose(input_file);
  return 0;
}

