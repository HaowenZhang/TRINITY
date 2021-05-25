// Calculate the mass probability distribution of galaxies hosting a quasar 
// with given SMBH mass and luminosity at a given redshift. 
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
  gsl_set_error_handler_off();
  if (argc<2+NUM_PARAMS) {
    fprintf(stderr, "Usage: %s mass_cache (mcmc output) z bh_mass(in Msun) bh_Lbol(in erg/s)\n", argv[0]);
    exit(1);
  }
  for (i=0; i<NUM_PARAMS; i++)
    smf.params[i] = atof(argv[i+2]);
  double z = atof(argv[i+4]);
  double Mbh = atof(argv[i+5]);
  double Lbol = atof(argv[i+6]);
  // Turn off the built-in GSL error handler that kills the program
  // when an error occurs. We handle the errors manually.
  gsl_set_error_handler_off();
  // We use the non-linear scaling relation between the total and radiative
  // Eddington ratios.
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
  printf("#z=%.2f, Mbh=%.6f, Lbol=%.6f\n", z, Mbh, Lbol);
  printf("#Mstar prob prob_Mbh\n");
  double t,m;
  int64_t step;
  double f;
  // Calculate the # of the snapshot that is the closest to
  // the input redshift.
  calc_step_at_z(z, &step, &f);
  // Scatters around the stellar mass--halo mass relation
  // and the black hole mass--bulge mass relation.
  double sm_scatter = steps[step].smhm.scatter;
  double bh_scatter = steps[step].smhm.bh_scatter;
  printf("#bh_scatter: %f\n", bh_scatter);

  double total = 0;
  // Calculate the Eddington ratio given the input luminosity
  // and BH mass. This is because we calculate and store 
  // Eddington ratio distributions.
  double eta = Lbol - 38.1 - Mbh;
  // Define the stellar mass bins where we calculate the 
  // probabilities.
  double sm_min = 8; double sm_max = 12;
  int sm_bpdex = 10;
  int sm_bins = (int)((sm_max - sm_min) * sm_bpdex);
  double sm_inv_bpdex = 1.0 /sm_bpdex;

  // Allocate the space to save probabilities.
  double *probs, *probs_eta, *probs_Mbh, *probs_Mstar;
  probs = malloc(sizeof(double) * sm_bins);
  memset(probs, 0, sizeof(double) * sm_bins);
  probs_eta = malloc(sizeof(double) * sm_bins);
  memset(probs_eta, 0, sizeof(double) * sm_bins);
  probs_Mbh = malloc(sizeof(double) * sm_bins);
  memset(probs_Mbh, 0, sizeof(double) * sm_bins);
  probs_Mstar = malloc(sizeof(double) * sm_bins);
  memset(probs_Mstar, 0, sizeof(double) * sm_bins);

  // iterate over each stellar mass bin.
  for (i = 0; i < sm_bins; i++)
  {
    // Calculate the median BH mass given the best-fitting black
    // hole mass--bulge mass relation.
    double sm = (sm_min + (i + 0.5) * sm_inv_bpdex);
    double bm = bulge_mass(sm + steps[step].smhm.mu, steps[step].scale);
    double log_bh_mass = calc_bh_at_bm(bm, steps[step].smhm);
    // Calculate the probability of the galaxies' having a SMBH with input
    // mass, given the median BH mass.
    double dMbh = Mbh - log_bh_mass;
    double prob_Mbh = 1 / (sqrt(2*M_PI) * bh_scatter) * exp(-dMbh*dMbh / (2*bh_scatter*bh_scatter));
    
    // Calculate the contribution to this stellar mass bin from
    // every halo mass bin and accumulate them.
    for (j = 0; j < M_BINS; j++)
    {
      // Calculate the probability of the halos' having the given
      // stellar mass, given the median stellar mass in this halo
      // mass bin.
      double dMstar = sm - steps[step].log_sm[j];
      double prob_Mstar = 1 / (sqrt(2*M_PI) * sm_scatter) * exp(-dMstar*dMstar / (2*sm_scatter*sm_scatter));
      // And calculate the probability of these halos' having black holes
      // with the needed Eddington ratio.
      double eta_frac = eta - steps[step].bh_eta[j];
      double prob_eta;
      if (eta_frac < steps[step].ledd_min[j] || eta_frac > steps[step].ledd_max[j]) 
        prob_eta = 0; 
      else
      {
        double bher_f = (eta_frac-steps[step].ledd_min[j])*steps[step].ledd_bpdex[j];
        int64_t bher_b = bher_f;
        bher_f -= bher_b;
        double p1 = steps[step].bher_dist[j*BHER_BINS+bher_b];
        double p2 = steps[step].bher_dist[j*BHER_BINS+bher_b+1];
        if (bher_b >= BHER_BINS-1) p2 = p1;
        prob_eta = p1 + bher_f * (p2 - p1);
      }
      // Add up their contributions.
      probs[i] += prob_Mstar * prob_eta * steps[step].t[j];
    }

    probs[i] *= prob_Mbh;
    probs_Mbh[i] = prob_Mbh;
    total += probs[i];
  }

  //normalize and convert it to dex^-1.
  if (total > 0)
  {
    for (i = 0; i < sm_bins; i++) probs[i] /= (total*sm_inv_bpdex); 
  }
  
  // Print them out.
  for (i = 0; i < sm_bins; i++) 
  {
    printf("%.3f %.6e %.6e\n", (sm_min + (i + 0.5) * sm_inv_bpdex), probs[i], probs_Mbh[i]);
  }
  
  return 0;
}
