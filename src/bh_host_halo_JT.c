// Calculate the mass probability distribution of halos hosting a quasar 
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
  int64_t i;
  struct smf_fit smf;
  if (argc<2+NUM_PARAMS) {
    fprintf(stderr, "Usage: %s mass_cache (mcmc output) z bh_mass(in Msun) bh_Lbol(in erg/s)\n", argv[0]);
    exit(1);
  }
  for (i=0; i<NUM_PARAMS; i++)
    smf.params[i] = atof(argv[i+2]);
  double z = atof(argv[i+4]);
  double Mbh = atof(argv[i+5]);
  double Lbol = atof(argv[i+6]);

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
  printf("#z=%.2f, Mbh=%.6f, Lbol=%.6f\n", z, Mbh, Lbol);
  printf("#Note: prob = Const.*prob_Mbh*prob_eta*nd_halo, where the constant is chosen so that the summation of prob equals unity.\n");
  printf("#Mh prob prob_Mbh prob_eta nd_halo\n");
  double t,m;
  int64_t step;
  double f;
  // Calculate the # of the snapshot that is the closest to
  // the input redshift.
  calc_step_at_z(z, &step, &f);
  // calculate the total scatter in BH mass at fixed ***halo mass***,
  // which is a quadratic sum of the scatter around the median black
  // hole mass--bulge mass relation, and that of the median stellar
  // mass--halo mass relation, enhanced by the slope of the black
  // hole mass--bulge mass relation.
  double sm_scatter = steps[step].smhm.scatter * steps[step].smhm.gamma;
  double scatter = sqrt(sm_scatter*sm_scatter 
                      + steps[step].smhm.bh_scatter*steps[step].smhm.bh_scatter);
  double probs[M_BINS] = {0};
  double probs_eta[M_BINS] = {0}; //The probability of halos' having a SMBH with the required
                                  //Eddington ratio.
  double probs_Mbh[M_BINS] = {0}; //The probability of halos' having a SMBH with the required
                                  //mass.
  double total = 0;
  double eta = Lbol - 38.1 - Mbh; //The required Eddington ratio given the luminosity and SMBH mass.
  
  for (i = 0; i < M_BINS; i++)
  {
    // Calculate the probability of halos' having the required SMBH mass, which is simply a
    // log-normal distribution.
    double dMbh = Mbh - steps[step].log_bh_mass[i];
    double prob_Mbh = 1 / (sqrt(2*M_PI) * scatter) * exp(-dMbh*dMbh / (2*scatter*scatter)); //dimension: dlogMbh^{-1}
    
    // Calculate the fractional Eddington ratio (relative to the typical Eddington ratio)
    double eta_frac = eta - steps[step].bh_eta[i];
    double prob_eta;

    // Get prob_eta by interpolating the Eddington ratio distributions.
    if (eta_frac < steps[step].ledd_min[i] || eta_frac > steps[step].ledd_max[i]) 
      prob_eta = 0; 
    else
    {
      double bher_f = (eta_frac-steps[step].ledd_min[i])*steps[step].ledd_bpdex[i];
      int64_t bher_b = bher_f;
      bher_f -= bher_b;
      double p1 = steps[step].bher_dist[i*BHER_BINS+bher_b];
      double p2 = steps[step].bher_dist[i*BHER_BINS+bher_b+1];
      if (bher_b >= BHER_BINS-1) p2 = p1;
      prob_eta = p1 + bher_f*(p2-p1);
    }

    // The total probability is the product of prob_eta, prob_Mbh, and the number density
    // of halos, steps[step].t[i].
    probs[i] = prob_eta * prob_Mbh * steps[step].t[i];
    probs_eta[i] = prob_eta;
    probs_Mbh[i] = prob_Mbh;
    total += probs[i]; //Count the normalization.
  }

  //normalize and convert it to dex^-1.
  if (total > 0)
  {
    for (i = 0; i < M_BINS; i++) probs[i] /= (total*INV_BPDEX); 
  }
  
  // Print out.
  for (i = 0; i < M_BINS; i++) 
  {
    printf("%.1f %.6e %.6e %.6e %.6e\n", (M_MIN + (i + 0.5) * INV_BPDEX), probs[i], probs_Mbh[i], probs_eta[i], steps[step].t[i]);
  }
  


  
  
  return 0;
}
