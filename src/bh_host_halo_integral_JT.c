// Calculate the mass probability distribution of halos hosting a quasar 
// with given SMBH mass ***range*** and Eddington ratio ***range***
// at a given redshift. 
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

// Swap numbers stored in a and b.
void swap(double* a, double *b)
{
  double tmp = *a;
  *a = *b;
  *b = tmp;
}

int main(int argc, char **argv)
{
  int64_t i, j;
  struct smf_fit smf;
  if (argc < 8bh_host_halo_inte) 
  {
    fprintf(stderr, "Usage: %s mass_cache parameter_file z Mbh_low(in Msun) Mbh_high log_eta_low log_eta_high\n", argv[0]);
    exit(1);
  }

  // Read in model parameters, redshift, BH mass, and lower and upper limits of Eddington ratio.
  FILE *param_input = check_fopen(argv[2], "r");
  char buffer[2048];
  fgets(buffer, 2048, param_input);
  read_params(buffer, smf.params, NUM_PARAMS);

  double z = atof(argv[3]);
  double Mbh_low = atof(argv[4]);
  double Mbh_high = atof(argv[5]);
  double eta_low = atof(argv[6]);
  double eta_high = atof(argv[7]);

  if (Mbh_low >= Mbh_high) swap(&Mbh_low, &Mbh_high);
  if (eta_low >= eta_high) swap(&eta_low, &eta_high);
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
  printf("#z=%.2f, Mbh_low=%.6f, Mbh_high=%.6f, eta_low=%.6f, eta_high=%.6f\n", z, Mbh_low, Mbh_high, eta_low, eta_high);
  // printf("#Note: prob = Const.*prob_Mbh*prob_eta*nd_halo, where the constant is chosen so that the summation of prob equals unity.\n");
  printf("#Mh prob_Mbh prob_eta nd_halo [Mpc^(-3) dex^(-1)]\n");
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

  // We have to count the SMBHs in the mass range (Mbh_low, Mbh_high). So a grid
  // must be defined.
  const int64_t mbh_bins = 10000;
  double f_mass[mbh_bins]; // The SMBH mass dependence of AGN duty cycle.
                           // This is needed when we consider a ***range***
                           // of SMBH mass.
  double Mbh_step = (Mbh_high - Mbh_low) / mbh_bins;

  // Calculate f_mass for each SMBH mass.
  for (i = 0; i < mbh_bins; i++)
  {
    double Mbh = Mbh_low + (i + 0.5) * Mbh_step;
    double fac_exp = exp((Mbh - steps[step].smhm.dc_mbh) / steps[step].smhm.dc_mbh_w);
    f_mass[i] = fac_exp / (1 + fac_exp);
  }

  // calculate the normalization of the Gaussian distribution.
  double fac_to_x = 1 / (sqrt(2) * scatter);
  double norm_gauss = fac_to_x / sqrt(M_PI);

  // The redshift-dependent component of AGN duty cycle.
  printf("#Duty cycle factor: %.6f\n", steps[step].smhm.bh_duty);

  for (i = 0; i < M_BINS; i++)
  {
    // Count the number of active SMBHs within the given mass range.
    double prob_Mbh = 0;
    for (j = 0; j < mbh_bins; j++)
    {
      double Mbh = Mbh_low + (j + 0.5) * Mbh_step;
      double dx = (Mbh - steps[step].log_bh_mass[i]) * fac_to_x;
      prob_Mbh += norm_gauss * exp(-dx*dx) * f_mass[i] * Mbh_step * steps[step].smhm.bh_duty;
    }

    // Count the number of SMBHs with Eddington ratios within the Eddington ratio range.
    // This is done by interpolating the Eddington ratio distribution.
    double eta_frac_low = eta_low - steps[step].bh_eta[i];
    double eta_frac_high = eta_high - steps[step].bh_eta[i];
    double prob_eta;

    if (eta_frac_high < steps[step].ledd_min[i] || eta_frac_low > steps[step].ledd_max[i]) 
      prob_eta = 0; 

    else
    {
      double f1 = (eta_frac_low-steps[step].ledd_min[i])*steps[step].ledd_bpdex[i];
      int64_t b1;
      if (f1 < 0)
      {
        b1 = 0;
        f1 = 0;
      }
      else
      {
        b1 = f1;
        f1 -= b1;
      }
      
      double f2 = (eta_frac_high-steps[step].ledd_min[i])*steps[step].ledd_bpdex[i];
      int64_t b2;
      if (f2 > BHER_BINS)
      {
        b2 = BHER_BINS-1;
        f2 = 1.0;
      }
      else
      {
        b2 = f2;
        f2 -= b2;
      }

      for (j=b1+1; j<b2; j++)
      {
        prob_eta += steps[step].bher_dist[i*BHER_BINS+j];
      }
      prob_eta += (1.0 - f1) * steps[step].bher_dist[i*BHER_BINS+b1];
      prob_eta += f2 * steps[step].bher_dist[i*BHER_BINS+b2];
      prob_eta /= steps[step].ledd_bpdex[i]; 
    }
    // The total probability is a product of prob_eta, prob_Mbh, and the
    // halo number density, steps[step].t[i].
    probs[i] = prob_eta * prob_Mbh * steps[step].t[i] * BPDEX;
    probs_eta[i] = prob_eta;
    probs_Mbh[i] = prob_Mbh;
    total += probs[i];
  }

  total *= INV_BPDEX;

  printf("#Total Number Density: %.6e Mpc^(-3)\n", total);
  // Print them out.
  for (i = 0; i < M_BINS; i++) 
  {
    printf("%.1f %.6e %.6e %.6e %.6e\n", (M_MIN + (i + 0.5) * INV_BPDEX), probs[i], probs_Mbh[i], probs_eta[i], steps[step].t[i] * BPDEX);
  }
  


  
  
  return 0;
}
