// Calculate the host galaxy stellar matter halo mass function for 
// a quasar population at a given redshift, with bolometric luminosities
// within the range (lbol_low, lbol_high) (in log10 units).
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

// Stellar mass bins over which the mass function is calculated.
#define MASS_START 7
#define MASS_STOP 12.5
#define MASS_BPDEX 10
#define MASS_BINS (int)((MASS_STOP-MASS_START)*(MASS_BPDEX))
#define MASS_STEP (1.0/(double)MASS_BPDEX)

// Evaluate the total galaxy stellar mass function at a given 
// mass m, between redshift (z_low, z_high).
float integrate_smf(float z_low, float z_high, double m, struct smf_fit *fit) 
{
  float smf_val, epsilon;
  double v_high = comoving_volume(z_high);
  double v_low = comoving_volume(z_low);
  double weight = fabs(v_high - v_low);

  if (z_low != z_high) 
  {
    // For more details of chi2_err_helper(), see observations.c.
    epsilon = chi2_err_helper((v_high+v_low)/2.0, &m)*weight*1e-5;
    smf_val = adaptiveSimpsons(chi2_err_helper, &m,
       v_low, v_high, epsilon, 10);
    smf_val /= weight;
  }
  else 
  {
    smf_val = chi2_err_helper(v_low, &m);
  }
  return (smf_val ? log10f(smf_val) : -15);
}

int main(int argc, char **argv)
{
  int64_t i, j, k, l;
  struct smf_fit smf;
  
  if (argc < 6) 
  {
    fprintf(stderr, "Usage: %s mass_cache param_file z lbol_low(in erg/s) lbol_high(in erg/s) (> output_file)\n", argv[0]);
    exit(1);
  }

  // Read in model parameters
  FILE *param_input = check_fopen(argv[2], "r");
  char buffer[2048];
  fgets(buffer, 2048, param_input);
  read_params(buffer, smf.params, NUM_PARAMS);

  double z = atof(argv[3]);
  double Lbol_low = atof(argv[4]);
  double Lbol_high = atof(argv[5]);

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
  // assert some model parameters to pre-defined values. See all_smf.c.
  assert_model(&smf);
  // Calculate the star-formation histories and black hole histories. See calc_sfh.c.
  calc_sfh(&smf);
   

  printf("#Is the model invalid? %e\n", INVALID(smf));
  printf("#z=%.2f, Lbol_low=%.6f, Lbol_high=%.6f\n", z, Lbol_low, Lbol_high);
  printf("#Mbh Phi(Mbh)\n");
  double t,m;
  int64_t step;
  double f;
  // Calculate the # of the snapshot that is the closest to
  // the input redshift.
  calc_step_at_z(z, &step, &f);
  //The probabilities of having luminosities between (Lbol_low, Lbol_high), as a function of BH mass.
  double prob_lbol[MBH_BINS] = {0}; 


  // double sm_scatter = steps[step].smhm.scatter * steps[step].smhm.bh_gamma;
  // double scatter = sqrt(sm_scatter*sm_scatter 
  //                     + steps[step].smhm.bh_scatter*steps[step].smhm.bh_scatter);

  // Define the SMBH mass range over which we calculate prob_lbol.
  double mbh_min = steps[step].bh_mass_min, mbh_max = steps[step].bh_mass_max;
  double mbh_inv_bpdex = (mbh_max - mbh_min) / MBH_BINS;
  double mu = steps[step].smhm.mu;
  // Calculate prob_lbol.
  for (i=0; i<MBH_BINS; i++)
  {
    double mbh = mbh_min + (i + 0.5) * mbh_inv_bpdex;
    // luminosity bin numbers corresponding to lbol_low and lbol_high.
    double lbol_f_low = (Lbol_low - LBOL_MIN) * LBOL_BPDEX;
    double lbol_f_high = (Lbol_high - LBOL_MIN) * LBOL_BPDEX;
    int64_t lbol_b_low = lbol_f_low; lbol_f_low -= lbol_b_low;
    int64_t lbol_b_high = lbol_f_high; lbol_f_high -= lbol_b_high;
    if (lbol_b_low < 0) {lbol_b_low = 0; lbol_f_low = 1;}
    // If even the lbol_high does not lie in the meaningful luminosity range,
    // just leave prob_lbol[j] = 0;
    if (lbol_b_high < 0) continue;
    if (lbol_b_low >= LBOL_BINS - 1) {lbol_b_low = LBOL_BINS - 2; lbol_f_low = 1;}
    if (lbol_b_high >= LBOL_BINS - 1) {lbol_b_high = LBOL_BINS - 2; lbol_f_high = 1;}
    // Add up the probabilities between lbol_low and lbol_high.
    double prob_lum = (1 - lbol_f_low) * steps[step].lum_dist_full[i*LBOL_BINS+lbol_b_low];
    double prob_tot = 0;
    for (j=0; j<LBOL_BINS; j++) prob_tot += steps[step].lum_dist_full[i*LBOL_BINS+j];
    for (j=lbol_b_low+1; j<lbol_b_high; j++) prob_lum += steps[step].lum_dist_full[i*LBOL_BINS+j];
    prob_lum += lbol_f_high * steps[step].lum_dist_full[i*LBOL_BINS+lbol_b_high];
    
    // Duty cycle, i.e., the fraction of SMBHs that are active.
    double dc = steps[step].smhm.bh_duty;
    double f_mass = exp((mbh - steps[step].smhm.dc_mbh) / steps[step].smhm.dc_mbh_w);
    f_mass = f_mass / (1 + f_mass);
    dc *= f_mass;
    if (dc < 1e-4) dc = 1e-4;

    // Normalize the probability, and scale it by the duty cycle.
    prob_lum /= prob_tot;
    prob_lbol[i] = prob_lum * dc;
    fprintf(stderr, "%f %e\n", mbh, prob_lbol[i]);
  }

  // Iterate over each ***observed*** stellar mass bin.
  for (i=0; i<MASS_BINS; i++) 
  {
    //observed stellar mass.
    m = MASS_START + i*MASS_STEP; 
    // Calculate the observed stellar mass function.
    double smf_tot = integrate_smf(z, z, m, &smf);
    // Calculate the median SMBH mass given the stellar mass.
    double mb = bulge_mass(m + mu, 1.0 / (1 + z));
    double mbh_med = calc_bh_at_bm(mb, steps[step].smhm);
    // The probability of this stellar mass bin's hosting a quasar between (lbol_low, lbol_high).
    double p_l = 0;

    // Add up the contribution from all SMBH mass bins to this halo mass bin.
    for (j=0; j<MBH_BINS; j++)
    {
      double mbh = mbh_min + (j + 0.5) * mbh_inv_bpdex;
      double dmbh = (mbh - mbh_med) / steps[step].smhm.bh_scatter;
      // The probability of having a BH mass.
      double prob_tmp = 1 / sqrt(2*M_PI) / steps[step].smhm.bh_scatter * exp(-0.5*dmbh*dmbh) * mbh_inv_bpdex;
      p_l += prob_tmp * prob_lbol[j];

    }
    
    // Multiply by the number density term.
    printf("%f %e\n", m, exp10(smf_tot) * p_l);

  }
  
  return 0;
}
