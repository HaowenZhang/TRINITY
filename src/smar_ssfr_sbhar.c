// Calculate the average specific halo mass accretion rate, specific
// star formation rate, and specific SMBH accretion rate at a given
// halo mass, as functions of redshift.
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "observations.h"
#include "smf.h"
#include "mah.h"
#include "all_smf.h"
#include "distance.h"
#include "integrate.h"
#include "mlist.h"
#include "calc_sfh.h"
#include "expcache2.h"

// Calculate dark matter halo accretion rates
// given the snapshot number, n, and the halo 
// mass bin number, j.
float _mar_from_mbins(int64_t n, int64_t j) 
{
  int64_t i;
  if (!n) return pow(10, M_MIN+(j+0.5)*INV_BPDEX)/steps[n].dt;
  if (n>=num_outputs-1) return _mar_from_mbins(num_outputs-2, j);
  if (j>=M_BINS-1) return 0;
  //Find out average progenitor mass:
  double sum = 0;
  double count = 0;
  for (i=0; i<M_BINS; i++) {
    if (!steps[n].mmp[i*M_BINS + j]) continue;
    sum += pow(10, M_MIN+(i+0.5)*INV_BPDEX)*steps[n].mmp[i*M_BINS+j];
    count += steps[n].mmp[i*M_BINS+j];
  }
  if (!count) return 0;
  sum /= count;
  return ((pow(10, M_MIN+(j+0.5)*INV_BPDEX) - sum)/steps[n].dt);
}

// Calculate dark matter halo accretion rates
// given the snapshot number, n, and the halo 
// mass bin number, j. The only difference with
// _mar_from_mbins is that we take the average
// value between the two snapshots.
float mar_from_mbins(int64_t n, int64_t j) 
{
  float mar1 = _mar_from_mbins(n,j);
  float mar2 = _mar_from_mbins(n+1,j);
  return (0.5*(mar1+mar2));
}

int main(int argc, char **argv)
{
  int64_t i;
  struct smf_fit smf;
 
  
  if (argc < 4) 
  {
    fprintf(stderr, "Usage: %s mass_cache param_file halo_mass (> output_file)\n", argv[0]);
    exit(1);
  }

  // Read in model parameters
  FILE *param_input = check_fopen(argv[2], "r");
  char buffer[2048];
  fgets(buffer, 2048, param_input);
  read_params(buffer, smf.params, NUM_PARAMS);

  // Find out the halo mass bin to do the interpolation.
  double mh = atof(argv[3]);
  double mf = (mh - M_MIN) * BPDEX - 0.5;
  int mb = mf;
  mf -= mb;

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
  
  
  printf("#z Mh SMAR SM SSFR Mbh SBHAR\n");


  for (i=0; i<num_outputs; i++) 
  {
    double mbh = steps[i].bh_mass_avg[mb] + mf * (steps[i].bh_mass_avg[mb+1] - steps[i].bh_mass_avg[mb]);
    double mstar = steps[i].sm_avg[mb] + mf * (steps[i].sm_avg[mb+1] - steps[i].sm_avg[mb]);
    double bhar = steps[i].bh_acc_rate[mb] + mf * (steps[i].bh_acc_rate[mb+1] - steps[i].bh_acc_rate[mb]);
    double sfr = steps[i].sfr[mb] + mf * (steps[i].sfr[mb+1] - steps[i].sfr[mb]);
    printf("%f %f %e %e %e %e %e\n", 1 / steps[i].scale - 1.0, mh, ma_rate_avg_mnow(mh, steps[i].scale) / exp10(mh), 
                  mstar, sfr/mstar,
                  mbh, bhar/mbh);
  }
  return 0;
}
