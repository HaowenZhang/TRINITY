// Calculate the galaxy stellar mass functions between redshifts z_low and z_high,
// as a function of ***observed*** stellar mass.
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

#define MASS_START 7
#define MASS_STOP 12.5
#define MASS_BPDEX 10
#define MASS_BINS (int)((MASS_STOP-MASS_START)*(MASS_BPDEX))
#define MASS_STEP (1.0/(double)MASS_BPDEX)

float integrate_smf(float z_low, float z_high, double m, struct smf_fit *fit) 
{
  float smf_val, epsilon;
  double v_high = comoving_volume(z_high);
  double v_low = comoving_volume(z_low);
  double weight = fabs(v_high - v_low);

  if (z_low != z_high) 
  {
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
  float z_low, z_high, m;
  struct smf_fit smfs[4];
  float smf_points[4][MASS_BINS];
  int i,j;
  if (argc<4+NUM_PARAMS) 
  {
    fprintf(stderr, "Usage: %s z_low z_high mass_cache (mcmc output)\n", argv[0]);
    exit(1);
  }
  // Read in the lower and higher redshifts and model parameters.
  z_low = atof(argv[1]);
  z_high = atof(argv[2]);
  for (i=0; i<NUM_PARAMS; i++)
    smfs[0].params[i] = atof(argv[i+4]);
  smfs[0].params[NUM_PARAMS] = 0;

  smfs[1] = smfs[2] = smfs[3] = smfs[0];

  // Turn off the built-in GSL error handler that kills the program
  // when an error occurs. We handle the errors manually.
  gsl_set_error_handler_off();
  // We use non-linear scaling relation between the radiative and total Eddington ratios.
  nonlinear_luminosity = 1;
  // Set up the PSF for stellar mass functions. See observations.c.
  setup_psf(1);
  // Load cached halo mass functions.
  load_mf_cache(argv[3]);
  // Initialize all the timesteps/snapshots.
  init_timesteps();
  // Fix some model parameters.
  assert_model(smfs + 3);  
  // Calculate the star-formation histories and black hole histories. See calc_sfh.c.
  calc_sfh(smfs + 3);
  // Calculate stellar mass functions.
  for (j=0; j<4; j++) 
  {
    for (i=0; i<MASS_BINS; i++) 
    {
      m = MASS_START + i*MASS_STEP;
      smf_points[j][i] = integrate_smf(z_low, z_high, m, smfs+j);
    }
    if (j==1) setup_psf(1);
  }
  // output
  for (i=0; i<MASS_BINS; i++) 
  {
    m = MASS_START + i*MASS_STEP;
    printf("%f %f %f %f %f\n", m , smf_points[0][i], smf_points[1][i], smf_points[2][i], smf_points[3][i]); 
  }

  return 0;
}
