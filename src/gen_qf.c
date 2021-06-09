// Calculate the galaxy quenched fractions between redshifts z_low and z_high,
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

float integrate_qf(float z_low, float z_high, double m, struct smf_fit *fit) 
{
  float smf_val, epsilon;
  double v_high = comoving_volume(z_high);
  double v_low = comoving_volume(z_low);
  double weight = fabs(v_high - v_low);

  if (z_low != z_high) 
  {
      epsilon = chi2_err_helper_qf((v_high+v_low)/2.0, &m)*weight*1e-5;
      //if (PHI_HIGH_Z < z_high) epsilon *= 1e1;
      smf_val = adaptiveSimpsons(chi2_err_helper_qf, &m,
				 v_low, v_high, epsilon, 10);
      smf_val /= weight;
  }
  else 
  {
    smf_val = chi2_err_helper_qf(v_low, &m);
  }
  return smf_val;
}

int main(int argc, char **argv)
{
  float z_low, z_high, m;
  struct smf_fit smfs[4];
  float qf_points[4][MASS_BINS];
  int i,j;

  if (argc < 5) 
  {
    fprintf(stderr, "Usage: %s z_low z_high mass_cache param_file (> output_file)\n", argv[0]);
    exit(1);
  }
  z_low = atof(argv[1]);
  z_high = atof(argv[2]);

  // Read in model parameters
  FILE *param_input = check_fopen(argv[4], "r");
  char buffer[2048];
  fgets(buffer, 2048, param_input);
  read_params(buffer, smfs[0].params, NUM_PARAMS);

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

  smfs[1] = smfs[2] = smfs[3] = smfs[0];
  // Fix some model parameters.
  assert_model(smfs + 3);  
  // Calculate star formation histories.
  calc_sfh(smfs + 3);

  // Calculate and output the quenched fractions.
  for (j=0; j<4; j++) 
  {
    for (i=0; i<MASS_BINS; i++) 
    {
      m = MASS_START + i*MASS_STEP;
      qf_points[j][i] = integrate_qf(z_low, z_high, m, smfs+j);
    }
    if (j==1) setup_psf(1);
  }

  for (i=0; i<MASS_BINS; i++) 
  {
    m = MASS_START + i*MASS_STEP;
    printf("%f %f %f %f %f\n", m , qf_points[0][i], qf_points[1][i], qf_points[2][i], qf_points[3][i]); 
  }

  return 0;
}
