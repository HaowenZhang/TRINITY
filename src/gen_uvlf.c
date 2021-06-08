// Calculate the galaxy UV luminosity functions between redshifts z_low and z_high,
// as a function of ***observed*** UV magnitude.
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

#define UV_START -25 
#define UV_STOP -15
#define MAG_BPMAG 10
#define MAG_BINS (int)((UV_STOP-UV_START)*(MAG_BPMAG))
#define MAG_STEP (1.0/(double)MAG_BPMAG)

float integrate_uvlf(float z_low, float z_high, double uv, struct smf_fit *fit) 
{
  float uvlf_val, epsilon;
  double v_high = comoving_volume(z_high);
  double v_low = comoving_volume(z_low);
  double weight = fabs(v_high - v_low);

  if (z_low != z_high) 
  {
      uvlf_val = chi2_err_helper_uv((v_high+v_low)/2.0, &uv);
  }
  else 
  {
    uvlf_val = chi2_err_helper_uv(v_low, &uv);
  }
  return (uvlf_val ? log10f(uvlf_val) : -15);
}

int main(int argc, char **argv)
{
  float z_low, z_high, uv;
  struct smf_fit smfs[4];
  float uvlf_points[4][MAG_BINS];
  int i,j;
  
  if (argc < 5) 
  {
    fprintf(stderr, "Usage: %s z_low z_high mass_cache param_file (> output_file)\n", argv[0]);
    exit(1);
  }
  // Read in the lower and higher redshifts and model parameters.
  z_low = atof(argv[1]);
  z_high = atof(argv[2]);

  // Read in model parameters
  FILE *param_input = check_fopen(argv[4], "r");
  char buffer[2048];
  fgets(buffer, 2048, param_input);
  read_params(buffer, smfs[0].params, NUM_PARAMS);
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
  // Calculate and output observed UV luminosity functions.
  for (j=0; j<4; j++) 
  {

    for (i=0; i<MAG_BINS; i++) 
    {
      uv = UV_START + i*MAG_STEP;
      uvlf_points[j][i] = integrate_uvlf(z_low, z_high, uv, smfs+j);
    }
    if (j==1) setup_psf(1);
  }

  for (i=0; i<MAG_BINS; i++) 
  {
    uv = UV_START + i*MAG_STEP;
    printf("%f %f %f %f %f\n", uv , uvlf_points[0][i], uvlf_points[1][i], uvlf_points[2][i], uvlf_points[3][i]); 
  }

  return 0;
}
