// Calculate a series galaxy stellar mass functions between redshifts z_low and z_high,
// as a function of ***observed*** stellar mass. The stellar mass functions are generated
// based on a series of steps in an MCMC chain.
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
#include "check_syscalls.h"

#define MASS_START 7
#define MASS_STOP 12.5
#define MASS_BPDEX 10
#define MASS_BINS (int)((MASS_STOP-MASS_START)*(MASS_BPDEX))
#define MASS_STEP (1.0/(double)MASS_BPDEX)
#define ONE_SIGMA 0.682689492137
// The indices corresponding to +/- 1-sigma range.
#define SIGMA_UP ((int)(num_smfs*(0.5+ONE_SIGMA/2.0)))
#define SIGMA_DOWN ((int)(num_smfs*(0.5-ONE_SIGMA/2.0)))


void read_params(char *buffer, double *params, int max_n);
int float_compare(const void *a, const void *b);

float integrate_smf(float z_low, float z_high, double m, struct smf_fit *fit) 
{
  float smf_val, epsilon;
  double v_high = comoving_volume(z_high);
  double v_low = comoving_volume(z_low);
  double weight = fabs(v_high - v_low);

  if (z_low != z_high) 
  {
      epsilon = chi2_err_helper((v_high+v_low)/2.0, &m)*weight*1e-5;
      //if (PHI_HIGH_Z < z_high) epsilon *= 1e1;
      smf_val = adaptiveSimpsons(chi2_err_helper, &m,
				 v_low, v_high, epsilon, 10);
      smf_val /= weight;
  }
  else {
    smf_val = chi2_err_helper(v_low, &m);
  }
  return (smf_val ? log10f(smf_val) : -15);
}

int main(int argc, char **argv)
{
  float z_low, z_high, m;
  struct smf_fit smf;
  float *smf_points[MASS_BINS];
  int64_t num_smfs = 0;
  int i,j;
  FILE *input;
  char buffer[1024];

  if (argc<6) 
  {
    fprintf(stderr, "Usage: %s z_low z_high mass_cache num_smfs mcmc_file\n", argv[0]);
    exit(1);
  }

  // Read in the lower and higher redshifts and model parameters.
  z_low = atof(argv[1]);
  z_high = atof(argv[2]);
  num_smfs = atof(argv[4]);

  // Open the file containing the MCMC chain, and allocate spaces
  // to store the stellar mass functions.
  input = check_fopen(argv[5], "r");
  for (i=0; i<MASS_BINS; i++)
    smf_points[i] = check_realloc(NULL, sizeof(float)*num_smfs, "ASMF");

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
  // Read in MCMC steps to calculate stellar mass functions. Note
  // that the first MCMC step must contain the best-fitting parameters.
  for (i=0; i<num_smfs; i++) 
  {
    if (!fgets(buffer, 1024, input)) break;
    read_params(buffer, smf.params, NUM_PARAMS);
    INVALID(smf) = 0;
    // For each MCMC step, calculate the star formation histories.
    calc_sfh(&smf);
    // Then calculate stellar mass functions.
    for (j=0; j<MASS_BINS; j++) 
    {
      m = MASS_START + j*MASS_STEP;
      smf_points[j][i] = integrate_smf(z_low, z_high, m, &smf);
    }
  }
  num_smfs = i;

  // Sort the stellar mass function values at each stellar mass value
  // to calculate the confidence intervals.
  for (i=0; i<MASS_BINS; i++) 
  {
    m = MASS_START + i*MASS_STEP;
    float best = smf_points[i][0];
    qsort(smf_points[i], num_smfs, sizeof(float), float_compare);    
    printf("%f %f %f %f\n", m, best, smf_points[i][SIGMA_UP]-best, best-smf_points[i][SIGMA_DOWN]);
  }
  return 0;
}

// Helper function to do sorting.
int float_compare(const void *a, const void *b) 
{
  float c = *((float *)a);
  float d = *((float *)b);
  if (c<d) return -1;
  if (c>d) return 1;
  return 0;
}

void read_params(char *buffer, double *params, int max_n) {
  int num_entries = 0;
  char *cur_pos = buffer, *end_pos;
  float val = strtod(cur_pos, &end_pos);
  while (cur_pos != end_pos && num_entries < max_n) {
    params[num_entries] = val;
    num_entries++;
    cur_pos=end_pos;
    while (*cur_pos==' ' || *cur_pos=='\t' || *cur_pos=='\n') cur_pos++;
    val = strtod(cur_pos, &end_pos);
  }
}
