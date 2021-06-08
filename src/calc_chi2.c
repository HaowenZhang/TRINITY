// Calculate the chi2 given the parameter set.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "all_smf.h"
#include "smf.h"
#include "observations.h"

#define FITTING_ROUNDS 4
#define CHI2_LIMIT 1e-2
#define STEP_LIMIT 8
#define MIN_FIT_LIMIT 10
#define STEP_RATIO (1/2.0) //How quickly step sizes shrink
extern int no_z_scaling;
extern int no_matching_scatter;
extern int no_systematics;
extern int no_obs_scatter;
extern int64_t num_obs_smfs;
extern struct obs_smf_point *obs_smfs;
int iterations = 0;

float fitter(float *params) 
{
  struct smf_fit test;
  int i;
  for (i=0; i<NUM_PARAMS; i++)
    test.params[i] = params[i];

  assert_model(&test);

  for (i=0; i<NUM_PARAMS; i++)
    params[i] = test.params[i];
  iterations++;
  float err = all_smf_chi2_err(test);
  if (!isfinite(err) || err<0) return 1e30;
  return err;
}

void fitting_round(float *params, float *steps);
void read_params_and_steps(float *params, float *steps);
float calc_chi2(float *params);

int main(int argc, char **argv) 
{
  int i;
  float params[NUM_PARAMS], steps[NUM_PARAMS];

  init_mcmc_from_args(argc, argv);
  char buffer[1024];
  fgets(buffer, 1024, stdin);
  read_params(buffer, params, NUM_PARAMS);
  double chi2_total, prior, chi2_data;
  chi2_total = calc_chi2(params);
  double chi2_type[7] = {0};
  for (i = 0; i < 7; i++)
  {
    for (int j = 0; j < num_obs_smfs; j++)
    {
      if (obs_smfs[j].type == i) chi2_type[i] += obs_smfs[j].chi2;
    }
    chi2_data += chi2_type[i];
    printf("chi2 type %d: %e\n", i, chi2_type[i]);
  }
  prior = chi2_total - chi2_data;
  return 0;
}


float calc_chi2(float *params) 
{
  return fitter(params);
}

// void read_params(char *buffer, float *data, int max_n) 
// {
//   int num_entries = 0;
//   char *cur_pos = buffer, *end_pos;
//   float val = strtod(cur_pos, &end_pos);
//   while (cur_pos != end_pos && num_entries < max_n) 
//   {
//     data[num_entries] = val;
//     num_entries++;
//     cur_pos=end_pos;
//     while (*cur_pos==' ' || *cur_pos=='\t' || *cur_pos=='\n') cur_pos++;
//     val = strtod(cur_pos, &end_pos);
//   }
// }

void read_params_and_steps(float *params, float *steps) 
{
  int i;
  char buffer[1024];
  float inputs[NUM_PARAMS];
  fgets(buffer, 1024, stdin);
  read_params(buffer, params, NUM_PARAMS);
}
