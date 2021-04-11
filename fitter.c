#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "all_smf.h"
#include "smf.h"
#include "observations.h"
#include "mlist.h"

#define FITTING_ROUNDS 4
#define CHI2_LIMIT 1e-2
#define STEP_LIMIT 8
#define MIN_FIT_LIMIT 10
//#define STEP_RATIO (1/2.0) //How quickly step sizes shrink
//#define STEP_RATIO (7/10.0)
#define STEP_RATIO 1.1

extern int no_z_scaling;
extern int no_matching_scatter;
extern int no_systematics;
extern int no_obs_scatter;

extern double frac_below8[MBH_BINS];
extern double frac_below11[MBH_BINS];

int iterations = 0;

extern int ind_var[4];


float fitter(float *params) {
  struct smf_fit test;
  int i;
  for (i=0; i<NUM_PARAMS; i++)
    test.params[i] = params[i];

  assert_model(&test);
  // chi2_type(test);

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

int main(int argc, char **argv) {
  int i;
  float params[NUM_PARAMS], steps[NUM_PARAMS];
  init_frac_below8();
  gsl_set_error_handler_off();
  init_mcmc_from_args(argc, argv);
  read_params_and_steps(params, steps);
  struct smf_fit test;
  for (i=0; i<NUM_PARAMS; i++)
    test.params[i] = params[i];

  //chi2_type(test);
  fprintf(stderr, "Initial Chi2 before asserting: %e\n", calc_chi2(params));
  assert_model(&test);
  chi2_type(test);
  // for (i=0; i<20; i++) 
  fprintf(stderr, "Initial Chi2 after asserting: %e\n", calc_chi2(params));

  // assert_model(&test);
  // chi2_type(test);

  //for (i=0; i<20; i++) chi2_type(test);
  //return 0;
  for (i=0; i<4; i++) 
  {
	  fprintf(stderr, "Starting of the %d-th round.\n", i);
	  fitting_round(params, steps);
  }
  for (i=0; i<NUM_PARAMS-1; i++) printf("%.12f ", params[i]);
  printf("%.12f\n", params[NUM_PARAMS-1]);
  printf("Chi2: %e\n", calc_chi2(params));
  for (i=0; i<NUM_PARAMS; i++)
    test.params[i] = params[i];
  chi2_type(test);
  printf("Iterations: %d\n", iterations);
  return 0;
}

//inline 
float minimum_delta(float chi2_l, float chi2, float chi2_r)
{
  float d_chi2_dx = (chi2_r - chi2_l) / 2.0;
  float d2_chi2_dx2 = (chi2_r + chi2_l - 2*chi2);
  if (chi2 > chi2_r && chi2 > chi2_l) // We are close to a local maximum
    return ((chi2_l < chi2_r) ? -1 : 1);
  return (d2_chi2_dx2 ? (-d_chi2_dx / d2_chi2_dx2) : 0);
}


float improve_fit(float *params, float *steps, int i, float chi2)
{
  float chi2_l, chi2_r, p, dx, chi2_new;
  p = params[i];
  params[i] = p + steps[i];
  chi2_r = calc_chi2(params);
  params[i] = p - steps[i];
  chi2_l = calc_chi2(params);

  dx = minimum_delta(chi2_l, chi2, chi2_r);
  if (fabs(dx) > STEP_LIMIT) dx = copysign(STEP_LIMIT, dx);
  params[i] = p + dx*steps[i];
  chi2_new = calc_chi2(params);

  if (chi2_l < chi2_new) { chi2_new = chi2_l; params[i] = p-steps[i]; }
  if (chi2_r < chi2_new) { chi2_new = chi2_r; params[i] = p+steps[i]; }
  if (chi2 < chi2_new) { chi2_new = chi2; params[i] = p; }
  //fprintf(stderr, "chi2_l: %e, chi2_r: %e, chi2_new: %e, \n", chi2_l, chi2_r, chi2_new);
  steps[i]*=STEP_RATIO;
  //  fprintf(stderr, ".");
  return chi2_new;
}


void fitting_round(float *params, float *steps)
{
  float last_chi2, chi2;
  float cur_steps[NUM_PARAMS];
  int i, j=0;
  
  memcpy(cur_steps, steps, sizeof(float)*NUM_PARAMS);
  chi2 = calc_chi2(params);
  //fprintf(stderr, "chi2 before the round: %.3f\n", chi2);
  last_chi2 = chi2*(2+CHI2_LIMIT);
  while (chi2*(1+CHI2_LIMIT) < last_chi2 || (j++ < MIN_FIT_LIMIT)) {
    last_chi2 = chi2;
    for (i=0; i<NUM_PARAMS; i++)
 
 //   for (i=0; i<4; i++)
{
      //if (i >= 63) continue;
      chi2 = improve_fit(params, cur_steps, i, chi2);
      fprintf(stderr, "chi2: %.3f\n", chi2);
}
  }
}


float calc_chi2(float *params) {
  return fitter(params);
}

void read_params(char *buffer, float *data, int max_n) {
  int num_entries = 0;
  char *cur_pos = buffer, *end_pos;
  float val = strtod(cur_pos, &end_pos);
  while (cur_pos != end_pos && num_entries < max_n) {
    data[num_entries] = val;
    num_entries++;
    cur_pos=end_pos;
    while (*cur_pos==' ' || *cur_pos=='\t' || *cur_pos=='\n') cur_pos++;
    val = strtod(cur_pos, &end_pos);
  }
}

void read_params_and_steps(float *params, float *steps) {
  int i;
  char buffer[2048];
  float inputs[NUM_PARAMS];
  fgets(buffer, 2048, stdin);
  read_params(buffer, params, NUM_PARAMS);
  fgets(buffer, 2048, stdin);
  read_params(buffer, steps, NUM_PARAMS);
  for (i=0; i<NUM_PARAMS; i++) {
    //steps[i]-=params[i];
    //steps[i] = fabs(steps[i]);
    steps[i] = 0.03 * fabs(steps[i]);
  }
  fgets(buffer, 2048, stdin);
  read_params(buffer, inputs, NUM_PARAMS);
  //for (i=0; i<NUM_PARAMS; i++) {
  //  inputs[i] = fabs(params[i] - inputs[i]);
  //  if (inputs[i] > steps[i]) steps[i] = inputs[i];
  //}
}
