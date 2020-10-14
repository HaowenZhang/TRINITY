#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include "all_smf.h"
#include "smf.h"
#include "observations.h"

extern int no_matching_scatter;
extern int no_systematics;
extern int no_obs_scatter;

float fitter(float *params) {
  struct smf_fit test;
  int i;
  for (i=0; i<NUM_PARAMS; i++)
    test.params[i] = params[i];
  
  if (no_matching_scatter)
    SCATTER(test) = 0;

  if (no_systematics) {
    KAPPA(test) = 0;
    MU(test) = 0;
    KAPPA_A(test) = 0;
    MU_A(test) = 0;
    BURST_DUST_AMP(test) = 0;
    BURST_DUST_Z(test) = 1;
    BURST_DUST(test) = 0;
    SCATTER_A(test) = 0;
  }

  if (no_obs_scatter)
    SIGMA_Z(test) = 0;

  INVALID(test) = 0;

  for (i=0; i<NUM_PARAMS; i++)
    params[i] = test.params[i];
  float err = all_smf_chi2_err(test);
  if (!isfinite(err) || err<0) return 1e30;
  return err;
}

void read_params_from_stdin(float *params);
float calc_chi2(float *params);

int main(int argc, char **argv) {
  float params[NUM_PARAMS];

  init_mcmc_from_args(argc, argv);
  while (1) {
    read_params_from_stdin(params);
    double chi2 = calc_chi2(params);
    fprintf(stdout, "%f\n", chi2);
    fflush(stdout);
  }
  return 0;
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

void read_params_from_stdin(float *params) {
  char buffer[1024];
  fgets(buffer, 1024, stdin);
  if (!strncmp(buffer, "quit", 4)) exit(0);
  read_params(buffer, params, NUM_PARAMS);
}
