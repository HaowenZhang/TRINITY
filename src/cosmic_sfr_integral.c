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
#include "mah.h"
#include "check_syscalls.h"
#include "expcache2.h"

#define NUM_ZS 9
extern int64_t num_outputs;
extern struct timestep *steps;

float total_sm(int64_t i, int64_t j) {
  int64_t k;
  double t = 0;
  for (k=0; k<i; k++) {
    t += steps[i].sm_hist[num_outputs*j+k];
  }
  return t;
}

int main(int argc, char **argv)
{
  struct smf_fit the_smf;
  int64_t i;

  if (argc < 3) 
  {
    fprintf(stderr, "Usage: %s mass_cache parameter_file (> output_file)\n", argv[0]);
    exit(1);
  }

  // Read in model parameters
  FILE *param_input = check_fopen(argv[2], "r");
  char buffer[2048];
  fgets(buffer, 2048, param_input);
  read_params(buffer, smf.params, NUM_PARAMS);

  gen_exp10cache();
  setup_psf(1);
  load_mf_cache(argv[1]);
  init_timesteps();
  INVALID(the_smf) = 0;
  calc_sfh(&the_smf);

  float total_sm = 0;
  float total_sm2 = 0;
  float total_sm3 = 0;
  float last_scale = 0;
  int64_t j = 0;
  for (i=0; i<num_outputs; i++) {
    total_sm = total_sm2 = 0;
    for (j=0; j<=i; j++)
      total_sm += steps[j].cosmic_sfr*steps[i].dt*steps[i].smloss[j];
    last_scale = 0;
    for (j=0; j<=i; j++) {
      float z = 1.0/(0.5*(steps[j].scale+last_scale))-1.0 - 1.243;
      total_sm2 += 0.180/(pow(10, -0.997*z)+pow(10,0.241*z))*steps[j].dt*(1-0.28);
      last_scale = steps[j].scale;
    }
    total_sm3 = 0;
    float total_sm4 = 0;
    for (j=0; j<M_BINS; j++) {
      if (!isfinite(steps[i].sm_avg[j])) continue;
      total_sm3 += steps[i].sm_avg[j]*steps[i].t[j];
      total_sm4 += steps[i].sm_icl[j]*steps[i].t[j];
    }
    total_sm4 += total_sm3;
    printf("%f %e %e %e %e\n", 1.0/steps[i].scale-1.0, total_sm2, total_sm, total_sm3, total_sm4);
  }
  return 0;
}
