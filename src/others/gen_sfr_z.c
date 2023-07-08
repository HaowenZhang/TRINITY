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

extern struct timestep *steps;
extern int64_t num_outputs;
#define NUM_ZS 4

int main(int argc, char **argv)
{
  struct smf_fit the_smf;
  int i,j, k;
  float z[NUM_ZS] = {4,2,1,0.1};
  
  if (argc<2+NUM_PARAMS) {
    fprintf(stderr, "Usage: %s mass_cache (mcmc output)\n", argv[0]);
    exit(1);
  }
  for (i=0; i<NUM_PARAMS; i++)
    the_smf.params[i] = atof(argv[i+2]);

  setup_psf(1);
  load_mf_cache(argv[1]);
  init_timesteps();

  calc_sfh(&the_smf);
  j=0;
  for (i=0; i<num_outputs; i++) {
    if (steps[i].scale > 1.0/(1+z[j])) {
      printf("#z = %f\n", z[j]);
      for (k=0; k<M_BINS; k++) {
	if (steps[i].sfr[k]>0)
	  printf("%f %e\n", M_MIN+k*INV_BPDEX, steps[i].sfr[k]);
      }
      printf("\n");
      j++;
      if (j==NUM_ZS) break;
    }
  }
  return 0;
}
