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

#define NUM_ZS 9
extern int64_t num_outputs;
extern struct timestep *steps;

int main(int argc, char **argv)
{
  struct smf_fit the_smf;
  int i, j;

  if (argc<3+NUM_PARAMS) {
    fprintf(stderr, "Usage: %s mass_cache (mcmc output)\n", argv[0]);
    exit(1);
  }
  for (i=0; i<NUM_PARAMS; i++)
    the_smf.params[i] = atof(argv[i+2]);

  setup_psf(1);
  load_mf_cache(argv[1]);
  init_timesteps();
  calc_sfh(&the_smf);
  
  for (i=1; i<num_outputs-1; i++) {
    printf("%f ", steps[i].scale);
    for (j=11; j<16; j++) {
      int b = (j-M_MIN)*BPDEX;
      printf("%e ", steps[i].merged_frac[b]);
    }
    printf("\n");
  }
  return 0;
}
