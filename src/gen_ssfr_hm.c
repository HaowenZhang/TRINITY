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

#define MASS_START 7.5
#define MASS_STOP 12.5
#define MASS_BPDEX 10
#define MASS_BINS (int)((MASS_STOP-MASS_START)*(MASS_BPDEX))
#define MASS_STEP (1.0/(double)MASS_BPDEX)

extern struct timestep *steps;
extern int64_t num_outputs;

int main(int argc, char **argv)
{
  float m;
  struct smf_fit the_smf;
  int i, j;
  if (argc<3+NUM_PARAMS) {
    fprintf(stderr, "Usage: %s mass_cache (mcmc output)\n", argv[0]);
    exit(1);
  }
  for (i=0; i<NUM_PARAMS; i++)
    the_smf.params[i] = atof(argv[i+2]);
  the_smf.params[NUM_PARAMS] = 0;


  setup_psf(1);
  load_mf_cache(argv[1]);
  init_timesteps();
  calc_sfh(&the_smf);

  printf("#z log10(M_halo)[Msun] SSFR[yr^-1]\n");
  for (j=0; j<num_outputs; j++) {
    if (steps[j].scale < 1.0/9.0) continue;
    float z = 1.0/steps[j].scale - 1.0;
    for (i=0; i<M_BINS; i++) {
      m = M_MIN + (i+0.5)*INV_BPDEX;
      if (m<9 || steps[j].t[i] < 1e-9) continue;
      printf("%f %f %g\n", z, m, steps[j].sfr[i]/steps[j].sm_avg[i]);
    }
  }
  return 0;
}
