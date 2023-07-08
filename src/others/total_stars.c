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

int main(int argc, char **argv)
{
  float z;
  struct smf_fit the_smf;
  int i;
  if (argc<2+NUM_PARAMS) {
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
  float total_stars = 0;
  for (i=0; i<num_outputs; i++)
    total_stars += steps[i].cosmic_sfr*steps[i].dt;
  printf("Total star density: %e\n", total_stars);
  return 0;
}
