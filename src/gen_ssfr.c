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

int main(int argc, char **argv)
{
  float z, m;
  struct smf_fit the_smf;
  int i;
  if (argc<3+NUM_PARAMS) {
    fprintf(stderr, "Usage: %s z mass_cache (mcmc output)\n", argv[0]);
    exit(1);
  }
  z = atof(argv[1]);
  for (i=0; i<NUM_PARAMS; i++)
    the_smf.params[i] = atof(argv[i+3]);
  the_smf.params[NUM_PARAMS] = 0;

  gsl_set_error_handler_off();
  setup_psf(1);
  load_mf_cache(argv[2]);
  init_timesteps();
  calc_sfh(&the_smf);

  for (i=0; i<MASS_BINS; i++) {
    m = MASS_START + i*MASS_STEP;
    printf("%f %g\n", m, log10(calc_ssfr(m, z)));
  }

  return 0;
}
