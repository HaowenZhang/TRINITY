#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include "../base/observations.h"
#include "../base/smf.h"
#include "../base/all_smf.h"
#include "../base/distance.h"
#include "../base/integrate.h"
#include "../base/mlist.h"
#include "../base/calc_sfh.h"
#include "../base/param_default.h"

extern double param_default[];

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
  // if (argc<3+NUM_PARAMS) {
  //   fprintf(stderr, "Usage: %s z mass_cache (mcmc output)\n", argv[0]);
  //   exit(1);
  // }
  if (argc<3) {
    fprintf(stderr, "Usage: %s z mass_cache (mcmc output)\n", argv[0]);
    exit(1);
  }
  z = atof(argv[1]);

  // for (i=0; i<NUM_PARAMS; i++)
  //   the_smf.params[i] = atof(argv[i+3]);

  // Read in the model parameter values if provided by the user
  if (argc >= 3+NUM_PARAMS)
    for (i=0; i<NUM_PARAMS; i++)
      the_smf.params[i] = atof(argv[i+3]);
  // Otherwise use our default values
  else
    for (i=0; i<NUM_PARAMS; i++)
      the_smf.params[i] = param_default[i];


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
