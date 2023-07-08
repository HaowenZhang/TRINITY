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

int main(int argc, char **argv)
{
  float z;
  struct smf_fit the_smf;
  int i;
  // if (argc<2+NUM_PARAMS) {
  //   fprintf(stderr, "Usage: %s mass_cache (mcmc output)\n", argv[0]);
  //   exit(1);
  // }
  if (argc<2) {
    fprintf(stderr, "Usage: %s mass_cache (mcmc output)\n", argv[0]);
    exit(1);
  }


  // for (i=0; i<NUM_PARAMS; i++)
  //   the_smf.params[i] = atof(argv[i+2]);

  // Read in the model parameter values if provided by the user
  if (argc >= 2+NUM_PARAMS)
    for (i=0; i<NUM_PARAMS; i++)
      the_smf.params[i] = atof(argv[i+2]);
  // Otherwise use our default values
  else
    for (i=0; i<NUM_PARAMS; i++)
      the_smf.params[i] = param_default[i];


  the_smf.params[NUM_PARAMS] = 0;

  setup_psf(1);
  load_mf_cache(argv[1]);
  init_timesteps();
  gsl_set_error_handler_off();

  assert_model(&the_smf);
  calc_sfh(&the_smf);
  printf("#z csfr_obs csfr cbhar_obs cbhar\n");
  for (i=0; i < num_outputs; i++)
  {
    z = 1 / steps[i].scale - 1;
    printf("%f %e %e %e %e\n", z, steps[i].observed_cosmic_sfr, steps[i].cosmic_sfr, steps[i].observed_cosmic_bhar, steps[i].cosmic_bhar);
  }
  
  // for (i=0; i<50; i++) {
  //   z = i*0.05;
  //   printf("%f %g %e\n", z, log10(calc_cosmic_sfr(z)), calc_cosmic_bhar(z));
  // }
  // for (i=0; i<85; i++) {
  //   z = 2.5+i*0.2;
  //   printf("%f %g %e\n", z, log10(calc_cosmic_sfr(z)), calc_cosmic_bhar(z));
  // }

  return 0;
}
