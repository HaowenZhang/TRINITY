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
