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
#include "expcache2.h"

int main(int argc, char **argv)
{
  int64_t i;
  struct smf_fit smf;
  double l;
  if (argc<3+NUM_PARAMS) {
    fprintf(stderr, "Usage: %s z mass_cache (mcmc output)\n", argv[0]);
    exit(1);
  }
  double z = atof(argv[1]);
  for (i=0; i<NUM_PARAMS; i++)
    smf.params[i] = atof(argv[i+3]);

  nonlinear_luminosity = 1;
gsl_set_error_handler_off();
  setup_psf(1);
  load_mf_cache(argv[2]);
  init_timesteps();
  INVALID(smf) = 0;
  calc_sfh(&smf);
  int64_t step;
  double f;
  calc_step_at_z(z, &step, &f);

  // printf("#Double power-law norm: %e\n", doublePL_norm(-0.6, -10, 2, NULL));
  // printf("#Schechter avg: %e\n", schechter_inv_avg(-0.6, -10, 2, NULL));
  // printf("#Double power-law norm: %e\n", doublePL_norm(-0.6, -4, 2, NULL));
  // printf("#Schechter avg: %e\n", schechter_inv_avg(-0.6, -4, 2, NULL));
  printf("#BH_Alpha: %f\n", steps[step].smhm.bh_alpha);
  printf("#BH_delpha: %f\n", steps[step].smhm.bh_delta);
  double SM_delim[14] = {5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5};
  for (i = 0; i < 13; i++)
  {
    for (l=-10; l>-34.1; l-=0.1) 
    {
    printf("%f %f %f %f\n", SM_delim[i], SM_delim[i + 1], l, log10(calc_quasar_lf_split_sm(l, z, SM_delim[i], SM_delim[i + 1]))); //, log10(calc_quasar_lf(l+6.25, z))-2.5);
     }
  }
  

  return 0;
}
