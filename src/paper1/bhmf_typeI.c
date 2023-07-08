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
#include "../base/expcache2.h"
#include "../base/param_default.h"

extern double param_default[];

int main(int argc, char **argv)
{
  int64_t i;
  struct smf_fit smf;
  double m;
  // if (argc<4+NUM_PARAMS) {
  //   fprintf(stderr, "Usage: %s z lbol_lim mass_cache (mcmc output)\n", argv[0]);
  //   exit(1);
  // }

  if (argc<4) {
    fprintf(stderr, "Usage: %s z lbol_lim mass_cache (mcmc output)\n", argv[0]);
    exit(1);
  }

  double z = atof(argv[1]);
  double lbol_lim = atof(argv[2]);

  // Read in the model parameter values if provided by the user
  if (argc >= 4+NUM_PARAMS)
    for (i=0; i<NUM_PARAMS; i++)
      smf.params[i] = atof(argv[i+4]);
  // Otherwise use our default values
  else
    for (i=0; i<NUM_PARAMS; i++)
      smf.params[i] = param_default[i];


gsl_set_error_handler_off();
  setup_psf(1);
  load_mf_cache(argv[3]);
  init_timesteps();
  INVALID(smf) = 0;
  calc_sfh(&smf);
  int64_t step;
  double f;
  calc_step_at_z(z, &step, &f);
  

  fprintf(stderr, "Mbh ND\n");
  // FILE * pfile;
  // pfile = fopen('./ABHMF_check.txt', 'w');
  // fprintf(pfile, "Mbh Mh ND f_active\n");
  for (m=4; m<11; m+=0.03) {
    printf("%f %e\n", m, calc_bhmf_typeI_new(m, z, lbol_lim));
  }
  // fclose(pfile);
  return 0;
}
