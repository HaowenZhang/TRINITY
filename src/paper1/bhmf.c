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
  // if (argc<3+NUM_PARAMS) {
  //   fprintf(stderr, "Usage: %s z mass_cache (mcmc output)\n", argv[0]);
  //   exit(1);
  // }

  if (argc<3) {
    fprintf(stderr, "Usage: %s z mass_cache (mcmc output)\n", argv[0]);
    exit(1);
  }
  double z = atof(argv[1]);


  // for (i=0; i<NUM_PARAMS; i++)
  //   smf.params[i] = atof(argv[i+3]);

  // Read in the model parameter values if provided by the user
  if (argc >= 3+NUM_PARAMS)
    for (i=0; i<NUM_PARAMS; i++)
      smf.params[i] = atof(argv[i+3]);
  // Otherwise use our default values
  else
    for (i=0; i<NUM_PARAMS; i++)
      smf.params[i] = param_default[i];

  gsl_set_error_handler_off();
  setup_psf(1);
  load_mf_cache(argv[2]);
  init_timesteps();
  INVALID(smf) = 0;
  calc_sfh(&smf);
  int64_t step;
  double f;
  calc_step_at_z(z, &step, &f);
  
  // printf("#double power-law norm: %e\n", doublePL_norm(-0.6, -10, 2, NULL));
  // printf("#Schechter avg: %e\n", schechter_inv_avg(-0.6, -10, 2, NULL));
  // printf("#double power-law norm: %e\n", doublePL_norm(-0.6, -4, 2, NULL));
  // printf("#Schechter avg: %e\n", schechter_inv_avg(-0.6, -4, 2, NULL));
  printf("#BH_Alpha: %f\n", steps[step].smhm.bh_alpha);
  printf("#BH_delpha: %f\n", steps[step].smhm.bh_delta);
  double s1 = steps[step].smhm.scatter;
  double s2 = steps[step].smhm.bh_scatter;
  double s = sqrt(s1*s1+s2*s2);
  printf("#BH_Scatter at fixed mass: %f\n", s);
  printf("#BH_Mass ND ND(Active)\n");
  fprintf(stderr, "Mbh Mh Mbh_med scatter ND ND_active\n");
  // FILE * pfile;
  // pfile = fopen('./ABHMF_check.txt', 'w');
  // fprintf(pfile, "Mbh Mh ND f_active\n");
  for (m=4; m<10.5; m+=0.1) {
    printf("%f %e %e\n", m, calc_bhmf(m, z), calc_active_bhmf(m,z));
  }
  // fclose(pfile);
  return 0;
}
