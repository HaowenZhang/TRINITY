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
#include "check_syscalls.h"
#include "expcache2.h"
#include "all_luminosities.h"
#include <string.h>

extern int64_t num_outputs;
extern struct timestep *steps;


int main(int argc, char **argv)
{
  struct smf_fit the_smf;
  int64_t i;

  if (argc<3+NUM_PARAMS) {
    fprintf(stderr, "Usage: %s mass_cache z (mcmc output)\n", argv[0]);
    exit(1);
  }

  for (i=0; i<NUM_PARAMS; i++)
    the_smf.params[i] = atof(argv[i+3]);

  float z = atof(argv[2]);
  int64_t step;
  double f;
  gen_exp10cache();
  setup_psf(1);
  load_mf_cache(argv[1]);
  init_timesteps();
  INVALID(the_smf) = 0;
  calc_sfh(&the_smf);
  calc_step_at_z(z, &step, &f);

  printf("#Log10(Mvir) [Msun] Log10(SM) Log10(SM+ICL) Log10(M200c)\n");
  double scatter_corr = exp(pow(steps[step].smhm.scatter*log(10), 2)/2.0);
  for (i=0; i<M_BINS; i++) {
    float m = M_MIN + (i+0.5)/(double)BPDEX;
    if (m < 10.5) continue;
    if (m > 15) continue;
    float sm = steps[step].sm[i] + steps[step].sm_icl[i]/scatter_corr;
    printf("%f %f %f\n", m, log10(steps[step].sm[i]), log10(sm));
  }
  return 0;
}
  
