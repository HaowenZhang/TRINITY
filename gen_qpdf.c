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
  if (argc<3+NUM_PARAMS) {
    fprintf(stderr, "Usage: %s z mass_cache (mcmc output)\n", argv[0]);
    exit(1);
  }
  double z = atof(argv[1]);
  for (i=0; i<NUM_PARAMS; i++)
    smf.params[i] = atof(argv[i+3]);
gsl_set_error_handler_off();
  nonlinear_luminosity = 1;
  setup_psf(1);
  load_mf_cache(argv[2]);
  init_timesteps();
  INVALID(smf) = 0;
  calc_sfh(&smf);
  int64_t step;
  double f;
  calc_step_at_z(z, &step, &f);

  printf("#SM LB Prob(LB|SM,z)\n");
  float sm, lb;
  for (sm=8; sm<13; sm+=0.1) {
    for (lb=40; lb<49; lb += 0.25) {
      double prob = calc_qpdf_at_l_m_z(lb, sm, z);
      if (prob <= 0) prob = -1000;
      else prob = log10(prob);
      printf("%f %f %f\n", sm, lb, prob);
    }
  }
  return 0;
}
