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

  printf("#SM eta Prob(eta|SM,z)\n");
  float sm, eta;
  for (sm=8; sm<13; sm+=0.1) {
    for (eta=-6; eta<2; eta += 0.2) {
      double pro = calc_qpdf_at_sBHAR_m_z(eta, sm, z);
      if (pro <= 0) pro = -1000;
      else pro = log10(pro);

      double pro_new = calc_qpdf_at_sBHAR_m_z_new(eta, sm, z);
      if (pro_new <= 0) pro_new = -1000;
      else pro_new = log10(pro_new);

      printf("%f %f %e\n", sm, eta, pro_new);
    }
  }
  return 0;
}
