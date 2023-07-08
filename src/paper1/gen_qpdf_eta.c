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
    for (eta=-6; eta<4; eta += 0.2) {
      double pro = calc_qpdf_at_sBHAR_m_z_new(eta, sm, z);
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
