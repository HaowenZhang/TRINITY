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
  printf("#BH_alpha: %f\n", steps[step].smhm.bh_alpha);
  printf("#BH_delpha: %f\n", steps[step].smhm.bh_delta);
  for (l=-10; l>-34.1; l-=0.1) {
    printf("%f %f\n", l, log10(calc_quasar_lf_new(l, z))); //, log10(calc_quasar_lf(l+6.25, z))-2.5);
  }

  /*printf("#Alpha Inv.Avg\n");
  for (i=0; i<10; i++) {
    float alpha = i/10.0;
    }*/
  
  /*  printf("#M_h M_bh dM_bh/dt Edd_r. Eta L_typ L_typ2\n");
  for (i=0; i<M_BINS; i++) {
    printf("%f %e %e %e %f %f\n", M_MIN+i*INV_BPDEX, steps[step].bh_mass[i], steps[step].bh_acc_rate[i], steps[step].bh_acc_rate[i]/steps[step].bh_mass[i]*4.5e7, steps[step].bh_eta[i], -5.26 -2.5*log10(steps[step].bh_acc_rate[i]*4.5e7));
    }*/
  return 0;
}
