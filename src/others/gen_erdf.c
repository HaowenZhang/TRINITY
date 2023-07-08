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
  double bher;
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
  double mbh_grid[6] = {5,6,7,8,9,10};
  // printf("#Double power-law norm: %e\n", doublePL_norm(-0.6, -10, 2, NULL));
  // printf("#Schechter avg: %e\n", schechter_inv_avg(-0.6, -10, 2, NULL));
  // printf("#Double power-law norm: %e\n", doublePL_norm(-0.6, -4, 2, NULL));
  // printf("#Schechter avg: %e\n", schechter_inv_avg(-0.6, -4, 2, NULL));
  // printf("#BH_alpha: %f\n", steps[step].smhm.bh_alpha);
  // printf("#BH_delpha: %f\n", steps[step].smhm.bh_delta);
  double mbh_min = steps[step].bh_mass_min;
  double mbh_max = steps[step].bh_mass_max;
  double mbh_bpdex = MBH_BINS / (mbh_max - mbh_min);
  double mbh_inv_bpdex = 1.0 / mbh_bpdex;

  for (bher=-6; bher<4; bher+=0.1) {
    double bher_dist = 0;
    double bher_dist_mbh[5] = {0};
    for (i=0; i<MBH_BINS; i++)
    {
	double mbh = mbh_min + (i + 0.5) * mbh_inv_bpdex;
        double lum = 38.1 + mbh + bher;
        if (lum < LBOL_MIN || lum > LBOL_MAX || mbh < mbh_grid[0]) continue;
        double lbol_f = (lum - LBOL_MIN) * LBOL_BPDEX;
        int64_t lbol_b = lbol_f;
        lbol_f -= lbol_b;
        if (lbol_b >= LBOL_BINS - 1) {lbol_b = LBOL_BINS - 2; lbol_f=1;}
        double p1 = steps[step].lum_func_full[i*LBOL_BINS+lbol_b];
        double p2 = steps[step].lum_func_full[i*LBOL_BINS+lbol_b+1];
        bher_dist += p1 + lbol_f * (p2 - p1);
        int64_t mbh_bin = mbh - mbh_grid[0];
        bher_dist_mbh[mbh_bin] += p1 + lbol_f * (p2 - p1);   
    }
    printf("%f %e", bher, bher_dist); //, log10(calc_quasar_lf(l+6.25, z))-2.5);
    for (int j=0; j<5; j++) printf(" %e", bher_dist_mbh[j]);
    printf("\n");
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
