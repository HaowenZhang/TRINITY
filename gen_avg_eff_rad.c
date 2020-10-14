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
#include "universe_time.h"
//#include "fitter.h"

#define INTERP(y,x) { double s1, s2; s1 = (1.0-f)*steps[i].x[j]+f*steps[i+1].x[j];     \
    s2 = (1.0-f)*steps[i].x[j+1]+f*steps[i+1].x[j+1];			               \
    y = s1+mf*(s2-s1); }

#define LINTERP(y,x) { double s1, s2; s1 = (1.0-f)*log10(steps[i].x[j])+f*log10(steps[i+1].x[j]); \
    s2 = (1.0-f)*log10(steps[i].x[j+1])+f*log10(steps[i+1].x[j+1]);	\
    y = s1+mf*(s2-s1); }

float fitter(float *params) {
  struct smf_fit test;
  int i;
  for (i=0; i<NUM_PARAMS; i++)
    test.params[i] = params[i];

  assert_model(&test);

  for (i=0; i<NUM_PARAMS; i++)
    params[i] = test.params[i];
  //iterations++;
  float err = all_smf_chi2_err(test);
  if (!isfinite(err) || err<0) return 1e30;
  return err;
}

float calc_chi2(float *params) {
  return fitter(params);
}

int main(int argc, char **argv)
{
  int64_t i;
  struct smf_fit smf;
  if (argc<2+NUM_PARAMS) {
    fprintf(stderr, "Usage: %s mass_cache (mcmc output)\n", argv[0]);
    exit(1);
  }
  for (i=0; i<NUM_PARAMS; i++)
    smf.params[i] = atof(argv[i+2]);

  nonlinear_luminosity = 1;
  setup_psf(1);
  load_mf_cache(argv[1]);
  init_timesteps();
  INVALID(smf) = 0;
  //double chi2 = calc_chi2(smf.params);
  //printf("Actual chi2=%e\n", chi2);
  calc_sfh(&smf);

  printf("#sn_id z total_accretion eff_rad\n");

  double weight = 0;
  double tot = 0;
  double bh_density = 0;
  double bh_unmerged = 0;

  for (i=0; i<num_outputs-1; i++) 
  {
    // i = t;
    // double f = t-i;
    // double zp1 = (1.0-f)/steps[i].scale + f/steps[i+1].scale;
    // double mu = (1.0 - f) * steps[i].smhm.mu + f * steps[i].smhm.mu;

    double eff_rad = steps[i].smhm.bh_efficiency_rad;
    double dt = steps[i].dt;
    double acc_tot = 0;


    for (int64_t j = 0; j < M_BINS; j++) 
    {
      double m;
      double bh_acc_rate = steps[i].bh_acc_rate[j] == 1e-8 ? 0 : steps[i].bh_acc_rate[j];
      acc_tot += bh_acc_rate * steps[i].t[j] * dt;
    }



    weight += acc_tot;
    tot += acc_tot * eff_rad;

    printf("%d %f %f %f\n", i, 1 / steps[i].scale - 1, acc_tot, eff_rad);

    if (i == num_outputs-2)
    {
      for (int64_t j = 0; j < M_BINS; j++) 
      {
        bh_density += steps[i].bh_mass_avg[j] * steps[i].t[j];
        bh_unmerged += steps[i].bh_unmerged[j] * steps[i].t[j];
      }
    }

  }

  printf("The mass averaged radiative efficiency is: %.6f\n", tot / weight);
  printf("The cosmic central BH density: : %.3e\n", bh_density);
  printf("The cosmic unmerged BH density: : %.3e\n", bh_unmerged);
  return 0;
}
