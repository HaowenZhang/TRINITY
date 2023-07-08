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
  int64_t i, j;
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

  int64_t n = 177;
  for (i=0; i<M_BINS; i++)
  {
    double mbh_unmerged_sum, mbh_unmerged_tot;
    mbh_unmerged_sum = 0;
    for (j=0; j<MBH_BINS; j++)
    {
      double mbh_um = exp10(MBH_MIN + (j + 0.5) * MBH_INV_BPDEX);
      //printf("i=%d, j=%d, bh_unmerged_dist[i*MBH_BINS+j]=%.6e\n", i, j, steps[n].bh_unmerged_dist[i*MBH_BINS+j]);
      if (isfinite(steps[n].bh_unmerged_dist[i*MBH_BINS+j])) mbh_unmerged_sum += steps[n].bh_unmerged_dist[i*MBH_BINS+j] * mbh_um;
    }
    mbh_unmerged_tot = steps[n].t[i] * steps[n].bh_unmerged_avail[i];
    printf("M: %.2f, Summation: %.6e, Total: %.6e, total/nd: %.6e\n", M_MIN + (i + 0.5) * INV_BPDEX, mbh_unmerged_sum, mbh_unmerged_tot, steps[n].bh_unmerged_avail[i]);
  }

  return 0;
}
