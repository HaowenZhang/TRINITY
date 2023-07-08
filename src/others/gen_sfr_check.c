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

#define MASS_START 7.5
#define MASS_STOP 12.5
#define MASS_BPDEX 10
#define MASS_BINS (int)((MASS_STOP-MASS_START)*(MASS_BPDEX))
#define MASS_STEP (1.0/(double)MASS_BPDEX)

#define INTERP(y,x) { double s1, s2; s1 = (1.0-f)*steps[i].x[j]+f*steps[i+1].x[j];     \
    s2 = (1.0-f)*steps[i].x[j+1]+f*steps[i+1].x[j+1];                    \
    y = s1+mf*(s2-s1); }

#define LINTERP(y,x) { double s1, s2; s1 = (1.0-f)*log10(steps[i].x[j])+f*log10(steps[i+1].x[j]); \
    s2 = (1.0-f)*log10(steps[i].x[j+1])+f*log10(steps[i+1].x[j+1]); \
    y = s1+mf*(s2-s1); }

int main(int argc, char **argv)
{
  struct smf_fit the_smf;
  int i, j;
  double t, m;
  if (argc<2+NUM_PARAMS) {
    fprintf(stderr, "Usage: %s mass_cache (mcmc output)\n", argv[0]);
    exit(1);
  }

  for (i=0; i<NUM_PARAMS; i++)
    the_smf.params[i] = atof(argv[i+2]);
  the_smf.params[NUM_PARAMS] = 0;


  setup_psf(1);
  load_mf_cache(argv[1]);
  init_timesteps();
  calc_sfh(&the_smf);

  printf("1+z Mh SM new_sm logV frac_SF SFR sm_icl icl_frac icl_frac_max\n");

  for (t=0; t<num_outputs-1; t+=1.0) {
    i = t;
    double f = t-i;
    double zp1 = (1.0-f)/steps[i].scale + f/steps[i+1].scale;

    for (j=0; j<M_BINS; j++) {
    // for (m=8; m<15; m+=0.05) {
      // double mf = (m-M_MIN)*BPDEX+0.5;
      // int64_t j = mf;
      // mf -= j;
      double lv, sfr, nd, sfrac, sm, new_sm, sm_icl, icl_frac, icl_frac_max;

      lv = steps[i].lv[j];
      sfr = steps[i].sfr[j];
      sfrac = steps[i].sfrac[j];
      sm = steps[i].sm_avg[j];
      new_sm = steps[i].new_sm[j];
      nd = steps[i].t[j];

      new_sm *= (steps[i].smloss[j] / (1.0 - steps[i].smhm.icl_frac));
      sm_icl = steps[i].sm_icl[j];
      icl_frac = steps[i].smhm.icl_frac * sm_icl / new_sm;
      icl_frac_max = sm_icl / new_sm;

      // INTERP(lv, lv);
      // INTERP(sfr, sfr);
      // INTERP(sfrac, sfrac);
      // INTERP(sm, sm_avg);
      // INTERP(new_sm, new_sm);
      // INTERP(nd, t);
      // //if (nd < 1e-8) continue;

      printf("%f %f %f %f %f %f %f %f %f %f\n", zp1, M_MIN + (j + 0.5) * INV_BPDEX, log10(sm), log10(new_sm), lv, log10(sfrac), log10(sfr), log10(sm_icl), icl_frac, icl_frac_max);
    }
  }

  return 0;
}
