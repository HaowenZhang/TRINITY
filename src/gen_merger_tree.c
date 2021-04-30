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
  int64_t i, j, k;
  struct smf_fit smf;
  if (argc<2) {
    fprintf(stderr, "Usage: %s mass_cache\n", argv[0]);
    exit(1);
  }
  // for (i=0; i<NUM_PARAMS; i++)
  //   smf.params[i] = atof(argv[i+2]);

  nonlinear_luminosity = 1;
  setup_psf(1);
  load_mf_cache(argv[1]);
  init_timesteps();
  
  // printf("#1+z yrs M_h M_b M_bh SFR dM_bh/dt Edd_r. bh_eta Obs.ER L_typ BH_Merge_Rate f_active acc_rate_obs\n");
 
 // for (t=0; t<num_outputs-1; t++)
 //  {
 //    for (i = 0; i < M_BINS; ++i)
 //    {
 //      printf("%f %f %f\n", 1 / steps[t].scale, M_MIN + (i + 0.5) / BPDEX, log10(steps[t].bh_mass_avg[i]));
 //    }
 //  }
  FILE *f_mf = fopen("/media/hzhang/DataBackup/UM_DR1/Trinity/HaloMF.dat", "w");
  FILE *f_mr = fopen("/media/hzhang/DataBackup/UM_DR1/Trinity/HaloMR.dat", "w");
  fprintf(f_mf, "#num_snapshot scale nd[cMpc^(-3)/bin], BPDEX=%d\n", BPDEX);
  fprintf(f_mr, "#num_snapshot scale mr[cMpc^(-3)/bin^2], BPDEX=%d\n", BPDEX);
  for (i=0; i<num_outputs; i++) 
  {
    fprintf(f_mf, "%d %.6f", i, steps[i].scale);
    fprintf(f_mr, "%d %.6f", i, steps[i].scale);
    for (j=0; j<M_BINS; j++)
    {
      fprintf(f_mf, " %.6e", steps[i].t[j]);
      for (k=0; k<M_BINS; k++)
        fprintf(f_mr, " %.6e", steps[i].merged[k*M_BINS+j]);
    }
    fprintf(f_mf, "\n");
    fprintf(f_mr, "\n");
  }
  fclose(f_mf);
  fclose(f_mr);
  return 0;
}
