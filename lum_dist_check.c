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



int main(int argc, char **argv)
{
  int64_t i, j;
  struct smf_fit smf;
  if (argc<2+NUM_PARAMS) {
    fprintf(stderr, "Usage: %s z mass_cache (mcmc output)\n", argv[0]);
    exit(1);
  }
  for (i=0; i<NUM_PARAMS; i++)
    smf.params[i] = atof(argv[i+3]);
  assert_model(&smf);
  gsl_set_error_handler_off();
  nonlinear_luminosity = 1;
  setup_psf(1);
  double z = atof(argv[1]);
  load_mf_cache(argv[2]);
  init_timesteps();
  INVALID(smf) = 0;
  //double chi2 = calc_chi2(smf.params);
  //printf("Actual chi2=%e\n", chi2);
  calc_sfh(&smf);

  int64_t step;
  double f;
  calc_step_at_z(z, &step, &f);

  printf("#z=%f, snapshot=%d\n", z, step);
  printf("#Mbh_min: %f\n", steps[step].bh_mass_min);
  printf("#Mbh_max: %f\n", steps[step].bh_mass_max);
  printf("#Mbh_bpdex: %f\n", (double)(MBH_BINS) / (steps[step].bh_mass_max - steps[step].bh_mass_min));
  printf("#Lbol_min: %f\n", (double)(LBOL_MIN));
  printf("#Lbol_max: %f\n", (double)(LBOL_MAX));
  printf("#Lbol_bpdex: %f\n", (double)(LBOL_BPDEX));

  fprintf(stderr, "#z=%f, snapshot=%d\n", z, step);
  fprintf(stderr, "#Mbh_min: %f\n", steps[step].bh_mass_min);
  fprintf(stderr, "#Mbh_max: %f\n", steps[step].bh_mass_max);
  fprintf(stderr, "#Mbh_bpdex: %f\n", (double)(MBH_BINS) / (steps[step].bh_mass_max - steps[step].bh_mass_min));
  fprintf(stderr, "#Lbol_min: %f\n", (double)(LBOL_MIN));
  fprintf(stderr, "#Lbol_max: %f\n", (double)(LBOL_MAX));
  fprintf(stderr, "#Lbol_bpdex: %f\n", (double)(LBOL_BPDEX));
  double mbh_bpdex = (double)(MBH_BINS) / (steps[step].bh_mass_max - steps[step].bh_mass_min);
  for (i=0; i<MBH_BINS; i++)
  {
    printf("%f ", steps[step].bh_mass_min + (i + 0.5) / mbh_bpdex);
    fprintf(stderr, "%f ", steps[step].bh_mass_min + (i + 0.5) / mbh_bpdex); 
    for (j=0; j<LBOL_BINS; j++)
    {
      printf("%e ", steps[step].lum_func_full[i*LBOL_BINS + j]);
      fprintf(stderr, "%e ", steps[step].lum_dist_full[i*LBOL_BINS + j]);
    }
    printf("\n");
    fprintf(stderr, "\n");
  } 

  return 0;
}
