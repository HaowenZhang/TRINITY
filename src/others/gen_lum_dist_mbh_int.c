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
  double m;
  if (argc<3+NUM_PARAMS) {
    fprintf(stderr, "Usage: %s z scatter_obs mass_cache (mcmc output)\n", argv[0]);
    exit(1);
  }
  double z = atof(argv[1]);
  double scatter_obs = atof(argv[2]);
  for (i=0; i<NUM_PARAMS; i++)
    smf.params[i] = atof(argv[i+4]);
  gsl_set_error_handler_off();
  setup_psf(1);
  load_mf_cache(argv[3]);
  init_timesteps();
  INVALID(smf) = 0;
  calc_sfh(&smf);
  int64_t step;
  double f;
  calc_step_at_z(z, &step, &f);
  
  
  printf("#BH_Alpha: %f\n", steps[step].smhm.bh_alpha);
  printf("#BH_delpha: %f\n", steps[step].smhm.bh_delta);
  double s1 = steps[step].smhm.scatter;
  double s2 = steps[step].smhm.bh_scatter;
  double s = sqrt(s1*s1+s2*s2);

  double mbh_min = steps[step].bh_mass_min;
  double mbh_max = steps[step].bh_mass_max;
  double mbh_bpdex = MBH_BINS / (steps[step].bh_mass_max - steps[step].bh_mass_min);
  double mbh_inv_bpdex = 1.0 / mbh_bpdex;


  for (i=0; i<MBH_BINS; i++)
  {
    if (mbh_min + (i + 0.5) * mbh_inv_bpdex < 8) continue;
    fprintf(stdout, "%.6f ", mbh_min + (i + 0.5) * mbh_inv_bpdex);
    for (int j=0; j<LBOL_BINS; j++)
      fprintf(stdout, "%.6e ", steps[step].lum_dist_full[i*LBOL_BINS+j]);
    fprintf(stdout, "\n");
  }
  return 0;
}
