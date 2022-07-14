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
  double mbh, mh;
  if (argc<3+NUM_PARAMS) {
    fprintf(stderr, "Usage: %s z mass_cache (mcmc output)\n", argv[0]);
    exit(1);
  }
  double z = atof(argv[1]);
  for (i=0; i<NUM_PARAMS; i++)
    smf.params[i] = atof(argv[i+3]);
  gsl_set_error_handler_off();
  setup_psf(1);
  load_mf_cache(argv[2]);
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
  printf("#BH_Scatter at fixed mass: %f\n", s);
  printf("#BH_Mass ND ND(Active)\n");
  fprintf(stderr, "Mbh ND ND(11<Mh<12) ND(12<Mh<13) ND(13<Mh<14) ND(14<Mh<15)\n");

  for (mbh=4; mbh<10.5; mbh+=0.1) 
  {
    printf("%f %e", mbh, calc_bhmf(mbh, z));
    for (mh=11; mh<15; mh+=1)
    {
      printf(" %e", calc_bhmf_mh(mbh, z, mh, mh+1));
    }
    printf("\n");
  }
  return 0;
}
