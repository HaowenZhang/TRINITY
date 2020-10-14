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
  
  
  printf("#BH_Mass BHAR_SFR_avg\n");


  for (m=5; m<10.5; m+=0.1) 
  {
    printf("%f %e\n", m, calc_bhar_sfr_mbh(m, z));
  }
  // fclose(pfile);
  return 0;
}
