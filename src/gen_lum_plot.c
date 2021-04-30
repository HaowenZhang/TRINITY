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
#include "mah.h"
#include "check_syscalls.h"
#include "expcache2.h"
#include "all_luminosities.h"

extern int64_t num_outputs;

int main(int argc, char **argv) {
  int64_t i,filter;
  struct smf_fit the_smf;
  FILE *output;
  double m,t;
  char buffer[1024];


  if (argc<3+NUM_PARAMS) {
    fprintf(stderr, "Usage: %s mass_cache filter (mcmc output)\n", argv[0]);
    exit(1);
  }
  
  for (i=0; i<NUM_PARAMS; i++)
    the_smf.params[i] = atof(argv[i+3]);

  gen_exp10cache();
  setup_psf(1);
  load_mf_cache(argv[1]);
  init_timesteps();
  INVALID(the_smf) = 0;
  calc_sfh(&the_smf);
  
  filter = atoi(argv[2]);
  sprintf(buffer, "lum_%s.dat", lum_names[filter]);
  output = check_fopen(buffer, "w");
  enum lum_types tp = (enum lum_types) filter;
  load_luminosities("../smf_mcmc2/data/mags_ab_zall_dust.dat");
  float *lums = NULL;
  gen_all_luminosities(tp, &lums, -1, 0);
  
  for (m=9; m<15.1; m+=0.05) {
    for (t=1; t<(num_outputs-1)*3; t++) {
      double ft = (float)(((int)t)%3) / 3.0;
      i = t/3;
      float tscale = steps[i].scale + ft*(steps[i+1].scale - steps[i].scale);
      float zp1 = 1.0/tscale;
      double lum = lum_at_hm_z(lums, m, zp1-1.0);
      int64_t j = (m - M_MIN)*BPDEX;
      double new_smd = 0.05*pow(10, mf_cache(tscale, m))*steps[i].new_sm[j];
      if (lum > 30) lum = 30;
      fprintf(output, "%f %f %f %e %e\n", zp1, m, lum, new_smd, steps[i].sm[j]);
    }
  }
  fclose(output);
  return 0;
}
