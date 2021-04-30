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

extern int64_t num_outputs;
extern struct timestep *steps;

int main(int argc, char **argv)
{
  struct smf_fit smf;
  int i,j;
  float sfr_thresh;
  if (argc<4+NUM_PARAMS) {
    fprintf(stderr, "Usage: %s sfr_thresh mass_cache (mcmc output)\n", argv[0]);
    exit(1);
  }
  sfr_thresh = atof(argv[1]);
  for (i=0; i<NUM_PARAMS; i++)
    smf.params[i] = atof(argv[i+3]);
  smf.params[NUM_PARAMS] = 0;
  INVALID(smf) = 0;

  setup_psf(1);
  load_mf_cache(argv[2]);
  init_timesteps();  
  calc_sfh(&smf);

  double lgt = log10(sfr_thresh);

  printf("#Z ND dCounts/dz (1 sq. deg)\n");
  for (i=0; i<num_outputs; i++) {
    double nd = 0;
    double sm_scatter_corr = exp(pow(steps[i].smhm.scatter*log(10), 2)/2.0);
    double scatter = sqrt(pow(steps[i].smhm.scatter, 2.0)+0.3*0.3);
    double norm = 1.0/(sqrt(2.0*M_PI)*scatter);
    for (j=0; j<M_BINS; j++) {
      if (!isfinite(steps[i].sfr[j]) || !(steps[i].sfr[j]>0)) continue;
      double dlgsfr = (log10(steps[i].sfr[j]/sm_scatter_corr) - lgt)/scatter;
      nd += exp(-0.5*dlgsfr*dlgsfr)*steps[i].t[j];
    }
    nd *= norm;
    float z =  1.0/steps[i].scale-1.0;
    printf("%f %e %e\n", z, nd, 0.00030461741*dVc(z)*nd);
  }

  return 0;
}
