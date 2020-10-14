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



int main(int argc, char **argv)
{
  int64_t i;
  struct smf_fit smf;
  gsl_set_error_handler_off();
  if (argc<2) {
    fprintf(stderr, "Usage: %s mass_cache\n", argv[0]);
    exit(1);
  }
  // FILE *input = fopen(argv[2], "r");
  char buffer[2048];

  
  
  // for (i=0; i<NUM_PARAMS; i++)
  //   smf.params[i] = atof(argv[i+2]);

  // nonlinear_luminosity = 1;
  setup_psf(1);
  //load_mf_cache(argv[1]);
  init_mcmc_from_args(argc, argv);
  fprintf(stderr, "#chi2_type0...9 chi2_smhm chi2_bhar chi2_sfr chi2_icl chi2_rad bhz0 bhz1\n");
  while (fgets(buffer, 2048, stdin))
  {
    //puts(buffer);
    read_params(buffer, smf.params, NUM_PARAMS);
    INVALID(smf) = 0;
    //init_timesteps();
    // calc_sfh(&smf);
    
   all_smf_chi2_err_write(smf);
  }
  
  //double chi2 = calc_chi2(smf.params);
  //printf("Actual chi2=%.3f\n", chi2);
  
  return 0;
}
