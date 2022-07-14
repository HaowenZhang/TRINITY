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
  if (argc<2+NUM_PARAMS) {
    fprintf(stderr, "Usage: %s mass_cache (mcmc output)\n", argv[0]);
    exit(1);
  }
  for (i=0; i<NUM_PARAMS; i++)
    smf.params[i] = atof(argv[i+2]);
gsl_set_error_handler_off();
  nonlinear_luminosity = 1;
  setup_psf(1);
  load_mf_cache(argv[1]);
  init_timesteps();
  INVALID(smf) = 0;
  calc_sfh(&smf);
  double t,m;
  printf("#1+z M_h SM SFR SM_new SM_new_SFR SM_ICL\n");
  for (i=0; i<num_outputs; i++) {
    double zp1 = steps[i].scale;
    double mu = steps[i].smhm.mu;
    double smloss = steps[i].smloss[i];
    double dt = steps[i].dt;


    for (int j=0; j<M_BINS; j++) {
      m = M_MIN + (j + 0.5) * INV_BPDEX;
      double  log_bm, log_sfr, log_sm, log_new_sm;
      double log_sm_icl;
      double icl_frac = steps[i].smhm.icl_frac;
      log_sfr = log10(steps[i].sfr[j]);
      log_new_sm = steps[i].new_sm[j];
      log_sm = log10(steps[i].sm_avg[j]);
      log_sm_icl = log10(steps[i].sm_icl[j]);
      double log_new_sm_sfr = log_sfr + log10(dt * smloss);

      //LINTERP(log_sfr, sfr);
      //LINTERP(log_new_sm, new_sm);
      //LINTERP(log_sm, sm_avg);
      //LINTERP(log_sm_icl, sm_icl);
      //log_sm += mu;
      log_new_sm = log10(log_new_sm + (icl_frac * exp10(log_sm_icl))); //In calc_sfh.c the smloss factor was canceled, so here we put it back.

      
	// float years = scale_to_years(1.0 / zp1);
      printf("%f %f %f %f %f %f %f\n", zp1, m, log_sm, log_sfr, log_new_sm, log_new_sm_sfr, log_sm_icl);
    }
  }
  return 0;
}
