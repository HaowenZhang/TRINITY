// Calculate the galaxy star formation rates, merger rates,
// and star-forming fractions (1.0 - quenched fractions) as
// functions of halo mass and redshift.
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


#define REAL_ND_CUTOFF 1e-9

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
  
  if (argc < 3) 
  {
    fprintf(stderr, "Usage: %s mass_cache param_file (> output_file)\n", argv[0]);
    exit(1);
  }
  // Read in model parameters
  FILE *param_input = check_fopen(argv[2], "r");
  char buffer[2048];
  fgets(buffer, 2048, param_input);
  read_params(buffer, smf.params, NUM_PARAMS);

  // Fix some model parameters.
  assert_model(&smf);
  // Turn off the built-in GSL error handler that kills the program
  // when an error occurs. We handle the errors manually.
  gsl_set_error_handler_off();
  // We use non-linear scaling relation between the radiative and total Eddington ratios.
  nonlinear_luminosity = 1;
  // Set up the PSF for stellar mass functions. See observations.c.
  setup_psf(1);
  // Load cached halo mass functions.
  load_mf_cache(argv[1]);
  // Initialize all the timesteps/snapshots.
  init_timesteps();
  INVALID(smf) = 0;
  // Calculate the star-formation histories and black hole histories. See calc_sfh.c.
  calc_sfh(&smf);

  printf("#Is the model invalid? %e\n", INVALID(smf));
  double t,m;
  printf("#1+z M_h SM SFR ND smf new_sm_sf new_sm_merger sfrac lv\n");
 

  for (t=0; t<num_outputs-1; t+=1.0) 
  {
    i = t;
    double f = t-i;
    double zp1 = (1.0-f)/steps[i].scale + f/steps[i+1].scale;
    double mu = (1.0 - f) * steps[i].smhm.mu + f * steps[i].smhm.mu;
    double csfr_completeness = (1.0 - f) * steps[i].smhm.csfr_completeness + f * steps[i].smhm.csfr_completeness;
    double smloss = steps[i].smloss[i];
    double icl_frac = steps[i].smhm.icl_frac;
    for (m=7.1; m<=16.1; m+=0.05) 
    {
      double mf = (m-M_MIN)*BPDEX+0.5;
      int64_t j = mf;
      mf -= j;
      double log_sm, log_sfr, nd, sfrac, sfr, new_sm_sf, new_sm_merger, sm_icl;

      LINTERP(log_sfr, sfr);
      INTERP(sfr, sfr);
      INTERP(log_sm, log_sm);
      INTERP(nd, t);
      INTERP(sfrac, sfrac);
      //sfr = steps[i].sfr[j];
      //log_sm = log10(steps[i].sm_avg[j]);
      nd = steps[i].t[j];
      sm_icl = steps[i].sm_icl[j];
      new_sm_sf = steps[i].new_sm[j] * smloss;
      new_sm_merger = steps[i].sm_icl[j] * icl_frac;
      sfrac = steps[i].sfrac[j];
      double lv = steps[i].lv[j];
      //log_sm += mu;
      if (nd < REAL_ND_CUTOFF) continue;

      printf("%f %f %f %e %e %e %e %e %e %f\n", zp1, m, log_sm, sfr, sfrac, sm_icl, new_sm_sf, new_sm_merger, sfrac, lv);
    }
  }
  return 0;
}
