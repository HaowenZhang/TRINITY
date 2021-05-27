// Calculate AGN duty cycle as a function of halo mass and redshift.
// The duty cycle is defined as the fraction of dark matter halos that
// contain active SMBHs, regardless of their (positive) Eddington ratios.
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

// Helper functions to interpolate between halo mass bins and snapshots.
#define INTERP(y,x) { double s1, s2; s1 = (1.0-f)*steps[i].x[j]+f*steps[i+1].x[j];     \
    s2 = (1.0-f)*steps[i].x[j+1]+f*steps[i+1].x[j+1];			               \
    y = s1+mf*(s2-s1); }

#define LINTERP(y,x) { double s1, s2; s1 = (1.0-f)*log10(steps[i].x[j])+f*log10(steps[i+1].x[j]); \
    s2 = (1.0-f)*log10(steps[i].x[j+1])+f*log10(steps[i+1].x[j+1]);	\
    y = s1+mf*(s2-s1); }


int main(int argc, char **argv)
{
  int64_t i,j;
  struct smf_fit smf;
  if (argc<2+NUM_PARAMS) 
  {
    fprintf(stderr, "Usage: %s mass_cache (mcmc output)\n", argv[0]);
    exit(1);
  }
  for (i=0; i<NUM_PARAMS; i++)
    smf.params[i] = atof(argv[i+2]);
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
  double t,m;
  printf("#1+z M_h SM M_bh duty_cycle\n");
  for (t=0; t<num_outputs-1; t+=1.0/3.0) 
  {
    i = t;
    double f = t-i;
    double zp1 = (1.0-f)/steps[i].scale + f/steps[i+1].scale;
    // the redshift-dependent component of BH duty cycle.
    double dc = (1.0-f)*steps[i].smhm.bh_duty + f*steps[i+1].smhm.bh_duty;
    double dc_mbh = (1.0-f)*steps[i].smhm.dc_mbh + f*steps[i+1].smhm.dc_mbh;
    double dc_mbh_w = (1.0-f)*steps[i].smhm.dc_mbh_w + f*steps[i+1].smhm.dc_mbh_w;

    double z = zp1 - 1;
    // mass_real is the halo mass threshold above which we deem the halos
    // as too rare to care about.
    double mass_real = 13.5351-0.23712*z+2.0187*exp(-z/4.48394);
    for (m=8; m<15; m+=0.05) 
    {
      if (m >= mass_real) continue;
      double mf = (m-M_MIN)*BPDEX+0.5;
      double sm;
      double dc = (1.0-f)*steps[i].smhm.bh_duty + f*steps[i+1].smhm.bh_duty;
      int64_t j = mf;
      mf -= j;
      
      double log_bh_mass, bh_eta, log_bh_acc_rate;
      LINTERP(log_bh_mass,bh_mass_avg);
      INTERP(sm, sm_avg);
      LINTERP(log_bh_acc_rate,bh_acc_rate);
      INTERP(bh_eta,bh_eta);
      if (!isfinite(log_bh_acc_rate)) continue;

      // Apply the mass-dependent component of duty cycle.
      double dc_mass_factor = exp((log_bh_mass - dc_mbh) / dc_mbh_w);
      dc_mass_factor = dc_mass_factor / (1 + dc_mass_factor);
      dc *= dc_mass_factor;

      printf("%f %f %f %f %f\n", zp1, m, log10(sm), log_bh_mass, dc);
    }
  }
  return 0;
}

