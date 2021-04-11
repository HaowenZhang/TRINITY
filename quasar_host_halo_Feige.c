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

float fitter(float *params) {
  struct smf_fit test;
  int i;
  for (i=0; i<NUM_PARAMS; i++)
    test.params[i] = params[i];

  assert_model(&test);

  for (i=0; i<NUM_PARAMS; i++)
    params[i] = test.params[i];
  //iterations++;
  float err = all_smf_chi2_err(test);
  if (!isfinite(err) || err<0) return 1e30;
  return err;
}

float calc_chi2(float *params) {
  return fitter(params);
}

int main(int argc, char **argv)
{
  int64_t i, j, k, l;
  struct smf_fit smf;
  // Calculate the halo mass distribution for a certain quasar population between (z_low, z_high), and
  // with BH mass and bolometric luminosity above bh_mass_low and bh_Lbol_low, respectively.
  if (argc<2+NUM_PARAMS) {
    fprintf(stderr, "Usage: %s mass_cache (mcmc output) z_low z_high bh_mass_low(in Msun) bh_Lbol_low(in erg/s)\n", argv[0]);
    exit(1);
  }
  for (i=0; i<NUM_PARAMS; i++)
    smf.params[i] = atof(argv[i+2]);
  double z_low = atof(argv[i+4]);
  double z_high = atof(argv[i+5]);
  double Mbh_low = atof(argv[i+6]);
  double Lbol_low = atof(argv[i+7]);

  nonlinear_luminosity = 1;
  gsl_set_error_handler_off();
  setup_psf(1);
  load_mf_cache(argv[1]);
  init_timesteps();
  INVALID(smf) = 0;
  //double chi2 = calc_chi2(smf.params);
  //printf("Actual chi2=%e\n", chi2);
  calc_sfh(&smf);
  printf("#Is the model invalid? %e\n", INVALID(smf));
  printf("#z_low=%.2f, z_high=%.2f, Mbh_low=%.6f, Lbol_low=%.6f\n", z_low, z_high, Mbh_low, Lbol_low);
  printf("#Note: prob = Const.*prob_Mbh*prob_eta*nd_halo, where the constant is chosen so that the summation of prob equals unity.\n");
  printf("#Mh prob prob_Mbh prob_eta nd_halo\n");
  double t,m;
  int64_t step_low, step_high;
  double f;
  calc_step_at_z(z_low, &step_high, &f);
  calc_step_at_z(z_high, &step_low, &f);

  double prob_Mh[M_BINS] = {0};
  
  // Calculate the probabilities integrated over
  for (i=step_low; i<=step_high; i++)
  {
    double sm_scatter = steps[i].smhm.scatter * steps[i].smhm.bh_gamma;
    double scatter = sqrt(sm_scatter*sm_scatter 
                        + steps[i].smhm.bh_scatter*steps[i].smhm.bh_scatter);
    double mbh_min = steps[i].bh_mass_min, mbh_max = steps[i].bh_mass_max;
    double mbh_inv_bpdex = (mbh_max - mbh_min) / MBH_BINS;

    double prob_lbol[MBH_BINS] = {0}; //The cumulative probabilities of being more luminous than Lbol_low,
                                      //as a function of black hole mass.

    for (k=0; k<MBH_BINS; k++)
    {
      double Mbh = mbh_min + (k + 0.5) * mbh_inv_bpdex;
      //fprintf(stderr, "i=%d, z=%f, mbh=%f\n", i, 1/steps[i].scale-1, Mbh);
      if (Mbh < Mbh_low) continue;
      double lbol_f = (Lbol_low - LBOL_MIN) * LBOL_BPDEX;
      int64_t lbol_b = lbol_f;
      lbol_f -= lbol_b;
      if (lbol_b >= LBOL_BINS - 1) {lbol_b = LBOL_BINS - 2; lbol_f = 1;}
      prob_lbol[k] = (1 - lbol_f) * steps[i].lum_dist_full[k*LBOL_BINS + lbol_b];
      for (l=lbol_b+1; l<LBOL_BINS; l++) prob_lbol[k] += steps[i].lum_dist_full[k*LBOL_BINS + l];
      //fprintf(stderr, "i=%d, z=%f, mbh=%f, prob_l=%e\n", i, 1/steps[i].scale-1, Mbh, prob_lbol[k]);
    }



    for (j=0; j<M_BINS; j++)
    {
      for (k=0; k<MBH_BINS; k++)
      {
        double Mbh = mbh_min + (k + 0.5) * mbh_inv_bpdex;
        if (Mbh < Mbh_low) continue;
        double dMbh = Mbh - steps[i].log_bh_mass[j];
        double prob_Mbh = 1 / (sqrt(2*M_PI) * scatter) * exp(-dMbh*dMbh / (2*scatter*scatter)) * mbh_inv_bpdex; //dimension: dlogMbh^{-1}
        prob_Mh[j] += prob_Mbh * prob_lbol[k];
      }
      prob_Mh[j] *= steps[i].t[j];
    }
  }

  
  for (i = 0; i < M_BINS; i++) 
  {
    printf("%.1f %.6e\n", (M_MIN + (i + 0.5) * INV_BPDEX), prob_Mh[i]);
  }
  
  return 0;
}
