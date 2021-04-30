// Here we only look at stellar mass of 10^12 Msun, trying to find out why it still has a
// strong amplitude.
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
  int64_t i, j;
  struct smf_fit smf;
  gsl_set_error_handler_off();
  if (argc<2+NUM_PARAMS) {
    fprintf(stderr, "Usage: %s mass_cache (mcmc output) z bh_mass(in Msun) bh_Lbol(in erg/s)\n", argv[0]);
    exit(1);
  }
  for (i=0; i<NUM_PARAMS; i++)
    smf.params[i] = atof(argv[i+2]);
  double z = atof(argv[i+2]);
  double Mbh = atof(argv[i+3]);
  double Lbol = atof(argv[i+4]);

  nonlinear_luminosity = 1;
  setup_psf(1);
  load_mf_cache(argv[1]);
  init_timesteps();
  INVALID(smf) = 0;
  //double chi2 = calc_chi2(smf.params);
  //printf("Actual chi2=%e\n", chi2);
  calc_sfh(&smf);
  printf("#Is the model invalid? %e\n", INVALID(smf));
  printf("#z=%.2f, Mbh=%.6f, Lbol=%.6f, SM=12\n", z, Mbh, Lbol);
  // printf("#Note: prob = Const.*prob_Mbh*prob_eta*nd_halo, where the constant is chosen so that the summation of prob equals unity.\n");
  printf("#Mh Mstar0 prob_Mstar prob_eta halo_ND\n");
  double t,m;
  int64_t step;
  double f;
  calc_step_at_z(z, &step, &f);

  double sm_scatter = steps[step].smhm.scatter;
  double bh_scatter = steps[step].smhm.bh_scatter;
  printf("#bh_scatter: %f\n", bh_scatter);
  double probs[SM_BINS] = {0};
  double probs_eta[M_BINS] = {0};
  // double probs_Mbh[SM_BINS] = {0};
  double probs_Mstar[M_BINS] = {0};
  double total = 0;
  double eta = Lbol - 38.1 - Mbh;
  // For now ignore the interpolation between different snapshots.

  double sm = 12;
  double bm = bulge_mass(sm, steps[step].scale);
  double log_bh_mass = calc_bh_at_bm(bm, steps[step].smhm);
  // double dMbh = Mbh - log_bh_mass;
  // double prob_Mbh = 1 / (sqrt(2*M_PI) * bh_scatter) * exp(-dMbh*dMbh / (2*bh_scatter*bh_scatter));
  // printf("sm: %.1f, prob_Mbh: %.6e\n", sm, prob_Mbh);
  for (j = 0; j < M_BINS; j++)
  {
    double dMstar = sm - steps[step].smhm.mu - steps[step].log_sm[j];
    double prob_Mstar = 1 / (sqrt(2*M_PI) * sm_scatter) * exp(-dMstar*dMstar / (2*sm_scatter*sm_scatter));
    double eta_frac = eta - steps[step].bh_eta[j];
    double prob_eta;
    if (eta_frac < steps[step].ledd_min[j] || eta_frac > steps[step].ledd_max[j]) 
      prob_eta = 0; 
    else
    {
      double bher_f = (eta_frac-steps[step].ledd_min[j])*steps[step].ledd_bpdex[j];
      int64_t bher_b = bher_f;
      bher_f -= bher_b;
      double p1 = steps[step].bher_dist[j*BHER_BINS+bher_b];
      double p2 = steps[step].bher_dist[j*BHER_BINS+bher_b+1];
      if (bher_b >= BHER_BINS-1) p2 = p1;
      prob_eta = p1 + bher_f*(p2-p1);

    }
    probs_Mstar[j] = prob_Mstar;
    probs_eta[j] = prob_eta;

    // probs[i] += prob_Mstar * prob_eta * steps[step].t[j];

  }
  // probs[i] *= prob_Mbh;
  // probs_Mbh[i] = prob_Mbh;
  // printf("SM: %.1f, probs_Mbh: %.6e\n", sm, probs_Mbh[i]);
  // total += probs[i];

  // if (total > 0)
  // {
  //   for (i = 0; i < SM_BINS; i++) probs[i] /= (total*SM_INV_BPDEX); //normalize and convert it to dex^-1.
  // }
  
  for (i = 0; i < M_BINS; i++) 
  {
    printf("%.1f %.6f %.6e %.6e %.6e\n", (M_MIN + (i + 0.5) * INV_BPDEX), steps[step].log_sm[i], probs_Mstar[i], probs_eta[i], steps[step].t[i]);
  }
  


  
  
  return 0;
}
