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
  int64_t i;
  struct smf_fit smf;
  if (argc<2+NUM_PARAMS) {
    fprintf(stderr, "Usage: %s mass_cache (mcmc output) z bh_mass(in Msun) bh_Lbol(in erg/s)\n", argv[0]);
    exit(1);
  }
  for (i=0; i<NUM_PARAMS; i++)
    smf.params[i] = atof(argv[i+2]);
  double z = atof(argv[i+4]);
  double Mbh = atof(argv[i+5]);
  double Lbol = atof(argv[i+6]);

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
  printf("#z=%.2f, Mbh=%.6f, Lbol=%.6f\n", z, Mbh, Lbol);
  printf("#Note: prob = Const.*prob_Mbh*prob_eta*nd_halo, where the constant is chosen so that the summation of prob equals unity.\n");
  printf("#Mh prob prob_Mbh prob_eta nd_halo\n");
  double t,m;
  int64_t step;
  double f;
  calc_step_at_z(z, &step, &f);

  double sm_scatter = steps[step].smhm.scatter * steps[step].smhm.gamma;
  double scatter = sqrt(sm_scatter*sm_scatter 
                      + steps[step].smhm.bh_scatter*steps[step].smhm.bh_scatter);
  double probs[M_BINS] = {0};
  double probs_eta[M_BINS] = {0};
  double probs_Mbh[M_BINS] = {0};
  double total = 0;
  double eta = Lbol - 38.1 - Mbh;
  // For now ignore the interpolation between different snapshots.
  for (i = 0; i < M_BINS; i++)
  {
    double dMbh = Mbh - steps[step].log_bh_mass[i];
    double prob_Mbh = 1 / (sqrt(2*M_PI) * scatter) * exp(-dMbh*dMbh / (2*scatter*scatter)); //dimension: dlogMbh^{-1}
    
    double eta_frac = eta - steps[step].bh_eta[i];
    double prob_eta;

    if (eta_frac < steps[step].ledd_min[i] || eta_frac > steps[step].ledd_max[i]) 
      prob_eta = 0; 
    else
    {
      double bher_f = (eta_frac-steps[step].ledd_min[i])*steps[step].ledd_bpdex[i];
      int64_t bher_b = bher_f;
      bher_f -= bher_b;
      double p1 = steps[step].bher_dist[i*BHER_BINS+bher_b];
      double p2 = steps[step].bher_dist[i*BHER_BINS+bher_b+1];
      if (bher_b >= BHER_BINS-1) p2 = p1;
      prob_eta = p1 + bher_f*(p2-p1);
    }
    probs[i] = prob_eta * prob_Mbh * steps[step].t[i];
    probs_eta[i] = prob_eta;
    probs_Mbh[i] = prob_Mbh;
    // printf("Mh: %f, eta_frac: %f, eta_min: %f, eta_max: %f, prob_eta: %e, prob_Mbh: %e, nd: %e\n", 
    //   (M_MIN + (i + 0.5) * INV_BPDEX), eta_frac, steps[step].ledd_min[i], steps[step].ledd_max[i],
    //   prob_eta, prob_Mbh, steps[step].t[i]);
    total += probs[i];
  }

  if (total > 0)
  {
    for (i = 0; i < M_BINS; i++) probs[i] /= (total*INV_BPDEX); //normalize and convert it to dex^-1.
  }
  
  for (i = 0; i < M_BINS; i++) 
  {
    printf("%.1f %.6e %.6e %.6e %.6e\n", (M_MIN + (i + 0.5) * INV_BPDEX), probs[i], probs_Mbh[i], probs_eta[i], steps[step].t[i]);
  }
  


  
  
  return 0;
}
