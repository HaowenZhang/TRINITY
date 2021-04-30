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

void swap(double* a, double *b)
{
  double tmp = *a;
  *a = *b;
  *b = tmp;
}

int main(int argc, char **argv)
{
  int64_t i, j;
  struct smf_fit smf;
  if (argc<2+NUM_PARAMS) {
    fprintf(stderr, "Usage: %s mass_cache (mcmc output) z Mbh_low(in Msun) Mbh_high log_eta_low log_eta_high\n", argv[0]);
    exit(1);
  }
  for (i=0; i<NUM_PARAMS; i++)
    smf.params[i] = atof(argv[i+2]);
  double z = atof(argv[i+2]);
  double Mbh_low = atof(argv[i+3]);
  double Mbh_high = atof(argv[i+4]);
  double eta_low = atof(argv[i+5]);
  double eta_high = atof(argv[i+6]);

  if (Mbh_low >= Mbh_high) swap(&Mbh_low, &Mbh_high);
  if (eta_low >= eta_high) swap(&eta_low, &eta_high);
  gsl_set_error_handler_off(); 

  nonlinear_luminosity = 1;
  setup_psf(1);
  load_mf_cache(argv[1]);
  init_timesteps();
  INVALID(smf) = 0;
  //double chi2 = calc_chi2(smf.params);
  //printf("Actual chi2=%e\n", chi2);
  calc_sfh(&smf);
  printf("#Is the model invalid? %e\n", INVALID(smf));
  printf("#z=%.2f, Mbh_low=%.6f, Mbh_high=%.6f, eta_low=%.6f, eta_high=%.6f\n", z, Mbh_low, Mbh_high, eta_low, eta_high);
  // printf("#Note: prob = Const.*prob_Mbh*prob_eta*nd_halo, where the constant is chosen so that the summation of prob equals unity.\n");
  printf("#Mh prob_Mbh prob_eta nd_halo [Mpc^(-3) dex^(-1)]\n");
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

  const int64_t mbh_bins = 10000;
  double f_mass[mbh_bins];
  double Mbh_step = (Mbh_high - Mbh_low) / mbh_bins;

  for (i = 0; i < mbh_bins; i++)
  {
    double Mbh = Mbh_low + (i + 0.5) * Mbh_step;
    double fac_exp = exp((Mbh - steps[step].smhm.dc_mbh) / steps[step].smhm.dc_mbh_w);
    f_mass[i] = fac_exp / (1 + fac_exp);
  }

  double fac_to_x = 1 / (sqrt(2) * scatter);
  double norm_gauss = fac_to_x / sqrt(M_PI);

  printf("#Duty cycle factor: %.6f\n", steps[step].smhm.bh_duty);

  // double eta = Lbol - 38.1 - Mbh;
  // For now ignore the interpolation between different snapshots.
  for (i = 0; i < M_BINS; i++)
  {
    
    double prob_Mbh = 0;

    for (j = 0; j < mbh_bins; j++)
    {
      double Mbh = Mbh_low + (j + 0.5) * Mbh_step;
      double dx = (Mbh - steps[step].log_bh_mass[i]) * fac_to_x;
      prob_Mbh += norm_gauss * exp(-dx*dx) * f_mass[i] * Mbh_step * steps[step].smhm.bh_duty;
    }



    // double dMbh_low = Mbh_low - steps[step].log_bh_mass[i];
    // double dMbh_high = Mbh_high - steps[step].log_bh_mass[i];
    // double x1 = dMbh_low * fac_to_x;
    // double x2 = dMbh_high * fac_to_x;

    // double prob_Mbh; //The integrated probability of hosting a BH within the mass range
    // if (x1 < 0 && x2 < 0) prob_Mbh = 0.5 * (erf(x1) - erf(x2));
    // else if (x1 < 0 && x2 > 0) prob_Mbh = 0.5 * (erf(x1) + erf(x2));
    // else prob_Mbh = 0.5 * (erf(x2) - erf(x1));

    double eta_frac_low = eta_low - steps[step].bh_eta[i];
    double eta_frac_high = eta_high - steps[step].bh_eta[i];
    double prob_eta;

    if (eta_frac_high < steps[step].ledd_min[i] || eta_frac_low > steps[step].ledd_max[i]) 
      prob_eta = 0; 

    else
    {
      double f1 = (eta_frac_low-steps[step].ledd_min[i])*steps[step].ledd_bpdex[i];
      int64_t b1;
      if (f1 < 0)
      {
        b1 = 0;
        f1 = 0;
      }
      else
      {
        b1 = f1;
        f1 -= b1;
      }
      
      double f2 = (eta_frac_high-steps[step].ledd_min[i])*steps[step].ledd_bpdex[i];
      int64_t b2;
      if (f2 > BHER_BINS)
      {
        b2 = BHER_BINS-1;
        f2 = 1.0;
      }
      else
      {
        b2 = f2;
        f2 -= b2;
      }

      for (j=b1+1; j<b2; j++)
      {
        prob_eta += steps[step].bher_dist[i*BHER_BINS+j];
      }
      prob_eta += (1.0 - f1) * steps[step].bher_dist[i*BHER_BINS+b1];
      prob_eta += f2 * steps[step].bher_dist[i*BHER_BINS+b2];
      prob_eta /= steps[step].ledd_bpdex[i];

      
    }
    probs[i] = prob_eta * prob_Mbh * steps[step].t[i] * BPDEX;
    probs_eta[i] = prob_eta;
    probs_Mbh[i] = prob_Mbh;
    // printf("Mh: %f, eta_frac: %f, eta_min: %f, eta_max: %f, prob_eta: %e, prob_Mbh: %e, nd: %e\n", 
    //   (M_MIN + (i + 0.5) * INV_BPDEX), eta_frac, steps[step].ledd_min[i], steps[step].ledd_max[i],
    //   prob_eta, prob_Mbh, steps[step].t[i]);
    total += probs[i];
  }

  total *= INV_BPDEX;

  printf("#Total Number Density: %.6e Mpc^(-3)\n", total);
  
  for (i = 0; i < M_BINS; i++) 
  {
    printf("%.1f %.6e %.6e %.6e %.6e\n", (M_MIN + (i + 0.5) * INV_BPDEX), probs[i], probs_Mbh[i], probs_eta[i], steps[step].t[i] * BPDEX);
  }
  


  
  
  return 0;
}
