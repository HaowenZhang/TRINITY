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

#define MASS_START 7
#define MASS_STOP 12.5
#define MASS_BPDEX 10
#define MASS_BINS (int)((MASS_STOP-MASS_START)*(MASS_BPDEX))
#define MASS_STEP (1.0/(double)MASS_BPDEX)


int main(int argc, char **argv)
{
  int64_t i, j, k, l;
  struct smf_fit smf;
  // Calculate the black hole mass function within a certain luminosity bin (lbol_low, lbol_high) at a certain redshift.
  if (argc<2+NUM_PARAMS) {
    fprintf(stderr, "Usage: %s mass_cache (mcmc output) z lbol_low(in erg/s) lbol_high(in erg/s)\n", argv[0]);
    exit(1);
  }
  for (i=0; i<NUM_PARAMS; i++)
    smf.params[i] = atof(argv[i+2]);
  double z = atof(argv[i+4]);
  double Lbol_low = atof(argv[i+5]);
  double Lbol_high = atof(argv[i+6]);

  nonlinear_luminosity = 1;
  gsl_set_error_handler_off();
  setup_psf(1);
  load_mf_cache(argv[1]);
  init_timesteps();
  INVALID(smf) = 0;
  //double chi2 = calc_chi2(smf.params);
  //printf("Actual chi2=%e\n", chi2);
  calc_sfh(&smf);
  assert_model(&smf); 

  printf("#Is the model invalid? %e\n", INVALID(smf));
  printf("#z=%.2f, Lbol_low=%.6f, Lbol_high=%.6f\n", z, Lbol_low, Lbol_high);
  printf("#Mbh Phi(Mbh)\n");
  double t,m;
  int64_t step;
  double f;
  calc_step_at_z(z, &step, &f);

  double prob_lbol[MBH_BINS] = {0}; //The probabilities of having luminosities between (Lbol_low, Lbol_high), as a function of BH mass.
  // double host_smf[MASS_BINS] = {0};


  double sm_scatter = steps[step].smhm.scatter * steps[step].smhm.bh_gamma;
  double scatter = sqrt(sm_scatter*sm_scatter 
                      + steps[step].smhm.bh_scatter*steps[step].smhm.bh_scatter);
  double mbh_min = steps[step].bh_mass_min, mbh_max = steps[step].bh_mass_max;
  double mbh_inv_bpdex = (mbh_max - mbh_min) / MBH_BINS;

  for (i=0; i<MBH_BINS; i++)
  {
    double mbh = mbh_min + (i + 0.5) * mbh_inv_bpdex;
    double lbol_f_low = (Lbol_low - LBOL_MIN) * LBOL_BPDEX;
    double lbol_f_high = (Lbol_high - LBOL_MIN) * LBOL_BPDEX;
    int64_t lbol_b_low = lbol_f_low; lbol_f_low -= lbol_b_low;
    int64_t lbol_b_high = lbol_f_high; lbol_f_high -= lbol_b_high;
    if (lbol_b_low >= LBOL_BINS - 1) {lbol_b_low = LBOL_BINS - 2; lbol_f_low = 1;}
    if (lbol_b_high >= LBOL_BINS - 1) {lbol_b_high = LBOL_BINS - 2; lbol_f_high = 1;}
    double prob_lum = (1 - lbol_f_low) * steps[step].lum_dist_full[i*LBOL_BINS+lbol_b_low];
    double prob_tot = 0;
    for (j=0; j<LBOL_BINS; j++) prob_tot += steps[step].lum_dist_full[i*LBOL_BINS+j];
    for (j=lbol_b_low+1; j<lbol_b_high; j++) prob_lum += steps[step].lum_dist_full[i*LBOL_BINS+j];
    prob_lum += lbol_f_high * steps[step].lum_dist_full[i*LBOL_BINS+lbol_b_high];
    
    double dc = steps[step].smhm.bh_duty;
    double f_mass = exp((mbh - steps[step].smhm.dc_mbh) / steps[step].smhm.dc_mbh_w);
    f_mass = f_mass / (1 + f_mass);
    dc *= f_mass;
    if (dc < 1e-4) dc = 1e-4;

    prob_lum /= prob_tot;
    prob_lbol[i] = prob_lum * dc;
    fprintf(stderr, "%f %e\n", mbh, prob_lbol[i]);
  }

  for (i=0; i<M_BINS; i++) 
  {
    m = M_MIN + (i + 0.5) * INV_BPDEX; //halo mass.
    double mbh_med = steps[step].log_bh_mass[i];
    double p_l = 0;

    for (j=0; j<MBH_BINS; j++)
    {
      double mbh = mbh_min + (j + 0.5) * mbh_inv_bpdex;
      double dmbh = (mbh - mbh_med) / scatter;
      // The probability of having a BH mass.
      double prob_tmp = 1 / sqrt(2*M_PI) / scatter * exp(-0.5*dmbh*dmbh) * mbh_inv_bpdex;
      p_l += prob_tmp * prob_lbol[j];
      //fprintf(stderr, "Mh=%.1f, SM=%.6f, mbh_med=%.6f, mbh=%.6f, scatter=%.6f, dmbh=%.6f, prob_tmp=%.6e\n",
        //     m, steps[step].log_sm[i], steps[step].log_bh_mass[i], mbh, scatter, dmbh, prob_tmp);
    }
    //fprintf(stderr, "Mh=%.1f, SM=%.6f, prob=%.6e\n", m, steps[step].log_sm[i], p_l);
    printf("%f %e\n", m, steps[step].t[i] * BPDEX * p_l);

  }




  
  // // Calculate the probabilities integrated over
  // for (i=step_low; i<=step_high; i++)
  // {
    

  //   double prob_lbol[MBH_BINS] = {0}; //The cumulative probabilities of being more luminous than Lbol_low,
  //                                     //as a function of black hole mass.

  //   for (k=0; k<MBH_BINS; k++)
  //   {
  //     double Mbh = mbh_min + (k + 0.5) * mbh_inv_bpdex;
  //     if (Mbh < Mbh_low) continue;
  //     double lbol_f = (Lbol_low - LBOL_MIN) * LBOL_BPDEX;
  //     int64_t lbol_b = lbol_f;
  //     lbol_f -= lbol_b;
  //     if (lbol_b >= LBOL_BINS - 1) {lbol_b = LBOL_BINS - 2; lbol_f = 1;}
  //     prob_lbol[k] = (1 - lbol_f) * steps[i].lum_dist_full[k*LBOL_BINS + lbol_b];
  //     for (l=lbol_b+1; l<LBOL_BINS; l++) prob_lbol[k] += steps[i].lum_dist_full[k*LBOL_BINS + l];
  //   }

  //   for (j=0; j<M_BINS; j++)
  //   {
  //     for (k=0; k<MBH_BINS; k++)
  //     {
  //       double Mbh = mbh_min + (k + 0.5) * mbh_inv_bpdex;
  //       if (Mbh < Mbh_low) continue;
  //       double dMbh = Mbh - steps[i].log_bh_mass[j];
  //       double prob_Mbh = 1 / (sqrt(2*M_PI) * scatter) * exp(-dMbh*dMbh / (2*scatter*scatter)) * mbh_inv_bpdex; //dimension: dlogMbh^{-1}
  //       prob_Mh[j] += prob_Mbh * prob_lbol[k];
  //     }
  //     prob_Mh[j] *= steps[i].t[j];
  //   }
  // }

  
  // for (i = 0; i < M_BINS; i++) 
  // {
  //   printf("%.1f %.6e\n", (M_MIN + (i + 0.5) * INV_BPDEX), prob_Mh[i]);
  // }
  
  return 0;
}
