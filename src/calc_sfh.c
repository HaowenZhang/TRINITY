#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include "mlist.h"
#include "smf.h"
#include "all_smf.h"
#include "universe_time.h"
#include "calc_sfh.h"
#include "mf_cache.h"
#include "expcache2.h"
#include "observations.h"
#include "check_syscalls.h"
#include "smoothing.h"
#include "distance.h"
#include "mah.h"
#include <omp.h>
#include <assert.h>

#define LogG -5.36456
#define SSFR_AVG_Q 2.5299810704438612E-12


//Notes: new ICL model!  icl_fract changed to return 1.0
//smf.h changed to include icl_m
//calc_sfh changed to include model based on ICL_m
//Altered ICL priors

extern int no_z_scaling;
extern float z_max;
extern float z_min;

static inline double doexp10(double x) 
{
  double a = exp(M_LN10*x);
  return a;
}

static inline double dolog10(double x) 
{
  if (x <= 0) return -1000;
  double a = log10(x);
  return a;
}

#define REAL_ND_CUTOFF 1e-9

extern struct timestep *steps;
extern int64_t num_outputs;
extern int no_obs_scatter;

extern double frac_below8[MBH_BINS];
extern double frac_below11[MBH_BINS];


void check_bhsm(struct smf_fit *f, int n)
{
  int i;
  #pragma omp for schedule(dynamic,5)
    for (i=0; i<M_BINS; i++)
    {
      if (steps[n].scale > 0.08 && (i > 0 && steps[n].log_bh_mass[i] > 0 && steps[n].log_bh_mass[i] < steps[n].log_bh_mass[i-1]))
      {
//            fprintf(stderr, "BHSM is not monotonic. a: %f, m: %f, mstar: %f, mbh: %f, mbh_i-1: %f, n=%d, i=%d\n", steps[n].scale, 
  //          steps[n].med_hm_at_a[i], steps[n].log_sm[i], steps[n].log_bh_mass[i], steps[n].log_bh_mass[i-1], n, i);
            INVALIDATE(f, "BHSM is not monotonic.");
      }
    }
}

void calc_sfh(struct smf_fit *f) 
{
  int64_t i,j;
  for (i=0; i<num_outputs; i++) 
  {
    if (!no_z_scaling || ((steps[i].scale < 1.0/(z_min+1.0)) &&
			  (steps[i].scale > 1.0/(z_max+1.0))))
    {
      calc_sm_hist(i, f);
      check_bhsm(f, i);
    }
    
  }

#pragma omp for schedule(guided,5)
  for (j=1; j<num_outputs; j++) {
    if (!no_z_scaling || ((steps[j].scale < 1.0/(z_min+1.0)) &&
			  (steps[j].scale > 1.0/(z_max+1.0)))) {
      calc_smf_and_ssfr(j, f);
      calc_uvlf(j, f);
      calc_total_sfr(j);

      calc_bh_acc_rate_distribution(j, f);
      calc_bh_lum_distribution_full(j, f);
      calc_active_bh_fraction(j, f);
      // calc_avg_eta_rad(j);
      //calc_total_bhar(j);
      // for (i=0; i<M_BINS; i++)
      //   calc_supEdd_frac_lum(j, i, 45);
    }
  }
  for (j=1; j<num_outputs; j++)
  {
    if (1 / steps[j].scale - 1 >= 6 && steps[j].cosmic_sfr < steps[j-1].cosmic_sfr)
    {
      //fprintf(stderr, "Increasing CSFR at z>6.\n");
      INVALIDATE(f, "decreasing CSFR at z>6.");
      break;
    }
    
  }
}


void calc_supEdd_frac_lum(int n, int i, double l)
{
  double ledd_min = steps[n].ledd_min[i];
  // double ledd_max = steps[n].ledd_max[i];
  double bpdex = steps[n].ledd_bpdex[i];

  
  double scatter_tot = sqrt(steps[n].smhm.bh_scatter * steps[n].smhm.bh_scatter +
                            steps[n].smhm.scatter * steps[n].smhm.bh_gamma * 
                            steps[n].smhm.scatter * steps[n].smhm.bh_gamma);
  double mbh_min = steps[n].log_bh_mass[i] - 8 * scatter_tot;
  double mbh_max = steps[n].log_bh_mass[i] + 8 * scatter_tot;
  double dmbh = 0.05;
  double norm_gauss = 1 / sqrt(2 * M_PI * scatter_tot * scatter_tot);
  double bher_norm = 0;

  for (int j=0; j<BHER_BINS; j++) bher_norm += steps[n].bher_dist[i * BHER_BINS + j];

  double prob_edd1 = 0;

  if (!bher_norm) 
  {
    steps[n].frac_supEdd_l[i] = 0;
    return;
  }

  double frac_above_l = 0.0;
  for (double mbh=mbh_min; mbh<=mbh_max; mbh+=dmbh)
  {
    double frac_above_l_tmp = 0.0;
    double eta_frac_l = l - 38.1 - mbh - steps[n].bh_eta[i];
    double bher_f_l = (eta_frac_l-ledd_min)*bpdex;
    int64_t bher_b_l = bher_f_l;
    bher_f_l -= bher_b_l;
    if (bher_b_l < 0) frac_above_l_tmp = 1.0;
    else if (bher_b_l > BHER_BINS - 1) frac_above_l_tmp = 0.0;
    else
    {
      frac_above_l_tmp = (1 - bher_f_l) * steps[n].bher_dist[i*BHER_BINS + bher_b_l];
      for (int j = bher_b_l + 1; j < BHER_BINS; j++) frac_above_l_tmp += steps[n].bher_dist[i*BHER_BINS + j];
      frac_above_l_tmp /= bher_norm;
    }

    double dmbh_in_sigma = (mbh - steps[n].log_bh_mass[i]) / scatter_tot;
    double prob = norm_gauss * exp(-0.5 * dmbh_in_sigma * dmbh_in_sigma);

    frac_above_l += frac_above_l_tmp * prob * dmbh;

    double eta_frac_edd1 = 0.0 - steps[n].bh_eta[i];
    double frac_above_edd1 = 0.0;

    double bher_f_edd1 = (eta_frac_edd1-ledd_min)*bpdex;
    int64_t bher_b_edd1 = bher_f_edd1;
    bher_f_edd1 -= bher_b_edd1;


    if (bher_f_edd1 + bher_b_edd1 <= bher_f_l + bher_b_l)
    {
      frac_above_edd1 = frac_above_l_tmp;
    }
    else
    {
      if (bher_b_edd1 < 0) frac_above_edd1 = 1.0;
      else if (bher_b_edd1 > BHER_BINS - 1) frac_above_edd1 = 0.0;
      
      else
      {
        frac_above_edd1 = (1 - bher_f_edd1) * steps[n].bher_dist[i*BHER_BINS + bher_b_edd1];
        for (int j = bher_b_edd1 + 1; j < BHER_BINS; j++) frac_above_edd1 += steps[n].bher_dist[i*BHER_BINS + j];
        frac_above_edd1 /= bher_norm;
      }
    } 
    prob_edd1 += prob * frac_above_edd1 * dmbh;

  }
  steps[n].frac_supEdd_l[i] = prob_edd1 / frac_above_l;
}

// This function removes the merger contribution to SMBH growth. Note that
// with this function, we are NOT assuming no mergers take place, but instead
// that mergers did take place, but we remove them to see how much does accretion
// contribute to different statistics, like BH mass function.
void remove_bh_merger()
{
  for (int i=1; i<num_outputs; i++)
  {
    for (int j=0; j<M_BINS; j++)
    {
      // Ignore halos that are too rare.
      if (!steps[i].t[j]) continue;
      
      if (isfinite(steps[i].new_bh_mass[j] * steps[i].bh_merge_rate[j] / \
                                (steps[i].bh_acc_rate[j] + steps[i].bh_merge_rate[j])))
        // Give back the merged BH mass to the unmerged (or wandering) BH reservoir.
        steps[i].bh_unmerged[j] += steps[i].new_bh_mass[j] * steps[i].bh_merge_rate[j] / \
                                  (steps[i].bh_acc_rate[j] + steps[i].bh_merge_rate[j]);
     
      if (isfinite(steps[i].bh_acc_rate[j] / \
                                (steps[i].bh_acc_rate[j] + steps[i].bh_merge_rate[j])))
        // New BH mass would be re-calculated to account for accretion only.
        steps[i].new_bh_mass[j] *= steps[i].bh_acc_rate[j] / \
                                  (steps[i].bh_acc_rate[j] + steps[i].bh_merge_rate[j]);

      // Update the average BH mass based on the old and new BH masses.
      double new_bh_mass_avg = steps[i].old_bh_mass[j] + steps[i].new_bh_mass[j];

      // The offset between the average vs. median BH mass, due to the log-normal scatter.
      double offset = log10(steps[i].bh_mass_avg[j]) - steps[i].log_bh_mass[j];

      // Update the median BH mass accordingly.
      steps[i].log_bh_mass[j] += log10(new_bh_mass_avg / steps[i].bh_mass_avg[j]);

      // In case the median BH mass is not finite, we assign an unphysically small value as a flag,
      // and mark other quantities as unphysical as well.
      if (!isfinite(steps[i].log_bh_mass[j])) steps[i].log_bh_mass[j] = -5;
      steps[i].bh_mass_avg[j] = new_bh_mass_avg;
      if (!finite(steps[i].bh_mass_avg[j])) steps[i].bh_mass_avg[j] = exp10(-5 + offset);
      if (!finite(steps[i].bh_unmerged[j])) steps[i].bh_unmerged[j] = 0;
      steps[i].bh_merge_rate[j] = 0;

    }

    // Since we updated the average and median BH masses for the i-th snapshot,
    // old BH mass for the (i+1)-th snapshot should be updated accordingly.
    if (i < num_outputs-1)
    {
      for (int j=0; j<M_BINS; j++)
      {
        if (!steps[i].t[j]) continue;
        steps[i+1].old_bh_mass[j] = 0;
        steps[i+1].bh_unmerged[j] = 0;
        for (int k=0; k<M_BINS; k++)
        {
          if (!steps[i].t[k]) continue;
          int bin = k*M_BINS + j;
          steps[i+1].old_bh_mass[j] += steps[i].bh_mass_avg[k]*steps[i+1].mmp[bin];
          steps[i+1].bh_unmerged[j] += steps[i+1].merged[bin]*(steps[i].bh_mass_avg[k] +
                              steps[i].bh_unmerged[k]) + steps[i+1].mmp[bin]*steps[i].bh_unmerged[k];
        }
        steps[i+1].old_bh_mass[j] /= steps[i+1].t[j];
        steps[i+1].bh_unmerged[j] /= steps[i+1].t[j];
      }
    }
  }
}

void calc_total_sfr(int n) 
{
  int64_t i;
  double sfr = 0, obs_sfr = 0;
  double mu = steps[n].smhm.mu;
  double kappa = steps[n].smhm.kappa;
  double z = 1.0 / steps[n].scale - 1;
  double sm_corr = pow(10, mu + kappa * exp(-(z - 2) * (z - 2) * 0.5));
  double sfr_minimum = BOUWENS_SFR_LIMIT;
  if (steps[n].scale > 1.0/(1.0+4.0))
    sfr_minimum *= (1.0/steps[n].scale - 0.9)/4.1;
  for (i=0; i<M_BINS; i++) 
  {
    sfr += steps[n].sfr[i]*steps[n].t[i];
    // if (steps[n].sfr[i]<sfr_minimum) continue;
    if (steps[n].sfr[i] > 0) 
    {
      float weight = 0.5+0.5*erf(log10(steps[n].sfr[i]/sfr_minimum)/(sqrt(2)*steps[n].smhm.scatter));
      obs_sfr += weight * steps[n].sfr[i]*steps[n].t[i]*sm_corr;
    }
  }
  steps[n].cosmic_sfr = sfr;
  steps[n].observed_cosmic_sfr = obs_sfr*steps[n].smhm.csfr_completeness;
}

void calc_total_bhar(int n)
{
 int64_t i;
 double bhar = 0;
 double bhar_obs = 0;

 double l_min = 41.5 - 2.5 * log10(steps[n].scale); // redshift-dependent lower limit of 2-10KeV QLF
 // double l_min = 42.5;
 double BC = 10.83 * pow(pow(10, l_min) / 3.86e43, 0.28) + 6.08 * pow(pow(10, l_min) / 3.86e43, -0.02);
 l_min += log10(BC); //Convert to bolometric
l_min = 90 - 2.5 * l_min; //Convert it to magnitude
 for (i = 0; i < M_BINS; i++)
 {
   if (steps[n].t[i] == 0) continue;
   //if (steps[n].bh_mass_avg[i] < 1e5) continue;

   // double eta = -(l_min + 5.26) / 2.5 - steps[n].log_bh_mass[i];
  double eta = -(l_min + 5.26) / 2.5 - log10(steps[n].bh_mass_avg[i]);
   double eta_frac = eta - steps[n].bh_eta[i];
   bhar += steps[n].bh_acc_rate[i] * steps[n].t[i] * steps[n].bh_f_occ[i];


   // if (eta_frac > steps[n].ledd_max[i]) continue;
   // if (eta_frac < steps[n].ledd_min[i]) 
   //  {
   //    bhar_obs += steps[n].bh_acc_rate[i] * steps[n].t[i] * steps[n].bh_f_occ[i];
   //  }
   // else
   // {
   //  double bher_f = (eta_frac - steps[n].ledd_min[i]) * steps[n].ledd_bpdex[i];
   //  int64_t bher_b = bher_f;
   //  bher_f -= bher_b;
   //  double p1 = 0, p2 = 0;
   //  if (i*BHER_BINS+bher_b+1 < M_BINS * BHER_BINS)
   //  {
   //    p1 = steps[n].bher_dist_full[i*BHER_BINS+bher_b];
   //    p2 = steps[n].bher_dist_full[i*BHER_BINS+bher_b+1];
   //  }
   //  else if (i*BHER_BINS+bher_b < M_BINS * BHER_BINS)
   //  {
   //    p1 = steps[n].bher_dist_full[i*BHER_BINS+bher_b];
   //    p2 = p1;
   //  }
    
   //  // if (bher_b >= BHER_BINS-1) p2 = p1;
   //  double prob = p1 + bher_f*(p2-p1);
   //  for (;bher_b + 1< BHER_BINS; bher_b++) prob += steps[n].bher_dist_full[i*BHER_BINS+bher_b+1];
   //  double total = 0;
   //  for (bher_b = 0; bher_b < BHER_BINS; bher_b++) total += steps[n].bher_dist_full[i*BHER_BINS+bher_b];
   //    //fprintf(stderr, "prob=%f, total=%f\n", prob, total);
   //  if (prob > 0) prob /= total;
    

   //  // double dc = steps[n].bh_duty[i];
   //  // if (dc < 1e-4) dc = 1e-4;
   //  // prob *= dc;
   //  bhar_obs += steps[n].bh_acc_rate[i] * steps[n].t[i] * prob;
   // }


 }
 steps[n].cosmic_bhar = bhar;
 steps[n].observed_cosmic_bhar = bhar;
}

// Also calculates observable duty cycle.
void calc_observed_bhar(int n)
{
 int64_t i;

 double l_min = 41.5 - 2.5 * log10(steps[n].scale); // redshift-dependent lower limit of 2-10KeV QLF
 // double l_min = 42.5;
 double BC = 10.83 * pow(pow(10, l_min) / 3.86e43, 0.28) + 6.08 * pow(pow(10, l_min) / 3.86e43, -0.02);
 l_min += log10(BC); //Convert to bolometric
l_min = 90 - 2.5 * l_min; //Convert it to magnitude
 for (i = 0; i < M_BINS; i++)
 {
   if (steps[n].t[i] == 0) continue;

   // double eta = -(l_min + 5.26) / 2.5 - steps[n].log_bh_mass[i];
   double dc = steps[n].bh_duty[i];
    if (dc < 1e-4) dc = 1e-4;
  double eta = -(l_min + 5.26) / 2.5 - log10(steps[n].bh_mass_avg[i]);
   double eta_frac = eta - steps[n].bh_eta[i];
   //fprintf(stderr, "eta_frac=%e\n", eta_frac);
   if (eta_frac > steps[n].ledd_max[i]) 
    {
      steps[n].bh_acc_rate_obs[i] = 0;
      steps[n].dc_obs[i] = 0;
    }
   if (eta_frac < steps[n].ledd_min[i]) 
    {
      steps[n].bh_acc_rate_obs[i] = steps[n].bh_acc_rate[i];
      steps[n].dc_obs[i] = dc;
    }
   else
   {
    double bher_f = (eta_frac - steps[n].ledd_min[i]) * steps[n].ledd_bpdex[i];
    int64_t bher_b = bher_f;
    bher_f -= bher_b;
    double p1 = steps[n].bher_dist_full[i*BHER_BINS+bher_b];
    double p2 = steps[n].bher_dist_full[i*BHER_BINS+bher_b+1];
    if (bher_b >= BHER_BINS-1) p2 = p1;
    double prob = p1 + bher_f*(p2-p1);
    for (;bher_b + 1< BHER_BINS; bher_b++) prob += steps[n].bher_dist_full[i*BHER_BINS+bher_b+1];
    double total = 0;
    for (bher_b = 0; bher_b < BHER_BINS; bher_b++) total += steps[n].bher_dist_full[i*BHER_BINS+bher_b];
      //fprintf(stderr, "prob=%f, total=%f\n", prob, total);
    if (prob > 0) prob /= total;
    steps[n].bh_acc_rate_obs[i] = steps[n].bh_acc_rate[i] * prob;
    steps[n].dc_obs[i] = dc * prob;
   }
 }
}

struct smf smhm_at_z(double z, struct smf_fit f) 
{
  struct smf c;
  double a = 1.0/(1.0+z);
  double a1 = -z/(1.0+z);
  double incompleteness;

  double z8 = z;
  double mass_real = 13.5351-0.23712*z+2.0187*exp(-z/4.48394);
  c.bin_real = (mass_real - M_MIN) * BPDEX - 1;

  if (z8 > 15) z8 = 15;
  a = 1.0/(1.0+z);
  a1 = a - 1.0;

  c.v_1 = M_1(f) + (a1*M_1_A(f) + log(1 + z) * M_1_A2(f)) + z8 * M_1_A3(f);
  c.sm_0 = 0;
  c.epsilon = doexp10(EFF_0(f) + EFF_0_A(f)*a1 + EFF_0_A2(f) * log(1 + z) + EFF_0_A3(f)*z8);
  c.alpha = ALPHA(f) + (a1*ALPHA_A(f) + log(1 + z) * ALPHA_A2(f)) + z8 * ALPHA_A3(f);
  c.delta = DELTA(f);
  c.beta = BETA(f) + a1*BETA_A(f) + z8*BETA_A2(f);
  c.gamma = 0;
  //c.lambda = doexp10(LAMBDA(f) + (a1*LAMBDA_A(f) + z8*LAMBDA_A2(f)));
  c.mu = MU(f) + a1*MU_A(f);
  c.kappa = KAPPA(f) + a1*KAPPA_A(f);
  c.passive_mass = 10.3 + z*0.5 - c.mu;
  c.scatter = SCATTER(f) + z*SCATTER_A(f);
  if (c.scatter > 0.4) c.scatter = 0.4;
  c.scatter_corr = exp(pow(c.scatter*log(10), 2)/2.0);
  c.obs_scatter = (no_obs_scatter) ? 0 : SIGMA_CENTER + SIGMA_Z(f)*(z);
  if (c.obs_scatter > 0.3) c.obs_scatter = 0.3;

  c.icl_frac = exp10(ICL_FRAC(f) + a1 * ICL_FRAC_E(f) + z8 * ICL_FRAC_Z(f));

  // incompleteness = BURST_DUST_AMP(f)/(1+exp(BURST_DUST_Z(f)-z));
  // if (z < 1.0) 
  incompleteness = 0;
  // if (z > 1.0) incompleteness -= BURST_DUST_AMP(f)/(1+exp(BURST_DUST_Z(f)-1.0));
  // if (incompleteness < 0) incompleteness = 0;
  c.sm_completeness = 1.0 - incompleteness;
  c.csfr_completeness = 1.0;
  c.sfr_sm_corr = 1.0 + (4.0*RHO_05(f)-3.23)*a + (2.46-4.0*RHO_05(f))*a*a;
  // if (c.csfr_completeness < 0.01) c.csfr_completeness = 0.01;
  // if (c.sm_completeness < 0.01) c.sm_completeness = 0.01;
  // c.ssfr_corr = 1.0/(1.0-BURST_DUST(f)*incompleteness);
  c.ssfr_corr = 1.0;
  c.sm_max = c.sm_min = 0;
  // c.f_1 = dolog10(2)-c.delta*pow(dolog10(2), c.gamma)/(1.0+exp(1));
  c.lm_slope = 1.0/c.alpha - 1.0;
  c.mpk = c.mu+c.kappa;
  c.combined_scatter = 0;
  c.valid = 1;

  c.qm = QM_0(f) + QM_1(f) *a1 + QM_2(f) *z8;
  c.qwidth = QWIDTH_0(f) + QWIDTH_1(f) *a1 + QWIDTH_2(f) * z8;

  // BH part
  c.bh_beta = BH_BETA_0(f) +BH_BETA_1(f)*a1 + BH_BETA_2(f)*z8;
  c.bh_gamma = BH_GAMMA_0(f) + BH_GAMMA_1(f)*a1 + BH_GAMMA_2(f)*z8;
  c.f_merge_bh = exp10(BH_MERGE_F_0(f) + BH_MERGE_F_1(f)*a1);
  // c.bh_merge_width = BH_MERGE_W_0(f) + BH_MERGE_W_1(f)*a1;
  c.bh_alpha = BH_ALPHA_0(f) + BH_ALPHA_1(f)*a1;
  c.bh_delta = BH_DELTA_0(f) + BH_DELTA_1(f)*a1;

  c.f_occ_min = exp10(F_OCC_MIN_0(f) + F_OCC_MIN_1(f)*log(1+z));
  c.bh_duty_m = (BH_DUTY_M_0(f) + BH_DUTY_M_1(f)*log(1+z));
  c.bh_duty_alpha = (BH_DUTY_ALPHA_0(f) + BH_DUTY_ALPHA_1(f)*log(1+z));
  // c.bh_duty_min = exp10(BH_DUTY_MIN_0(f) + BH_DUTY_MIN_1(f)*log(1+z));
  //fprintf(stderr, "z=%f, bh_duty_min=%e\n", z, c.bh_duty_min);

  c.abhmf_shift = ABHMF_SHIFT(f);
  c.bh_efficiency_rad = exp10(BH_EFFICIENCY_0(f) + a1 * BH_ETA_CRIT_1(f));
  c.bh_eta_crit = BH_ETA_CRIT_0(f);
  // c.bh_duty = BH_DUTY_0(f)+BH_DUTY_1(f)*a1;
  // if (c.bh_duty < 1e-4) c.bh_duty = 1e-4;
  c.bh_scatter = BH_SCATTER_0(f) + BH_SCATTER_1(f)*a1;
  c.rho_bh = RHO_BH_0(f) + a1*RHO_BH_1(f) + z8*RHO_BH_2(f);
  c.log_bh_scatter_corr = c.bh_gamma*c.bh_gamma*c.scatter*c.scatter + c.bh_scatter*c.bh_scatter;
  c.log_bh_scatter_corr = 0.5 * M_LN10 * c.log_bh_scatter_corr;

  c.dc_mbh = DC_MBH_0(f) + DC_MBH_1(f) * log(1+z);
  c.dc_mbh_w = DC_MBH_W_0(f) + DC_MBH_W_1(f) * log(1+z);
  c.eta_mu = ETA_MU_0(f);

  return c;
}

// Calculate MEDIAN BH mass as a function of bulge mass and redshift.
double calc_bh_at_bm(double bm, struct smf c) 
{
  return c.bh_beta + c.bh_gamma * (bm - 11.0);
}

// Calculate MEDIAN BH mass as a function of velocity dispersion and redshift.
double calc_bh_at_vdisp(double vd, struct smf c) {
  return c.bh_beta + c.bh_gamma * (vd - 2.30103);
}

// Calculate the MEDIAN stellar mass as a function of halo mass and redshift,
// by simple linear interpolations (in log10 units)
double calc_sm_at_m(double m, struct timestep steps) 
{
  double mf = (m - M_MIN) * BPDEX;
  int64_t mb;

  if (mf >= M_BINS - 1)
  {
    mb = M_BINS - 2;
    mf = 1;
  }
  else
  {
    mb = mf;
    mf -= mb;
  }

  double sm = steps.log_sm[mb] + mf * (steps.log_sm[mb+1] - steps.log_sm[mb]);
  return sm;
}

// Calculate MEDIAN SFR as a function of maximum circular velocity of halos and redshift.
// See Section 2.2 of Zhang et al. (2021)
double calc_sfr_at_lv(double m, double lv, struct smf c) 
{
  double vd = lv - c.v_1;
  double vd2 = vd/c.delta;
  double sfrac = 1 / (1 + exp((lv - c.qm) / c.qwidth));
  return (sfrac*c.epsilon * (1.0/(doexp10(c.alpha*vd) + doexp10(c.beta*vd)) + c.gamma*exp(-0.5*vd2*vd2)));
}

// Calculate galaxy quenched fraction as a function of maximum circular velocity of halos and redshift.
// See Section 2.2 of Zhang et al. (2021)
double calc_sfrac_at_m(double m, struct smf c) {
  double sfrac = 1 / (1 + exp((m - c.qm) / c.qwidth));
  return (sfrac);
}

// Calculate the intrinsic galaxy stellar mass function (SMF) 
// at a given halo mass and redshift: Phi(sm) = Phi(mh(sm)) * dlogmh / dlogsm.
double calc_smf_at_m(double m, double sm, int64_t n, int64_t n_spline, struct smf_fit *fit, gsl_interp_accel *ga) 
{
  // Ignore the stellar masses that exceed the range spanned by the interpolation.
  if (sm < steps[n].smhm.sm_min) return 1e-17;
  if (sm > steps[n].smhm.sm_max) return 1e-17;
  // Skip the calculation if the model is invalid to start with.
  if (!steps[n].smhm.valid) return 1e-17;

  if (m < M_MIN || m > M_MAX) return 1e-17;

  // Calculate the Jacobian: dlogmh/dlogsm
  double dlm_dlsm;
  int err = 0; 
  if (n_spline == 1) err = gsl_spline_eval_deriv_e(steps[n].spline, sm, ga, &dlm_dlsm);
  else if (n_spline == 2) err = gsl_spline_eval_deriv_e(steps[n].spline2, sm, ga, &dlm_dlsm);
  if (err || !isfinite(dlm_dlsm))
  {
    return 0;
  }
  dlm_dlsm = fabs(dlm_dlsm);

  // dndlogm is the halo mass function (HMF), obtained from the cache initialized at the
  // beginning.
  double dndlogm = mf_cache(steps[n].scale, m);
  double phi = doexp10(dndlogm)*dlm_dlsm;
  if (!isfinite(phi)) phi = 0;
  return phi;
}

// Calculate the galaxy stellar mass function (SMF) as a function
// of stellar mass and redshift.
double calc_smf_at_sm(int64_t n, double sm) 
{
  // Ignore the stellar masses that exceed the range spanned by the interpolation.
  if (sm < steps[n].smhm.sm_min) return 1e-17;
  if (sm > steps[n].smhm.sm_max) return 1e-17;
  // Skip the calculation if the model is invalid to start with.
  if (!steps[n].smhm.valid) return 1e-17;
  gsl_interp_accel ga = {0};

  // Find the halo mass that matches the stellar mass.
  double m;
  int err = gsl_spline_eval_e(steps[n].spline, sm, &ga, &m);
  // return 0 if GSL interpolation breaks down.
  if (err)
  {
    return 0;
  }

  // Based on the calculated halo mass, calculate the SMF.
  double result = calc_smf_at_m(m, sm, n, 1, NULL, &ga);

  // Reset the interpolation accelerator to avoid weird random bugs.
  gsl_interp_accel_reset(&ga);

  // Repeat the process if we have 2-segment scaling relation between stellar mass and halo mass.
  if (steps[n].alloc2smf)
  {
    err = gsl_spline_eval_e(steps[n].spline2, sm, &ga, &m); //Second segment.
    if (!err) result += calc_smf_at_m(m, sm, n, 2, NULL, &ga);
  }
  
  return result;
}



// Calculate the observed galaxy UV luminosity function (UVLF) 
// at a given halo mass and redshift: Phi(M_UV) = Phi(mh(M_UV)) * dlogmh / dM_UV.
double calc_uvlf_at_m(double m, double uv, int64_t n, int64_t n_spline, struct smf_fit *fit, gsl_interp_accel *ga) 
{
  // Ignore too bright/faint magnitudes, invalid models, or too small/big halo masses.
  if (uv < steps[n].smhm.uv_min) return 1e-17;
  if (uv > steps[n].smhm.uv_max) return 1e-17;
  if (!steps[n].smhm.valid) return 1e-17;
  if (m < M_MIN || m > M_MAX) return 1e-17;

  // Calculate the Jacobian: dlogmh / dM_UV
  double dlm_dMuv = 0;
  int err = 0;
  if (n_spline == 1) err = gsl_spline_eval_deriv_e(steps[n].spline_uv, uv, ga, &dlm_dMuv);
  else if (n_spline == 2) err = gsl_spline_eval_deriv_e(steps[n].spline_uv2, uv, ga, &dlm_dMuv);
  dlm_dMuv = fabs(dlm_dMuv);

  if (err || !isfinite(dlm_dMuv))
  {
    return 0;
  }

  // Calculate the halo mass function
  double dndlogm = mf_cache(steps[n].scale, m);
  // Phi(M_UV) = Phi(mh(M_UV)) * dlogmh / dM_UV.
  double phi = doexp10(dndlogm)*dlm_dMuv;
  if (!isfinite(phi)) phi = 0;
  return phi;
}

// Calculate the galaxy UV luminosity function (UVLF) as a function
// of UV magnitude and redshift.
double calc_uvlf_at_uv(int64_t n, double uv) 
{
  // Ignore too bright/faint magnitudes or invalid models.
  if (uv < steps[n].smhm.uv_min || 
      uv > steps[n].smhm.uv_max ||
      !steps[n].smhm.valid) 
    return 1e-17;


  gsl_interp_accel ga = {0};
  // Calculate the median halo mass given the UV magnitude
  double m;
  int err = gsl_spline_eval_e(steps[n].spline_uv, uv, &ga, &m);
  if (err)
  {
    return 0;
  }
  // Calculate UVLF by converting the halo mass function with Jacobian: dlogmh / dM_UV
  double result = calc_uvlf_at_m(m, uv, n, 1, NULL, &ga);
  // Reset the interpolation accelerator to avoid random weird bugs
  gsl_interp_accel_reset(&ga);

  // Repeat the process if we have a 2-segment scaling relation between M_UV and halo mass.
  if (steps[n].alloc2uvlf)
  {
    err = gsl_spline_eval_e(steps[n].spline_uv2, uv, &ga, &m); //Remember to add the second segment.
    if (!err) result += calc_uvlf_at_m(m, uv, n, 2, NULL, &ga);
  }
  return result;
}

// Calculate the active BH fraction, defined as the fraction of BHs above
// the Eddington ratio == 0.01. This definition is taken from Schulze & Wisotzki (2010)
// and Schulze et al. (2015).
void calc_active_bh_fraction(int n, struct smf_fit *fit) 
{
  int64_t i, j;
  for (i=0; i<M_BINS; i++) 
  {
    // ledd_min and ledd_max are the minimum and maximum ***SCALED*** Eddington ratios,
    // relative to the typical Eddington ratio, steps[n].bh_eta[i]
    double ledd_min = steps[n].ledd_min[i];
    double ledd_max = steps[n].ledd_max[i];

    // bins_per_dex and its reciprocal for the Eddington ratio distributions.
    double bpdex = steps[n].ledd_bpdex[i];

    // Scaled critical Eddington ratio (log10(0.01) == -2) relative to 
    // the typical value.
   double eta_frac = -2.0 - (steps[n].bh_eta[i] + log10(steps[n].bh_f_occ[i])); 

    // If the scaled Eddington ratio is too small or too large,
    // we can directly assign zero or one to the active fraction.
    if (eta_frac > ledd_max || !steps[n].t[i])
    {
      steps[n].f_active[i] = 0;
      continue;
    }
    else if (eta_frac < ledd_min)
    {
      steps[n].f_active[i] = 1;
      continue;
    } 
    
    // For all the other cases, we should traverse the whole Eddington
    // ratio distribution to count active BHs.
    double bher_f = (eta_frac-ledd_min)*bpdex;
    int64_t bher_b = bher_f;
    bher_f -= bher_b;

    double p_good = (1 - bher_f) * steps[n].bher_dist[i*BHER_BINS + bher_b]; //The area of the histogram with eta > -2
    double p_total = steps[n].bher_dist_norm[i]; //Total area

    if (!isfinite(p_good) || !isfinite(p_total) || p_total == 0)
    {
      steps[n].f_active[i] = 0;
      continue;
    }

    // Add up the Eddington ratio bins that are active
    for (j = bher_b + 1; j < BHER_BINS; j++)
    {
      p_good += steps[n].bher_dist[i*BHER_BINS + j];
    }

    // The Eddington ratio distributions are only for those that are active.
    // To account for the fact that not every SMBH is active, we need to 
    // fold in the factor of duty cycle.
    double dc = steps[n].bh_duty[i];
    if (dc < 1e-4) dc = 1e-4;
    // f_active is the fraction of active SMBHs among ***HOST*** halos.
    steps[n].f_active[i] = (p_good / p_total) * dc / steps[n].bh_f_occ[i];
  }
}

// Calculate the active BH fraction, defined as the fraction of BHs above
// certain either Eddington ratio or bolometric lumionsity. 

void calc_active_bh_fraction_lim(int n, struct smf_fit *fit, int ledd_or_lum, double lim) 
// Parameters:
  // n: the number of snapshot
  // fit: the smf model that are used in the calculation.
  // ledd_or_lum: flag that indicates whether the limit is in Eddington ratio or luminosity.
  //              ledd_or_lum = 1 if the limit is in Eddington ratios,
  //              ledd_or_lum = 0 if the limit is in luminosities.
  // lim: the lower limit in Eddington ratio or luminosity to calculate the active fraction.
{
  int64_t i, j;
  for (i=0; i<M_BINS; i++) {
    // Calculate the limit Eddington ratio. If the limit is already in Eddington ratio,
    // we don't need to do anything. Otherwise, we should calculate the corresponding
    // limit Eddington ratio based on the limit luminosity.
    double ledd_lim = ledd_or_lum ? lim : lim - 38.1 - steps[n].log_bh_mass[i];

    // For the code below, see void calc_active_bh_fraction() for detailed comments.
    double ledd_min = steps[n].ledd_min[i];
    double ledd_max = steps[n].ledd_max[i];

    double bpdex = steps[n].ledd_bpdex[i];

    double eta_frac = ledd_lim - steps[n].bh_eta[i];

    if (eta_frac < ledd_min || eta_frac > ledd_max) 
    {
      steps[n].f_active[i] = 0;
      continue;
    } 

    double bher_f = (eta_frac-ledd_min)*bpdex;
    int64_t bher_b = bher_f;
    bher_f -= bher_b;

    double p_good = (1 - bher_f) * steps[n].bher_dist[i*BHER_BINS + bher_b]; //The area of the histogram with eta > -2
    double p_total = steps[n].bher_dist_norm[i]; //Total area

    if (!isfinite(p_good) || !isfinite(p_total) || p_total == 0)
    {
      steps[n].f_active[i] = 0;
      continue;
    }

    for (j = bher_b + 1; j < BHER_BINS; j++)
    {
      p_good += steps[n].bher_dist[i*BHER_BINS + j];
    }

    double dc = steps[n].bh_duty[i];
    if (dc < 1e-4) dc = 1e-4;

    steps[n].f_active[i] = p_good / p_total * dc;
  }
}

// Convert the average BH Eddington ratio into the ***TYPICAL*** Eddington ratio
// (or the breaking point Eddington ratio of double power-law Eddington ratio
// distributions) for each halo mass bin, and calculate BH Eddington ratio 
// distributions for all mass bins. See Section 2.6 of Zhang et al. (2021).
void calc_bh_acc_rate_distribution(int n, struct smf_fit *fit) 
{
  int64_t i, j, k;
  memset(steps[n].bher_dist, 0, sizeof(double)*BHER_BINS*M_BINS);
  memset(steps[n].bher_dist_norm, 0, sizeof(double)*M_BINS);
  steps[n].ledd_min_abs = 100; 
  steps[n].ledd_max_abs = -100;

  double bher_min = BHER_EFF_MIN;

  // The critical total (or radiative) Eddington ratio where the scaling between total
  // and radiative Eddington ratios start to change.
  double bh_eta_crit = steps[n].smhm.bh_eta_crit;

  double (*bher_prob)(double, double, double, double, double, double) = &_prob_of_ledd_linear;
  if (nonlinear_luminosity) bher_prob = &_prob_of_ledd_nonlinear;

  for (i=0; i<M_BINS; i++)
  {
    // mass-dependent modulation of duty cycle. Note here that we ASSUME that the duty cycle is dependent on
    // the average BH mass of ***ALL*** halos, not SMBH host halos. Since bh_mass_avg is now the average BH mass
    // for host halos, we have to give back the bh_f_occ factor.
    double dc = (steps[n].bh_duty[i]);
    if (dc < 1e-4) dc = 1e-4;

    // Note that here we are going to calculate integrals of 1 / (x**a + x**b) dx and 1 / (x**(a+1) + x**(b+1)) dx,
    // while the doublePL_norm is implemented to calculate the norm assuming dP/dlogx = 1 / (x**a + x**b), so there is
    // an additional -1 offset in the following power-law indices!!!!!!
    double nom = doublePL_norm(steps[n].smhm.bh_alpha, steps[n].smhm.bh_delta, bher_min, BHER_EFF_MAX, fit);
    double dnom = doublePL_norm(steps[n].smhm.bh_alpha - 1, steps[n].smhm.bh_delta - 1, bher_min, BHER_EFF_MAX, fit);
    // Note that the duty cycle is defined as the fraction of active SMBHs among ***ALL*** halos, but we need the fraction
    // of active SMBHs among all SMBHs in the denominator of the correction factor. So a factor of 1 / bh_f_occ is needed.
    // double bh_eta_corr = log10(nom/dnom/(dc / steps[n].bh_f_occ[i]));
    double bh_eta_corr = log10(nom/dnom/(dc));
    //double bh_eta_corr = log10(nom/dnom/(dc));
    steps[n].bh_eta[i] += bh_eta_corr;


    // Update the minimum/maximum ABSOULTE RADIATIVE EDDINGTON ratios for this snapshot.
    double ledd_min = BHER_EFF_MIN + steps[n].bh_eta[i];
    double ledd_max = BHER_MAX + steps[n].bh_eta[i];
    if (ledd_min < steps[n].smhm.bh_eta_crit)
    {
      ledd_min = 2 * (ledd_min - 0.5 * steps[n].smhm.bh_eta_crit);
    }
    if (ledd_max < steps[n].smhm.bh_eta_crit)
    {
      ledd_max = 2 * (ledd_max - 0.5 * steps[n].smhm.bh_eta_crit);
    }
    if (ledd_min < steps[n].ledd_min_abs) steps[n].ledd_min_abs = ledd_min;
    if (ledd_max > steps[n].ledd_max_abs) steps[n].ledd_max_abs = ledd_max;

    ledd_min = steps[n].bh_eta[i] + BHER_MIN;
    //ledd_min = steps[n].bh_eta[i] - log10(steps[n].bh_f_occ[i]) + BHER_MIN;
    if (nonlinear_luminosity && ledd_min < bh_eta_crit)
    ledd_min = (ledd_min - 0.5 * bh_eta_crit)*2.0;
    ledd_min -= steps[n].bh_eta[i];
    //ledd_min -= (steps[n].bh_eta[i] - log10(steps[n].bh_f_occ[i]));

    //double ledd_eff_min = steps[n].bh_eta[i] - log10(steps[n].bh_f_occ[i]) + BHER_EFF_MIN;
    double ledd_eff_min = steps[n].bh_eta[i] + BHER_EFF_MIN;
    if (nonlinear_luminosity && ledd_eff_min < bh_eta_crit)
    ledd_eff_min = (ledd_eff_min - 0.5 * bh_eta_crit)*2.0;
    ledd_eff_min -= steps[n].bh_eta[i];
    //ledd_eff_min -= (steps[n].bh_eta[i] - log10(steps[n].bh_f_occ[i]));

    ledd_max = steps[n].bh_eta[i] + BHER_MAX;
    if (nonlinear_luminosity && ledd_max < bh_eta_crit)
    ledd_max = (ledd_max - 0.5 * bh_eta_crit)*2.0;
    ledd_max -= steps[n].bh_eta[i];

    double bpdex = ((double)BHER_BINS-1)/(ledd_max - ledd_min);
    double inv_bpdex = 1.0/bpdex;


    steps[n].ledd_min[i] = ledd_min;
    steps[n].ledd_eff_min[i] = ledd_eff_min;
    steps[n].ledd_max[i] = ledd_max;
    steps[n].ledd_bpdex[i] = bpdex;


    int64_t jmin = (ledd_eff_min-ledd_min)*bpdex;
    int64_t jmax = (ledd_max-ledd_min)*bpdex;

    for (j=jmin; j<jmax; j++) 
    {
      float bher = ledd_min+j*inv_bpdex;
      float prob = bher_prob(bher+steps[n].bh_eta[i], steps[n].bh_eta[i], steps[n].smhm.bh_alpha, steps[n].smhm.bh_delta, 1.0, bh_eta_crit);
      steps[n].bher_dist[i*BHER_BINS+j] = prob;
      steps[n].bher_dist_norm[i] += prob;

    }

    double total = 0;
    total = steps[n].bher_dist_norm[i] * inv_bpdex;
    if (total > 0)
      total = 1.0/total;
    for (k=0; k<BHER_BINS; k++) steps[n].bher_dist[i*BHER_BINS+k] *= total; 
    steps[n].bher_dist_norm[i] = bpdex;

  }
  steps[n].bh_mass_min = 2;
  steps[n].bh_mass_max = 11;
}

// Calculate the luminosity functions/distributions for all ***BH mass*** bins.
void calc_bh_lum_distribution_full(int n, struct smf_fit *fit)
{
  int64_t i, j, k;
  memset(steps[n].lum_dist_full, 0, sizeof(double) * MBH_BINS * LBOL_BINS);
  memset(steps[n].lum_func_full, 0, sizeof(double) * MBH_BINS * LBOL_BINS);

  double mbh_bpdex = MBH_BINS / (steps[n].bh_mass_max - steps[n].bh_mass_min);
  double mbh_inv_bpdex = 1.0 / mbh_bpdex;

  // The total scatter in BH mass at fixed ***halo mass*** is a quadratic sum
  // of the scatter in BH mass at fixed ***stellar mass*** and that in stellar mass
  // at fixed halo mass, scaled by the slope of the black hole mass--bulge mass relation.
  double sm_scatter = steps[n].smhm.scatter*steps[n].smhm.bh_gamma;
  double scatter = sqrt(sm_scatter*sm_scatter + steps[n].smhm.bh_scatter*steps[n].smhm.bh_scatter);
  // Pre-calculate the normalization factor before the exponential in Gaussian distribution.
  double gauss_norm = 1 / sqrt(2 * M_PI) / scatter;
  
  for (i=0; i<MBH_BINS; i++)
  {
    double mbh = steps[n].bh_mass_min + (i + 0.5) * mbh_inv_bpdex;
    
    double tnd = 0; //Total ND of BHs with this mass, including dormant and active BHs.
                    //This is gonna be used to calculate the luminosity distributions 
                    //(not functions) of active BHs, and also the total BH mass functions.
    // Traverse all the halo mass bins to count their contribution to each BH mass bin.
    for (j=0; j < M_BINS; j++)
    {
      
      // The mass-dependent component of AGN duty cycle.
      double dc = steps[n].bh_duty[j];
      if (dc < 1e-4) dc = 1e-4;


      // The difference between the median BH mass of the halo mass bin
      // and the BH mass we're interested in.
      double dmbh = (mbh - steps[n].log_bh_mass[j]) / scatter;
      
      // w_mbh is the fraction of the halos that host BHs of mass mbh.
      double w_mbh = gauss_norm * exp(-0.5 * dmbh * dmbh) * mbh_inv_bpdex; //Here I chose not to include
                                                            //mbh_inv_bpdex, which is "dlogMbh"
                                                          //in the integral of QLF/QPDF.
                                                      //The reason is that in this way
                                                  //the calculation of total luminosity function
                                              //is simply the summation of lum_func_full array along the Mbh direction,
                                              //because the width of each bin is already accounted for in
                                              //mbh_inv_bpdex; And the calculation of QPDF would be the weighted
                                              //summation of lum_dist_full along the Mbh direction, with the weight
                                              //being EXACTLY (i.e., including the normalization) the gaussian PDF of Mbh as given by the BHBM.

      // So the contribution to the BH mass functions from this halo mass bin
      // is simply w_mbh * halo number density (steps[n].t[j]) * SMBH occupation 
      // fraction (steps[n].bh_f_occ[j]).
      tnd += w_mbh * steps[n].t[j] * steps[n].bh_f_occ[j];

      // Calculate the halo mass bin's contribution to every bolometric luminosity bin.
      for (k=0; k<LBOL_BINS; k++)
      {
        double lbol = LBOL_MIN + (k + 0.5) * LBOL_INV_BPDEX;
        // Convert the luminosity into Eddington ratio...
        double bh_eta = lbol - 38.1 - mbh;
        // ... and scale it with the typical Eddington ratio.
        double eta_frac = bh_eta - (steps[n].bh_eta[j] + log10(steps[n].bh_f_occ[j]) 
          - (steps[n].smhm.rho_bh - 1) * (mbh - steps[n].log_bh_mass[j] - steps[n].smhm.log_bh_scatter_corr));
        // If the fractional (or scaled) Eddington ratio is too small/big or even
        // not a finite number, just skip.
        if ((!isfinite(eta_frac)) || eta_frac < steps[n].ledd_min[j] || eta_frac > steps[n].ledd_max[j])
	      {  
          continue;
        } 

        // Otherwise we just do linear interpolations.
        double bher_f = (eta_frac-steps[n].ledd_min[j])*steps[n].ledd_bpdex[j];
        int64_t bher_b = bher_f;
        bher_f -= bher_b;

        double p1 = steps[n].bher_dist[j*BHER_BINS+bher_b];
        double p2 = steps[n].bher_dist[j*BHER_BINS+bher_b+1];
        if (bher_b >= BHER_BINS-1) p2 = p1;
        // Note that the steps[n].bher_dist is not normalized to match the duty cycle,
        // so we need to put the duty cycle (dc) into the calculation of nd_l below.
        // We also have to account for the occupation fraction.
        double nd_l = (p1 + bher_f*(p2-p1)) * w_mbh * steps[n].t[j] * steps[n].bh_f_occ[j] * dc / steps[n].bh_f_occ[j]; //The number density (Mpc^-3)
                                                                                                                        //of AGN with this luminosity
                                                                                                                        //and this Mbh in this Mh bin.
        steps[n].lum_func_full[i*LBOL_BINS+k] += nd_l;
      }
    }

    // Count the total BH mass function here.
    steps[n].bhmf[i] = tnd * mbh_bpdex;

    // the distribution of luminosity distributions is simply the luminosity function
    // divided by the total number of black holes of this mass. Note that we already
     // included the duty cycle when calculating the lum_func_full, so it also lies
    // in this luminosity distribution.
    // Note that since we renormalize the luminosity distribution for each Mbh bin
    // with the total number density of BHs in this Mbh bin, we need to put the Mbh
    // bin width back into the probability of Mbh given M* when calculating qpdf at
    // different stellar masses.
    for (k=0; k<LBOL_BINS; k++) steps[n].lum_dist_full[i*LBOL_BINS+k] = steps[n].lum_func_full[i*LBOL_BINS+k] / tnd;

  }
}




// The new function to calculate the kinetic Eddington ratio distributions.
// This implementation is outdated since it still uses the scatter to smooth
// the kinetic Eddington ratio distributions.
void calc_bh_acc_rate_distribution_kinetic(int n, struct smf_fit *fit) 
{
  // No matter whether we need to calculate the kinetic Eddington ratio distribution,
  // the distributions must be initialized with zeros in the first place.
  memset(steps[n].bher_dist_kin, 0, sizeof(float)*BHER_BINS*M_BINS);

  if (!nonlinear_luminosity) return; //The linear luminosity ansatz implies that all the
                                      //gravitational energy from accretion gets converted
                                      //into radiative energy, so no kinetic energy.

  int64_t i, j, k;
  

  double bh_eta_crit = steps[n].smhm.bh_eta_crit;
  
  // The total scatter in BH mass at fixed ***halo mass*** is the quadratic 
  // sum of the scatter in the stellar mass--halo mass relation (scaled by the slope
  // of the BH mass--bulge mass relation), and that around the bulge mass--BH 
  // mass relation.
  double sm_scatter = steps[n].smhm.scatter*steps[n].smhm.bh_gamma;
  double scatter = sqrt(sm_scatter*sm_scatter + steps[n].smhm.bh_scatter*steps[n].smhm.bh_scatter);


  // Since this calculation only occurs when the luminosity is nonlinear, there is no
  // need to differentiate between linear and non-linear cases.
  double (*bher_prob)(double, double, double, double, double, double) = &_prob_of_ledd_kinetic;



  if (scatter <= 0) 
  {
    for (i=0; i<M_BINS; i++) 
    {
      for (j=0; j<BHER_BINS; j++) 
      {
        double bher = BHER_MIN+j*BHER_INV_BPDEX;
        steps[n].bher_dist_kin[i*BHER_BINS+j] = bher_prob(bher+steps[n].bh_eta[i], steps[n].bh_eta[i], steps[n].smhm.bh_alpha, steps[n].smhm.bh_delta, 1.0, bh_eta_crit)*steps[n].smhm.bh_prob_norm;
        steps[n].ledd_min[i] = BHER_MIN;
        steps[n].ledd_max[i] = BHER_MAX;
        steps[n].ledd_bpdex[i] = BHER_BPDEX;
      }
    }
  } 
  else 
  {
    const int64_t cmin = -BHER_BINS-6*BHER_BPDEX;
    const int64_t cmax = 6*BHER_BPDEX + BHER_BINS;
    float *ecache = check_realloc(NULL, sizeof(float)*(cmax-cmin+1), "Allocating ecache");

    for (i=0; i<M_BINS; i++) 
    {
      float *bher_dist = steps[n].bher_dist_kin+i*BHER_BINS;


      // Since the conversion between the ***kinetic*** and the ***total*** Eddington ratios is more complicated,
      // Here we set up the lower and upper limits in the ***kinetic*** Eddington ratio, instead of the total one.
      double ledd_min = BHER_MIN;
      double ledd_eff_min = BHER_EFF_MIN;
      double ledd_max = BHER_MAX;


      double bpdex = ((double)BHER_BINS-1)/(ledd_max - ledd_min);
      double inv_bpdex = 1.0/bpdex;
      int64_t kmax = 6.7*scatter*bpdex;
      int64_t kmin = -kmax;
      if (kmax - kmin > cmax - cmin + 1)
      {
        for (k=0; k<BHER_BINS; k++) bher_dist[k] = 0;
        free(ecache);
        return;
      }
      double inv_scatter = 1.0/scatter;

      for (k=kmin; k<kmax; k++) 
      {
        double db = inv_scatter*k*inv_bpdex;
        ecache[k-kmin] = exp(-0.5*db*db);
      }

      int64_t jmin = (ledd_eff_min-ledd_min)*bpdex;
      int64_t jmax = (ledd_max-ledd_min)*bpdex;

      for (j=jmin; j<jmax; j++) 
      {
        float bher = ledd_min+j*inv_bpdex;
        int64_t emin = j+kmin;
        int64_t emax = j+kmax;
        float *ec = ecache - emin;

        float prob = bher_prob(bher+steps[n].bh_eta[i], steps[n].bh_eta[i], steps[n].smhm.bh_alpha, steps[n].smhm.bh_delta, 1.0, bh_eta_crit);
        if (prob < 1e-15) 
        {
          continue;
        }


  
        for (k=emin; k<emax; k++) 
        {
          if (k < 0 || k > BHER_BINS - 1) continue;
          float weight = ec[k]; //exp(-0.5*db*db);
          bher_dist[k] += weight*prob;
        } 
      }

      //Normalize;
      double total = 0;
      for (k=0; k<BHER_BINS; k++) total += bher_dist[k];
      total *= inv_bpdex;
      if (!(total>0)) 
      {
      }
      if (total > 0)
  total = 1.0/total;
      for (k=0; k<BHER_BINS; k++) bher_dist[k] *= total;
    }
    free(ecache);
  }
}

// Calculate average radiative Eddington ratios for all halo
// mass bins in a certain snapshot.
void calc_avg_eta_rad(int n)
{
  int64_t i, j;
  for (i = 0; i < M_BINS; i++)
  {
    double eta0 = steps[n].bh_eta[i];
    double ledd_min = steps[n].ledd_min[i];
    double ledd_bpdex = steps[n].ledd_bpdex[i];
    double ledd_inv_bpdex = 1 / ledd_bpdex;

    double tot = 0;

    // mass- and redshift-dependent AGN duty cycle.
    double dc = steps[n].bh_duty[i];
    if (dc < 1e-4) dc = 1e-4; 

    // Count the contribution from all radiative Eddington ratio bins.
    for (j = 0; j < BHER_BINS; j++)
    {
      double prob = steps[n].bher_dist[i*BHER_BINS + j];
      if (isfinite(prob))
      {  
          tot += prob * exp10(eta0 + ledd_min + (j + 0.5) * ledd_inv_bpdex) * ledd_inv_bpdex * dc;
      }
    }
    // If we have a well-behaved Eddington ratio distribution, then
    // we adopt the calculation result.
    if (isfinite(tot)) steps[n].bh_eta_rad_avg[i] = tot;
  }
}

// Make the GSL interpolation objects for the calculation of 
// galaxy stellar mass functions (SMFs) and average specific
// star formation rates.
void calc_smf_and_ssfr(int n, struct smf_fit *fit) 
{
  int64_t i, j;
  char buffer[1024];
  double m, m2, sm, sm_max=0, sm_min=1000;
  // bin_start and bin_end are the bin # where we should do interpolations.
  // bin_peak is the bin # where the stellar mass--halo mass (SMHM) relation shows
  // a turnover. This is allowed to happen once because the Universe Machine (UM)
  // produces such turnovers at high-z, and we try to have this flexibility
  // to keep consistency with UM. 
  int64_t bin_start=-1, bin_end=-1, bin_peak = 0;
  gsl_interp_accel ga = {0};
  steps[n].smhm.sm_max = 0;
  // count_falling is the flag indicating whether we have already encountered
  // a turnover in the SMHM relation before. 
  int64_t count_falling = 0; //
  steps[n].alloc2smf = 0; //Flag that indicates whether we have a 
                          //2-segment stellar mass--halo mass relation.
  if (INVALID(*fit)) 
  {
    // printf("invalid!\n");
    return;
  }

  for (i=0; i<M_BINS; i++) 
  {
    // bin_start/bin_end is the first/last bin where stellar mass is positive.
    if (steps[n].log_sm[i]>0) 
    {
      bin_end = i;
      if (bin_start < 0) bin_start = i;
    }

    if (count_falling && i >= bin_peak 
        && steps[n].log_sm[i] > 0
        && steps[n].log_sm[i] >= steps[n].log_sm[i-1]
        && (i >= steps[n].smhm.bin_real || steps[n].log_sm[i] - steps[n].log_sm[i-1] < 0.01))
      // If it rises up again and is not really important (i > bin_real),
      // do an interpolation to keep it down.
    {
      steps[n].log_sm[i] = 2 * steps[n].log_sm[i - 1] - steps[n].log_sm[i - 2];
    }
    
     
    // Update the maximum and minium stellar masses for the snapshot, as needed.
    if (steps[n].log_sm[i]>steps[n].smhm.sm_max) 
    {
      sm_max = steps[n].smhm.sm_max = steps[n].log_sm[i];
      // bin_max = i;
    }
    if (steps[n].log_sm[i]>0 && steps[n].log_sm[i]<sm_min) 
    {
      sm_min = steps[n].smhm.sm_min = steps[n].log_sm[i];
    }

    // If we see a falling in the SMHM relation, and have not had this situation
    // before (count_falling == 0), we set count_falling = 1 and record the bin #
    // in bin_peak.
    if (i && (steps[n].log_sm[i] <= steps[n].log_sm[i-1]) && steps[n].log_sm[i] > 0 && (!count_falling)) 
    {
          count_falling = 1;
          bin_peak = i - 1;
    }

  }

  // If all the stellar masses are zero, there are too few positive stellar masses, or 
  // the SMHM relation starts to fall in the very beginning, we label the model
  // as invalid.
  if (bin_start < 0 || bin_end - bin_start + 1 < 10 || (count_falling && bin_peak <= bin_start + 1))
  {
    //sprintf(buffer, "All the stellar masses are zero.\n");
    //fprintf(stderr, "All the stellar masses are zero at z=%f.\n", 1.0 / steps[n].scale - 1);
    if (n > 0) 
    {
      //fprintf(stderr, "All the stellar masses are zero at z=%f.\n", 1.0 / steps[n].scale - 1);
      INVALIDATE(fit, buffer);
    }
    return;
  }


  // The helper arrays to do interpolations.
  double *log_sm_tmp = NULL;
  double  *hm_tmp = NULL;
  double  *sfr_tmp = NULL;

  // If the SMHM relation is monotonically rising, the interpolation
  // is very simple. We only need to do a 1-segment interpolation.
  if (!count_falling)
  {
    log_sm_tmp = malloc(sizeof(double)*(bin_end - bin_start + 1));
    hm_tmp = malloc(sizeof(double)*(bin_end - bin_start + 1));
    sfr_tmp = malloc(sizeof(double)*(bin_end - bin_start + 1));
    for (j=bin_start; j <= bin_end; j++)
    {
      log_sm_tmp[j - bin_start] = steps[n].log_sm[j];
      hm_tmp[j - bin_start] = steps[n].med_hm_at_a[j];
      sfr_tmp[j - bin_start] = steps[n].sfr[j];
    }
    
    steps[n].spline = gsl_spline_alloc(gsl_interp_linear, bin_end - bin_start + 1);
    steps[n].spline_sfr = gsl_spline_alloc(gsl_interp_linear, bin_end - bin_start + 1);
    steps[n].flag_alloc = 1;
    gsl_spline_init(steps[n].spline, log_sm_tmp, hm_tmp, bin_end - bin_start + 1);
    gsl_spline_init(steps[n].spline_sfr, log_sm_tmp, sfr_tmp, bin_end - bin_start + 1);
    free(log_sm_tmp); free(hm_tmp); free(sfr_tmp);
  }

  // If we see a turnover in the SMHM relation, segmented interpolations have to 
  // be done, because GSL interpolations require rigorously increasing independent
  // variables, which are, in this case, stellar masses.
  else
  {
    
    // For the first segment, we do interpolations as in the 1-segment 
    // interpolation case.
    log_sm_tmp = malloc(sizeof(double)*(bin_peak - bin_start + 1));
    hm_tmp = malloc(sizeof(double)*(bin_peak - bin_start + 1));
    sfr_tmp = malloc(sizeof(double)*(bin_peak - bin_start + 1));
    for (j=bin_start; j <= bin_peak; j++)
    {
      log_sm_tmp[j - bin_start] = steps[n].log_sm[j];
      hm_tmp[j - bin_start] = steps[n].med_hm_at_a[j];
      sfr_tmp[j - bin_start] = steps[n].sfr[j];
    }

    steps[n].spline = gsl_spline_alloc(gsl_interp_linear, bin_peak - bin_start + 1);
    steps[n].spline_sfr = gsl_spline_alloc(gsl_interp_linear, bin_peak - bin_start + 1);
    
    // Invalidate the model if GSL fails to make the interpolation.
    if ((!steps[n].spline) || (!steps[n].spline_sfr))
    {
      //fprintf(stderr, "Too few data points to do akima interpolation for the SMHM, segment 1. scale=%f\n", steps[n].scale);
      INVALIDATE(fit, buffer);
      return;
    }

    // For the second segment, things are slightly trickier because we need to
    // reverse the order of stellar masses (and also halo masses) on this segment, 
    // to ensure that the stellar masses are increasing.
    gsl_spline_init(steps[n].spline, log_sm_tmp, hm_tmp, bin_peak - bin_start + 1);
    gsl_spline_init(steps[n].spline_sfr, log_sm_tmp, sfr_tmp, bin_peak - bin_start + 1);
    steps[n].flag_alloc = 1;
    free(log_sm_tmp); free(hm_tmp); free(sfr_tmp);

    int n_seg = bin_end - bin_peak + 1;

    log_sm_tmp = malloc(sizeof(double)*n_seg);
    hm_tmp = malloc(sizeof(double)*n_seg);
    sfr_tmp = malloc(sizeof(double)*n_seg);

    for (j=bin_peak; j <= bin_end; j++)
    {
      log_sm_tmp[j - bin_peak] = steps[n].log_sm[bin_end - (j - bin_peak)];
      hm_tmp[j - bin_peak] = steps[n].med_hm_at_a[bin_end - (j - bin_peak)];
      sfr_tmp[j - bin_peak] = steps[n].sfr[bin_end - (j - bin_peak)];
    }

    steps[n].spline2 = gsl_spline_alloc(gsl_interp_linear, n_seg);
    steps[n].spline_sfr2 = gsl_spline_alloc(gsl_interp_linear, n_seg);
    int err_spline_init1 = gsl_spline_init(steps[n].spline2, log_sm_tmp, hm_tmp, n_seg);
    int err_spline_init2 = gsl_spline_init(steps[n].spline_sfr2, log_sm_tmp, sfr_tmp, n_seg);
    free(log_sm_tmp); free(hm_tmp); free(sfr_tmp);

    // Again, invalidate the model in case interpolation creation fails.
    if ((err_spline_init1) || (err_spline_init2))
    {
      //fprintf(stderr, "More than 1 turning point in the SMHM.\n");
      if (1.0 / steps[n].scale - 1 < 5) 
      {
        INVALIDATE(fit, buffer);
        //fprintf(stderr, "More than 1 turning point in the SMHM. at the %d-th snapshot.\n", n);
      }

      gsl_spline_free(steps[n].spline2); gsl_spline_free(steps[n].spline_sfr2);
      steps[n].flag_alloc = 1;
      return;
    }
    // In this case, we should set these flags as 1 to indicate that
    // we have allocated and made 2-segment interpolations.
    steps[n].flag_alloc = 1;
    steps[n].alloc2smf = 1;
  }

  // Calculate the intrinsic stellar mass function and average
  // star formation rate (as a function of stellar mass) based
  // on the interpolations built above.
  i = M_BINS-1;
  for (j=SM_BINS-1; j>=0; j--) 
  {
    // Skip too big/small stellar masses
    sm = SM_MIN + (double)j*SM_INV_BPDEX;
    if (sm > sm_max || sm < sm_min) 
    {
      steps[n].sfr_sm[j+SM_EXTRA] = 0;
      steps[n].sfr_sm_sf[j+SM_EXTRA] = 0;
      steps[n].sfr_sm_q[j+SM_EXTRA] = 0;
      steps[n].smf[j+SM_EXTRA] = 0;
    }
    else 
    {
      // Calculate the corresponding halo mass from the interpolation
      gsl_interp_accel_reset(&ga);
      int err = gsl_spline_eval_e(steps[n].spline, sm, &ga, &m);
      if (err || !isfinite(m))
      {
        INVALIDATE(fit, buffer);
        return;
      }

      // find out the sfrac to weight the smf_sf and smf_q.
      double f = (m - M_MIN) * BPDEX;
      int64_t b;
      if (f >= M_BINS - 1)
      {
        b = M_BINS - 2;
        f = 1;
      }
      else
      {
        b = f;
        f -= b;
      }

      double sfrac_tmp = steps[n].sfrac[b] + f * (steps[n].sfrac[b+1] - steps[n].sfrac[b]);
      steps[n].sfrac_sm[j+SM_EXTRA] = sfrac_tmp;

      // The intrinsic SMF is calculated by converting the halo mass function with the Jacobian.
      steps[n].smf[j+SM_EXTRA] = calc_smf_at_m(m,sm,n, 1, fit,&ga);
      gsl_interp_accel_reset(&ga);

      // And also calculate the average SFR.
      err = gsl_spline_eval_e(steps[n].spline_sfr, sm, &ga, &(steps[n].sfr_sm[j + SM_EXTRA]));
      
      // We need to calculate the average SFR for both quenched and star-forming galaxies.
      // Things are easier for quenched galaxies, because they (are assumed to) have
      // a fixed average specific SFR, SSFR_AVG_Q.
      steps[n].sfr_sm_q[j + SM_EXTRA] = exp10(sm) * SSFR_AVG_Q;
      // For the average SFR for star-forming galaxies, we simply need to subtract
      // quenched populations' contribution from the total value.
      steps[n].sfr_sm_sf[j + SM_EXTRA] = (steps[n].sfr_sm[j + SM_EXTRA] - (1 - sfrac_tmp) * steps[n].sfr_sm_q[j + SM_EXTRA]) / sfrac_tmp;
      // Don't forget to multiply these with stellar mass functions, as this form
      // will be easier to use in the comparisons with observations (see observations.c).
      steps[n].sfr_sm[j + SM_EXTRA] *= steps[n].smf[j+SM_EXTRA];
      steps[n].sfr_sm_q[j + SM_EXTRA] *= steps[n].smf[j+SM_EXTRA];
      steps[n].sfr_sm_sf[j + SM_EXTRA] *= steps[n].smf[j+SM_EXTRA];

      // Invalidate the model if the interpolation goes wrong.
      if (err || !isfinite(steps[n].sfr_sm[j + SM_EXTRA]))
      {
        //sprintf(buffer, "Error in GSL spline interpolation #3.\n");
        //fprintf(stderr, "Error in GSL spline interpolation #3. sfr=%e\n", steps[n].sfr_sm[j+SM_EXTRA]);
        INVALIDATE(fit, buffer);
        return;
      }

      // If we made 2-segment interpolations, simply repeat the process above for
      // the second segment.
      if (steps[n].alloc2smf)
      {
        gsl_interp_accel_reset(&ga);
        int err2 = gsl_spline_eval_e(steps[n].spline2, sm, &ga, &m2);
        //fprintf(stderr, "%d %f %f\n", n, sm, m2);
        if (!err2 && isfinite(m2))
        {
          double f2 = (m2 - M_MIN) * BPDEX;
          int64_t b2;

          if (f2 >= M_BINS - 1)
          {
            b2 = M_BINS - 2;
            f2 = 1;
          }
          else
          {
            b2 = f2;
            f2 -= b2;
          }

          if (b2 < 0) continue;

          // Given the second branch of the SMHM, we need to modify the sfrac_sm value
          // s.t. it is the weighted average of the two HM bins.
          steps[n].sfrac_sm[j+SM_EXTRA] *= steps[n].t[b] + f * (steps[n].t[b+1] - steps[n].t[b]);
          steps[n].sfrac_sm[j+SM_EXTRA] += (steps[n].sfrac[b2] + f2 * (steps[n].sfrac[b2+1] - steps[n].sfrac[b2])) *
                                            (steps[n].t[b2] + f2 * (steps[n].t[b2+1] - steps[n].t[b2]));
          steps[n].sfrac_sm[j+SM_EXTRA] /= (steps[n].t[b] + f * (steps[n].t[b+1] - steps[n].t[b]) + 
                                            steps[n].t[b2] + f2 * (steps[n].t[b2+1] - steps[n].t[b2]));
          sfrac_tmp = steps[n].sfrac[b2] + f2 * (steps[n].sfrac[b2+1] - steps[n].sfrac[b2]);
          double smf2 = calc_smf_at_m(m2,sm,n,2,fit,&ga);
          steps[n].smf[j+SM_EXTRA] += smf2;
          gsl_interp_accel_reset(&ga);
          double sfr_sm_tmp;
          double err_sfr = gsl_spline_eval_e(steps[n].spline_sfr2, sm, &ga, &(sfr_sm_tmp));
          if (!err_sfr && isfinite(sfr_sm_tmp))
          {
            steps[n].sfr_sm[j + SM_EXTRA] += sfr_sm_tmp * smf2;
            steps[n].sfr_sm_q[j + SM_EXTRA] += exp10(sm) * SSFR_AVG_Q * smf2;
            steps[n].sfr_sm_sf[j + SM_EXTRA] += (sfr_sm_tmp - (1 - sfrac_tmp) * exp10(sm) * SSFR_AVG_Q) / sfrac_tmp * smf2;
          }
        }
      }
    }
  } 

  

  //Check status of SMF bins. If the intrinsic SMFs from interpolations do
  // vary smoothly, we will directly use these pre-calculated values when
  // comparing with observations (see _interp_from_sm() in observations.c).
  steps[n].smf_ok[SM_BINS+SM_EXTRA-1] = 0;
  for (j=0; j<SM_BINS-1; j++) 
  {
    sm = SM_MIN + (double)j*SM_INV_BPDEX;
    double avg = 0.5*(steps[n].smf[j+SM_EXTRA-1]+steps[n].smf[j+SM_EXTRA+1]);
    if (fabs(avg-steps[n].smf[j+SM_EXTRA]) > steps[n].smf[j+SM_EXTRA]*5e-4) 
    {
      steps[n].smf_ok[j+SM_EXTRA] = 0;
    } 
    else 
    {
      steps[n].smf_ok[j+SM_EXTRA] = 1;
    }
  }
}

// Calculate the observed UV luminosity functions.
// The logic is very similar to calc_smf_and_ssfr().
// The only major difference is that the UV magnitude
// (M_UV) generally falls with increasing halo mass. So when
// building GSL interpolation objects, we by default
// need to reverse the order of M_UV for the first segment,
// and keep the second segment as it is, if needed.
void calc_uvlf(int n, struct smf_fit *fit) 
{
  int64_t i, j;
  char buffer[1024];
  double m, m2, uv, uv_max=-1000, uv_min=1000;
  int64_t bin_start=-1, bin_end=-1, bin_peak=0;
  gsl_interp_accel ga = {0};
  steps[n].smhm.uv_max = -1000;
  int64_t count_falling = 0;
  steps[n].alloc2uvlf = 0;
  if (INVALID(*fit)) 
  {
    return;
  }

  for (i=0; i<M_BINS; i++) 
  {

    if (steps[n].obs_uv[i] > -1000 && steps[n].obs_uv[i] < 1000) 
    {
      bin_end = i;
      if (bin_start < 0) bin_start = i;
    }

    if (count_falling && i >= bin_peak 
        && steps[n].obs_uv[i] < 10000
        && steps[n].obs_uv[i] <= steps[n].obs_uv[i-1]
        && (i >= steps[n].smhm.bin_real || steps[n].obs_uv[i-1] - steps[n].obs_uv[i] < 0.01))
      // If it rises up again and is not really important (i > bin_real),
      // do an interpolation to keep it down.
    {
      steps[n].obs_uv[i] = 2 * steps[n].obs_uv[i - 1] - steps[n].obs_uv[i - 2];
    }
    
     
    // INPLEMENTAION 2: bin_end - bin_min + 1 < M_BINS

    if (steps[n].obs_uv[i]>steps[n].smhm.uv_max && steps[n].obs_uv[i] < 1000) {
      uv_max = steps[n].smhm.uv_max = steps[n].obs_uv[i];
    }

    if (steps[n].obs_uv[i]<steps[n].smhm.uv_min) {
      uv_min = steps[n].smhm.uv_min = steps[n].obs_uv[i];
    }

    if (i && (steps[n].obs_uv[i] >= steps[n].obs_uv[i-1])
        && steps[n].obs_uv[i]<1000 && (!count_falling)) 
    {

        count_falling = 1;
        bin_peak = i - 1;
      
    }
  }

  if (bin_start < 0 || bin_end - bin_start + 1 < 10 || (count_falling && bin_peak <= bin_start + 1))
  {
    //fprintf(stderr, "All the UV magnitudes are zero at z=%f (%d-th snapshot). bin_min=%d, bin_max=%d\n", 1.0 / steps[n].scale - 1, n, bin_min, bin_max);
    if (1.0 / steps[n].scale - 1 >= 7.5) 
    {
      //fprintf(stderr, "All the UV magnitudes are zero at z=%f (%d-th snapshot). bin_end=%d, bin_start=%d\n", 1.0 / steps[n].scale - 1, n, bin_end, bin_start);
      INVALIDATE(fit, buffer);
    }
    return;
  }

  double *obs_uv_tmp = NULL;
  double  *hm_tmp = NULL;

  if (!count_falling)
  {
    obs_uv_tmp = malloc(sizeof(double)*(bin_end - bin_start + 1));
    hm_tmp = malloc(sizeof(double)*(bin_end - bin_start + 1));

    for (j=bin_start; j <= bin_end; j++)
    {
      obs_uv_tmp[j - bin_start] = steps[n].obs_uv[bin_end - (j - bin_start)];
      hm_tmp[j - bin_start] = steps[n].med_hm_at_a[bin_end - (j - bin_start)];
    }

    // The spline does not need to be modified. What should be changed is that the dn/dlogm
    // should be weighted by sfrac and 1 - sfrac for smf_sf and smf_q, respectively.
    //steps[n].spline_uv = gsl_spline_alloc(gsl_interp_cspline, bin_end - bin_start + 1);
    steps[n].spline_uv = gsl_spline_alloc(gsl_interp_linear, bin_end - bin_start + 1);

    steps[n].flag_alloc = 2;
    gsl_spline_init(steps[n].spline_uv, obs_uv_tmp, hm_tmp, bin_end - bin_start + 1);
    free(obs_uv_tmp); free(hm_tmp);
  }

  else
  {
    obs_uv_tmp = malloc(sizeof(double)*(bin_peak - bin_start + 1));
    hm_tmp = malloc(sizeof(double)*(bin_peak - bin_start + 1));

    for (j=bin_start; j <= bin_peak; j++)
    {
      obs_uv_tmp[j - bin_start] = steps[n].obs_uv[bin_peak - (j - bin_start)];
      hm_tmp[j - bin_start] = steps[n].med_hm_at_a[bin_peak - (j - bin_start)];
    }

    
    // The spline does not need to be modified. What should be changed is that the dn/dlogm
    // should be weighted by sfrac and 1 - sfrac for smf_sf and smf_q, respectively.
    //steps[n].spline_uv = gsl_spline_alloc(gsl_interp_cspline, bin_peak - bin_start + 1);
    steps[n].spline_uv = gsl_spline_alloc(gsl_interp_linear, bin_peak - bin_start + 1);
    
    if (!steps[n].spline_uv)
    {
      //fprintf(stderr, "Too few data points to do akima interpolation for the UVHM, segment 1. scale=%f\n", steps[n].scale);
      INVALIDATE(fit, buffer);
      return;
    }

    gsl_spline_init(steps[n].spline_uv, obs_uv_tmp, hm_tmp, bin_peak - bin_start + 1);

    free(obs_uv_tmp); free(hm_tmp);

    // Here we need to account for the fact that the segment of UVHM relation after
    // the turning over point may be fewer than 2 points. In this case
    // cubic spline interpolation won't work. And my solution is to do a linear
    // extrapolation out to 10^16.3 Msun to get the third point.
    // int n_seg = bin_end - bin_peak + 1 > 4 ? bin_end - bin_peak + 1 : 5;
    int n_seg = bin_end - bin_peak + 1;

    obs_uv_tmp = malloc(sizeof(double)*n_seg);
    hm_tmp = malloc(sizeof(double)*n_seg);

    for (j=bin_peak; j <= bin_end; j++)
    {
      obs_uv_tmp[j - bin_peak] = steps[n].obs_uv[j];
      hm_tmp[j - bin_peak] = steps[n].med_hm_at_a[j];
    }

    steps[n].spline_uv2 = gsl_spline_alloc(gsl_interp_linear, n_seg);
    int err_spline_init = gsl_spline_init(steps[n].spline_uv2, obs_uv_tmp, hm_tmp, n_seg);
        
    free(obs_uv_tmp); free(hm_tmp);
    if (err_spline_init)
    {
      //fprintf(stderr, "More than 1 turning point in the UVHM at %d-th snapshot.\n", n);
      if (1.0/steps[n].scale - 1 > 7.5) 
      {
        //fprintf(stderr, "More than 1 turning point in the UVHM at %d-th snapshot.\n", n);
        INVALIDATE(fit, buffer);
      }
      gsl_spline_free(steps[n].spline_uv2);
      steps[n].flag_alloc = 2;
      return;
    }
    steps[n].alloc2uvlf = 1;
    steps[n].flag_alloc = 2;   
  }


  
  i = M_BINS-1;
  for (j=UV_BINS-1; j>=0; j--) 
  {
    uv = UV_MIN + (double)j*UV_INV_BPMAG;
    if (uv > uv_max || uv < uv_min) 
    {
      steps[n].uvlf[j+UV_EXTRA] = 0;
      steps[n].std_uvlf[j+UV_EXTRA] = 0;
    }
    else 
    {
      gsl_interp_accel_reset(&ga);
      int err = gsl_spline_eval_e(steps[n].spline_uv, uv, &ga, &m);
      // fprintf(stdout, "n=%d, z=%f, uv=%f, interpolated m=%f, err=%d\n", n, 1.0 / steps[n].scale - 1, uv, m, err);
      if (err || !isfinite(m))
      //if (err || m < M_MIN - 1.0 || m > M_MAX + 1.0)
      {
        // sprintf(buffer, "Error in GSL spline interpolation #2.\n");
        //fprintf(stderr, "Error in GSL spline interpolation #2. n=%d, uv=%f, m=%f\n", n, uv, m);
        //for (int ii=0; ii<M_BINS; ii++) fprintf(stderr, "obs_uv[%d]=%f\n", ii, steps[n].obs_uv[ii]);
        INVALIDATE(fit, buffer);
        // free(obs_uv_tmp); free(hm_tmp);
        return;
      }
      // find out the sfrac to weight the smf_sf and smf_q.
      double f = (m - M_MIN) * BPDEX;
      int64_t b;
      if (f >= M_BINS - 1)
      {
        b = M_BINS - 2;
        f = 1;
      }
      else
      {
        b = f;
        f -= b;
      }
      steps[n].std_uvlf[j+SM_EXTRA] = steps[n].std_uv[b] + f * (steps[n].std_uv[b+1] - steps[n].std_uv[b]);
      steps[n].uvlf[j+SM_EXTRA] = calc_uvlf_at_m(m,uv,n,1,fit,&ga);
      if (steps[n].alloc2uvlf)
      {
        gsl_interp_accel_reset(&ga);
        int err2 = gsl_spline_eval_e(steps[n].spline_uv2, uv, &ga, &m2);
        if (!err2 && isfinite(m2))
        {
          steps[n].uvlf[j+SM_EXTRA] += calc_uvlf_at_m(m2,uv,n,2,fit,&ga);
        }
      }
    }
  }

  //Check status of UVLF bins
  steps[n].uvlf_ok[UV_BINS+UV_EXTRA-1] = 0;
  for (j=0; j<UV_BINS-1; j++) 
  {
    uv = UV_MIN + (double)j*UV_INV_BPMAG;
    double avg = 0.5*(steps[n].uvlf[j+UV_EXTRA-1]+steps[n].uvlf[j+UV_EXTRA+1]);
    if (fabs(avg-steps[n].uvlf[j+UV_EXTRA]) > steps[n].uvlf[j+UV_EXTRA]*5e-4) 
    {
      steps[n].uvlf_ok[j+UV_EXTRA] = 0;
    } 
    else 
    {
      steps[n].uvlf_ok[j+UV_EXTRA] = 1;
    }
  }
}

// Create fake star formation rate histories
// for halos that only emerge after some time
// in the simulation.
void create_fake_sfr_hist(int n, int i) 
{
  double frac_lost = 0, weight = 0;
  double sm, total_time, time, prev_scale;
  int k;

  for (k=0; k<=n; k++) 
  {
    frac_lost += steps[n].smloss[k]*steps[k].dt;
    weight += steps[k].dt;
  }
  frac_lost /= weight; //Average frac lost for constant sfr

  // Apply the stellar mass loss.
  sm = steps[n].sm[i] / frac_lost;

  // allocate the star formation uniformly over time.
  total_time = scale_to_years(steps[n].scale)-scale_to_years(0);
  for (k=0; k<=n; k++) 
  {
    prev_scale = (n) ? steps[n-1].scale : 0;
    time = scale_to_years(steps[n].scale) - scale_to_years(prev_scale);
    steps[n].sm_hist[i*num_outputs + k] = sm*time/total_time;
  }
  // And the new stellar mass is naturally the latest star formation history.
  steps[n].new_sm[i] = steps[n].sm_hist[i*num_outputs + n];

}

// Deprecated.
double icl_fract(int64_t i, int64_t j, int64_t n, double ICL_RATIO) 
{
  //double icl_frac;
  /*if (j>i) icl_frac = 1;
    else icl_frac = 1.0 - (2.0)/(1.0+doexp10((i-j+0.25)*INV_BPDEX*ICL_RATIO));*/
  /*  double exponent = steps[n].icl_exp; //0.3*pow(steps[n].scale, -0.7);
  if (j>i) icl_frac = 1.0-0.16;
  else icl_frac = 1.0 - 0.16*doexp10(-1.0*(i-j+0.25)*INV_BPDEX*exponent);
  if (n> 0 && steps[n-1].n[i]<0.5*steps[n].n[i]) icl_frac = 1;
  */
  return 1.0;
}

// Calculate the star formation histories for a certain snapshot.
void calc_sm_hist(int n, struct smf_fit *fit) 
{
  int64_t i, j, k, bin;
  double t;
  // sm_hist_i/j is the stellar mass histories of the i/j-th halo mass bin.
  // sm_inc is the amount of stellar mass coming from infalling galaxies,
  // and sm_mmp is the stellar mass inherited from the most-massive progenitors.
  double *sm_hist_i, *sm_hist_j, sm_inc, sm_mmp;

  // The intracluster light (ICL) mass histories of the i/j-th halo mass bin.
  double *icl_hist_i, *icl_hist_j;
  // The incoming stellar mass histories of the i/j-th halo mass bin.
  double *sm_inc_hist_i;
  double ejec_frac;
  steps[n].smhm = smhm_at_z((1.0/steps[n].scale)-1.0, *fit);
  if (n && steps[n].smhm.rho_bh > steps[n-1].smhm.rho_bh)
  {
	fprintf(stderr, "Increasing rho_bh.\n");
	INVALIDATE(fit, "Increasing rho_bh");
  }
  // Invalidate the model if the fraction that the incoming satellite galaxies
  // get merged into central galaxies is too big or too small.
  if (steps[n].smhm.icl_frac < 1e-20 || steps[n].smhm.icl_frac > 1)
  {
        //fprintf(stderr, "Bad ICL fraction: %e at z=%f\n", steps[n].smhm.icl_frac, 1/steps[n].scale-1);
        INVALIDATE(fit, "Bad ICL fraction");
  }

  // Calculate the scaling relation between SFR and UV magnitudes, based on the
  // fitting results shown in Appendix D of Zhang et al. (2021).
  // We calculate the k_uv and b_uv at the bin "bin_real", k_uv_real and b_uv_real,
  // before OMP takes effect, otherwise for some mass bins both of these 
  // values could be zero and cause errors.

  double z, a1, mh, k_uv_real, b_uv_real, std_uv_real, k_std_uv, b_std_uv, a_k, b_k, c_k, a_b, b_b, c_b;
  // double a_k_std, b_k_std, a_b_std, b_b_std;
  z = 1.0 / steps[n].scale - 1;
  a1 = steps[n].scale - 1;
  a_k = 0.15367887; b_k = -2.87608013; c_k = 9.4778494 + (-2.37851257) * a1;
  a_b = -0.34749622; b_b = 6.85270974; c_b = -50.34457998 + 1.99341162 * a1;

  mh = steps[n].med_hm_at_a[steps[n].smhm.bin_real];
  
  k_uv_real = a_k * mh * mh + b_k * mh;
  b_uv_real = a_b * mh * mh + b_b * mh;
  k_uv_real = mh < - 0.5 * b_k / a_k ? -0.25 * b_k * b_k / a_k + c_k : k_uv_real + c_k;
  b_uv_real = mh < - 0.5 * b_b / a_b ? -0.25 * b_b * b_b / a_b + c_b : b_uv_real + c_b;
  
  k_std_uv = -0.03098716 * z + 0.04258996;
  b_std_uv = 0.31872631 * z + 0.24110353;

  std_uv_real = k_std_uv * mh + b_std_uv;

  steps[n].flag_alloc = 0;

  // a200, m200, and v200 are necessary for converting between halo mass and maximum
  // rotational velocity.
  double a200 = steps[n].scale/0.3782;
  double m200 = log10(1.115e12/0.68 / (pow(a200, -0.142441) + pow(a200, -1.78959)));
  double v200 = log10(200);
  
#pragma omp for schedule(dynamic,5) private(j,k,sm_hist_i,sm_hist_j,bin,t)
  for (i=0; i<M_BINS; i++) 
  {
    double lv = v200 + (steps[n].med_hm_at_a[i]-m200)/3.0;
    steps[n].lv[i] = lv;
    // Here the SFR only includes the star-forming ones, which is not complete.
    // The reason why we should also account for the SFR in the quenched galaxies
    // is that they may have high enough merger fraction and thus even a small SFR
    // still implies a decent total growth. 
    // But given that we assume a constant ***SSFR*** for quenched galaxies, we cannot 
    // correct for that here. Instead we have to do this after calculating the 
    // old_sm.
    // This shouldn't affect anything, because the sm[i] calculated in line 962 is useful
    // only for generating fake SFH, where all the galaxies should be SF-ing.
    // steps[n].sfrac[i] = calc_sfrac_at_lv(lv, steps[n].smhm);
    steps[n].sfrac[i] = calc_sfrac_at_m(steps[n].lv[i], steps[n].smhm);
    steps[n].sfr[i] = calc_sfr_at_lv(steps[n].med_hm_at_a[i], lv, steps[n].smhm)*exp(pow(0.30*log(10), 2)/2.0);
    //if (n >= 111 && n <= 152 && i == 9) fprintf(stderr, "n=%d, i=%d, calling calc_sm_hist().\n", n, i);
    // Ignore the SFR in the first snapshot.
    if (n == 0) steps[n].sfr[i] = 0;
    
    // Note that we shouldn't calculate the BH mass here, as the real stellar mass will not
    // be calculated until calc_new_sm_and_sfr().
    steps[n].sm[i] = steps[n].sfr[i]*steps[n].dt;
    steps[n].sm_from_icl[i] = 0;
    steps[n].old_bh_mass[i] = 0;
    steps[n].bh_unmerged[i] = 0;
    if (no_z_scaling) continue;

    // steps[n].n[i] is the number density of the halos in the i-th
    // bin that are inherited from their most-massive progenitors.
    // steps[n].c[i] is the number density of the halos in the i-th
    // bin that emerge in the n-th snapshot. So if steps[n].n[i] == 0 
    // but steps[n].c[i] > 0, we need to take care of these emerging
    // halos by creating fake star formation histories for them.
    if (!steps[n].n[i]) 
    {
      if (steps[n].c[i]) 
      {
	      create_fake_sfr_hist(n, i);
      }
    }

    // If this is not the case, we should have these halos
    // inherit the stellar mass, ICL, and incoming satellite
    // galaxy mass histories from their most-massive progenitors (MMPs).
    else 
    {
      sm_hist_i = &(steps[n].sm_hist[i*num_outputs]);
      memset(sm_hist_i, 0, sizeof(double)*num_outputs);
      icl_hist_i = &(steps[n].icl_stars[i*num_outputs]);
      memset(icl_hist_i, 0, sizeof(double)*num_outputs);
      sm_inc_hist_i = &(steps[n].sm_inc_hist[i*num_outputs]);
      memset(sm_inc_hist_i, 0, sizeof(double)*num_outputs);
      if (n>0) 
      {
      	sm_inc = sm_mmp = 0;
      	for (j=0; j<M_BINS; j++) 
        {
      	  bin = j*M_BINS + i;

          // Calculate the old BH mass inherited from their MMPs.
          steps[n].old_bh_mass[i] += steps[n-1].bh_mass_avg[j]*steps[n-1].bh_f_occ[j]*steps[n].mmp[bin];
          // Unmerged (or wandering) BH masses come from: 1) unmerged BH masses
          // that have not been used up by the MMPs; 2) the BH masses from incoming
          // satellite galaxies.
      	  steps[n].bh_unmerged[i] += steps[n-1].bh_f_occ[j]*(steps[n].merged[bin]*(steps[n-1].bh_mass_avg[j] + 
                                    steps[n-1].bh_unmerged[j]) + steps[n].mmp[bin]*steps[n-1].bh_unmerged[j]);
          
          // steps[n].bh_merged[i] += steps[n].merged[bin]*(steps[n-1].bh_mass_avg[j]);
          // // The following lines are used to calculate the mass distributions of unmerged BHs for the i-th mass bin
          // // in the n-th snapshot.
          // for (k=0; k < MBH_BINS; k++)
          // {
          //   // Just take over all the distributions from the merged BHs, up to a transfer fraction that is the fraction of 
          //   // the j-th mass bin halos that got merged/become MMP
          //   if (steps[n-1].t[j] && isfinite(steps[n-1].t[j]) && isfinite(steps[n].merged[bin]) && isfinite(steps[n].mmp[bin]))
          //     steps[n].bh_unmerged_dist[i*MBH_BINS+k] += steps[n].merged[bin] / steps[n-1].t[j] * steps[n-1].bh_unmerged_dist[j*MBH_BINS+k] 
          //                                             + steps[n].mmp[bin] / steps[n-1].t[j] * steps[n-1].bh_unmerged_dist[j*MBH_BINS+k];
          // }

          // // Another part is the average BH mass from the merged BHs.
          // // Here we assume that the mass distribution ***within each BH mass bin*** is flat.
          // // double f_bh = (steps[n-1].log_bh_mass[j] - MBH_MIN) / MBH_INV_BPDEX;
          // double mbh = log10(steps[n-1].bh_mass_avg[j]);
          // cloud_in_cell(mbh, steps[n].merged[bin], steps[n].bh_unmerged_dist+i*MBH_BINS, MBH_MIN, MBH_BPDEX, MBH_BINS);
          
      	  sm_inc += steps[n].merged[bin]*steps[n-1].sm_avg[j];
      	}

      	steps[n].sm_inc[i] = sm_inc;
      	sm_inc = 0;
      	ejec_frac = 0;

        // ejec_frac is the fraction of the stellar mass from the
        // most-massive progenitors that is ejected into ICL during 
        // mergers. In the current model, we assume ejec_frac = 0.0,
        // i.e., stellar masses in the most-massive progenitors are
        // not ejected at all during mergers.
      	if (ejec_frac > 1) ejec_frac = 1;
      	if (ejec_frac < 0) ejec_frac = 0;

        // Transfer the stellar mass from the j-th halo mass bin 
        // (the most-massive progenitor) to the i-th halo mass bin.
      	for (j=0; j<M_BINS; j++) 
        {
      	  bin = j*M_BINS + i;
          // incoming_frac is the fraction of incoming satellite
          // stellar masses that is merged into galaxy mergers.
          // This is set to zero here, i.e., all the incoming
          // satellite stellar masses firstly go into the
          // intracluster light (ICL). We will take care of 
          // the merging of incoming satellite mass into cenral
          // galaxies later.
      	  double incoming_frac = 0;
          // t is the number density of most-massive progenitors (MMPs),
          // which determines the amount of stellar masses inherited
          // from MMPs.
      	  t = steps[n].mmp[bin]*(1.0-ejec_frac) + 
      	    incoming_frac*steps[n].merged[bin];
          // iclt1 and iclt2 are the contributions to ICL from MMPs
          // and satellites, respectively.
      	  double iclt1 = (ejec_frac*steps[n].mmp[bin]+
      			 (1.0-incoming_frac)*steps[n].merged[bin]);
      	  double iclt2 = (steps[n].mmp[bin]);
      	  if (!t && !iclt1 && !iclt2) continue;
      	  steps[n].sm_from_icl[i] += steps[n-1].sm_from_icl[j]*steps[n].mmp[bin];
      	  sm_hist_j = &(steps[n-1].sm_hist[j*num_outputs]);
      	  icl_hist_j = &(steps[n-1].icl_stars[j*num_outputs]);
          // transfer stellar mass, ICL mass, and incoming satellite mass histories.
      	  for (k=0; k<n; k++) 
          {
      	    sm_hist_i[k] += t*sm_hist_j[k];
      	    icl_hist_i[k] += iclt1*sm_hist_j[k];
      	    icl_hist_i[k] += iclt2*icl_hist_j[k];
            // sm_inc_hist stores the ICL that come from the stars of satellite galaxies,
            // excluding the ICL from satellites.
            sm_inc_hist_i[k] += steps[n].merged[bin] * sm_hist_j[k];
      	  }
      	}
      }
      // We want the average mass, so dividing by the halo number density is needed.
      for (k=0; k<n; k++) sm_hist_i[k] /= steps[n].t[i];
      for (k=0; k<n; k++) icl_hist_i[k] /= steps[n].t[i];
      for (k=0; k<n; k++) sm_inc_hist_i[k] /= steps[n].t[i];
      steps[n].sm_from_icl[i] /= steps[n].t[i];
      steps[n].sm_inc[i] /= steps[n].t[i];
      steps[n].old_bh_mass[i] /= steps[n].t[i];
      // steps[n].old_bh_mass_self[i] /= steps[n].t[i];
      // steps[n].old_bh_mass_next[i] /= steps[n].t[i];
      steps[n].bh_unmerged[i] /= steps[n].t[i];
      // steps[n].bh_merged[i] /= steps[n].t[i];
      // steps[n].bh_unmerged_avail[i] = steps[n].bh_unmerged[i];
      //for (k=0; k<MBH_BINS; k++) steps[n].bh_unmerged_dist /= steps[n].t[i]; // Don't forget the normalization.
    }
    // smooth_bins(steps[n].bh_unmerged_dist + i*MBH_BINS, steps[n].bh_unmerged_dist + i*MBH_BINS,
    //             combined_scatter, MBH_MIN, MBH_BPDEX, MBH_BINS,
    //                               MBH_MIN, MBH_BPDEX, MBH_BINS, 0, 0);

    // Calculate the old stellar masses by integrating over the star formation histories.
    calc_old_sm(n,i);

    // correct the SFR by accounting for the contribution from quenched galaxies.
    steps[n].sfr[i] += (1 - steps[n].sfrac[i]) * steps[n].old_sm[i] * SSFR_AVG_Q;
    if (n == 0) steps[n].sfr[i] = 0;


    // Calculate observed UV magnitudes. UV magnitudes are found to have a linear
    // relation with log10(SFR), with redshift and halo mass dependent slope (k_uv) and
    // intercept (b_uv). The log-normal scatter around this relation (std_uv) is 
    // also a function of redshift and halo mass.
    double k_uv, b_uv, std_uv;
    mh = steps[n].med_hm_at_a[i];
    k_uv = a_k * mh * mh + b_k * mh;
    b_uv = a_b * mh * mh + b_b * mh;
    k_uv = mh < - 0.5 * b_k / a_k ? -0.25 * b_k * b_k / a_k + c_k : k_uv + c_k;
    b_uv = mh < - 0.5 * b_b / a_b ? -0.25 * b_b * b_b / a_b + c_b : b_uv + c_b;

    std_uv = k_std_uv * mh + b_std_uv;

    if (i > steps[n].smhm.bin_real)
    {
      k_uv = k_uv_real;
      b_uv = b_uv_real;
      std_uv = std_uv_real;
    }

    steps[n].k_uv[i] = k_uv;
    steps[n].b_uv[i] = b_uv;

    double lgSFR = log10(steps[n].sfr[i]);
    steps[n].obs_uv[i] = k_uv * lgSFR + b_uv;
    steps[n].std_uv[i] = std_uv;
    if (steps[n].obs_uv[i] > 0) steps[n].obs_uv[i] = 10000;
    if (steps[n].std_uv[i] < 0) steps[n].std_uv[i] = 0.001; 
    if (steps[n].std_uv[i] > 1) steps[n].std_uv[i] = 1.000; 

    // Calculate the newly formed stellar masses, SFRs, BH occupation fractions,
    // new BH masses, and BH accretion and merger rates.
    
    //if (n >= 111 && n <= 152 && i == 9) fprintf(stderr, "n=%d, i=%d, before calling calc_new_sm_and_sfr().\n", n, i);
    calc_new_sm_and_sfr(n,i,fit);

    

  }


  

}

// Calculate the old stellar masses by integrating over the star formation histories.
// Note that steps[n].smloss[k] is the fraction of ***remaining*** stellar mass that
// has formed in the k-th snapshot by the n-th snapshot.
void calc_old_sm(int n, int j) 
{
  int k;
  steps[n].old_sm[j] = 0;
  for (k=0; k<n; k++)
    steps[n].old_sm[j] += steps[n].smloss[k]*steps[n].sm_hist[j*num_outputs+k];
  steps[n].sm_icl[j] = 0;
  for (k=0; k<n; k++)
    steps[n].sm_icl[j] += steps[n].smloss[k]*steps[n].icl_stars[j*num_outputs+k];
}

double prior_focc()
{
  int n, i;
  double chi2_prior = 0;
  for (n=0; n<num_outputs; n++)
  {
    for (i=0; i<M_BINS; i++)
    {
      if (steps[n].scale > 0.16 && steps[n].med_hm_at_a[i] >= 11)
      {
        
        if (steps[n].bh_f_occ[i] > steps[n].bh_f_occ_max[i])
        {
          chi2_prior += 10000 * log10(steps[n].bh_f_occ[i] / steps[n].bh_f_occ_max[i]);
        }
      }
    }
  }
  
  return chi2_prior;
}

// Calculate the number of quasars that exist between (z_low, z_high), and
// above a certain luminosity threshold, lbol, from a survey that covers 
// frac_area of the WHOLE SKY. For example, SDSS has an area of ~14000 deg^2,
// so frac_area = 14000 / 41252.
double number_high_z_low_mass_qso(double z_low, double z_high, double lbol, double frac_area)
{
  int n, i, j;
  // Find out the relevant snapshots
  double a_max = 1.0 / (1 + z_low);
  double a_min = 1.0 / (1 + z_high);
  double n_qso = 0;

  // Find out the relevant luminosity bins
  double lbol_f = (lbol - LBOL_MIN) * LBOL_BPDEX;
  int64_t lbol_b = lbol_f; lbol_f -= lbol_b;
  if (lbol_b >= LBOL_BINS - 1) {lbol_b = LBOL_BINS - 2; lbol_f = 1;}

  // go to the first snapshot within (z_low, z_high)
  for (n=0; n<num_outputs; n++) if (steps[n].scale >= a_min) break;
  for (; n<num_outputs; n++)
  {
    double n_qso_n = 0;
    // Calculate the comoving volume around this snapshot
    double z_max = 0.5 * (1.0 / steps[n-1].scale + 1.0 / steps[n].scale - 2);
    double z_min = 0.5 * (1.0 / steps[n+1].scale + 1.0 / steps[n].scale - 2);
    double V_comoving = comoving_volume(z_max) - comoving_volume(z_min);
    // traverse every MBH bin.
    for (i=0; i<MBH_BINS; i++)
    {
      // Only count BHs that are above the luminosity cut.
      double n_qso_mbh = (1 - lbol_f) * steps[n].lum_func_full[i*LBOL_BINS+lbol_b];
      for (j=lbol_b+1; j<LBOL_BINS; j++)
      {
        n_qso_mbh += steps[n].lum_func_full[i*LBOL_BINS+j];
      }
      // weight every MBH bin by the fraction that BHs in this bin get scattered
      // below 10^8 Msun due to virial estimate. This gives the number ***density***,
      // i.e., the Mpc^-3
      n_qso_n += n_qso_mbh * frac_below8[i];
    }
    // Multiply the numeber density with comoving volume.
    n_qso += n_qso_n * V_comoving;
    // Stop if the redshift is lower than the threshold.
    if (steps[n].scale > a_max) break;
  }
  // The calculated comoving volume is for the 4pi surface area, so we have to
  // multiply it by the fractional survey area relative to the whole sky.
  // Also, the luminosity functions are in the unit of dex^-1, we have to convert
  // them back by multiplying LBOL_INV_BPDEX.
  n_qso *= frac_area * LBOL_INV_BPDEX;
  return n_qso;
}

// Calculate the number ratio of low-mass (below 10^8 Msun) vs. massive (above 10^8 Msun)
// bright quasars between (z_low, z_high).
double ratio_high_z_qso(double z_low, double z_high, double lbol, double frac_area)
{
  int n, i, j;
  // Find out the relevant snapshots
  double a_max = 1.0 / (1 + z_low);
  double a_min = 1.0 / (1 + z_high);

  double n_qso_above = 0;
  double n_qso_below = 0;

  // Find out the relevant luminosity bins
  double lbol_f = (lbol - LBOL_MIN) * LBOL_BPDEX;
  int64_t lbol_b = lbol_f; lbol_f -= lbol_b;
  if (lbol_b >= LBOL_BINS - 1) {lbol_b = LBOL_BINS - 2; lbol_f = 1;}

  // go to the first snapshot within (z_low, z_high)
  for (n=0; n<num_outputs; n++) if (steps[n].scale >= a_min) break;
  for (; n<num_outputs; n++)
  {
    // double n_qso_n = 0;
    double n_qso_n_above = 0;
    double n_qso_n_below = 0;

    // Calculate the comoving volume around this snapshot
    double z_max = 0.5 * (1.0 / steps[n-1].scale + 1.0 / steps[n].scale - 2);
    double z_min = 0.5 * (1.0 / steps[n+1].scale + 1.0 / steps[n].scale - 2);
    double V_comoving = comoving_volume(z_max) - comoving_volume(z_min);
    for (i=0; i<MBH_BINS; i++)
    {
      double n_qso_mbh = (1 - lbol_f) * steps[n].lum_func_full[i*LBOL_BINS+lbol_b];
      for (j=lbol_b+1; j<LBOL_BINS; j++)
      {
        n_qso_mbh += steps[n].lum_func_full[i*LBOL_BINS+j];
      }
      // weight every MBH bin by the fraction that BHs in this bin get scattered
      // below and above 10^8 Msun due to virial estimate. This gives the number ***density***,
      // i.e., the Mpc^-3
      n_qso_n_below += n_qso_mbh * frac_below8[i];
      n_qso_n_above += n_qso_mbh * (frac_below11[i] - frac_below8[i]);
    }
    // Multiply the numeber density with comoving volume.
    n_qso_below += n_qso_n_below * V_comoving;
    n_qso_above += n_qso_n_above * V_comoving;
    // Stop if the redshift is lower than the threshold.
    if (steps[n].scale > a_max) break;
  }

  // The calculated comoving volume is for the 4pi survey area, but
  // we are calculating the number ratio, so the area does not matter.
  double ratio = n_qso_below / n_qso_above;
  return ratio;
}

// Penalize the rising SFH histories among massive halos.
double rising_sfh_penalty(void)
{
  double penalty = 0.0;
  int64_t n, i;
  int64_t imin = (RISING_SFH_M_CONSTRAINT - M_MIN) * BPDEX;
  for (n=0; n<num_outputs; n++) if (steps[n].scale > RISING_SFH_A_CONSTRAINT) break;
  int64_t nmin = n;
  for (i=imin; i<M_BINS; i++)
  {
    double sfr = steps[nmin].sfr[i] - (1 - steps[nmin].sfrac[i]) * steps[n].old_sm[i] * SSFR_AVG_Q;
    double sfr_next;
    for (n=nmin; n<num_outputs - 1; n++)
    {
      sfr_next = steps[n+1].sfr[i] - (1 - steps[n+1].sfrac[i]) * steps[n+1].old_sm[i] * SSFR_AVG_Q;
      if (sfr_next > sfr && sfr > 0) penalty += 10 * log10(sfr_next / sfr);
      sfr = sfr_next;
    }  
  }
  return penalty;
}

// Calculate the average SFR among massive halos
double recent_sfh_in_massive_halos(void) 
{
  int64_t n, j, count=0;
  double sfr = 0;
  // find the redshift below which the SFRs should be counted.
  for (n=0; n<num_outputs; n++) if (steps[n].scale > SFR_A_CONSTRAINT) break;
  for (; n<num_outputs; n++) 
  {
    double z = 1.0 / steps[n].scale - 1;
    // offset between the observed and true SFRs.
    double corr = doexp10(steps[n].smhm.mu + steps[n].smhm.kappa * exp(-0.5 * (z - 2) * (z - 2)));
    for (j=(SFR_M_CONSTRAINT-M_MIN)*BPDEX; j<M_BINS; j++) 
    {
      if (steps[n].t[j] <= 0) continue;
      sfr += steps[n].sfr[j]*corr;
      count++;
    }
  }
  sfr /= count;
  return sfr;
}

// Calculate the recent AGN kinetic power in massive halos.
double recent_kinetic_power_in_massive_halos(void) 
{
  int64_t n, j, k;
  double L_kin = 0;
  double nd_tot = 0;
  double norm = BHER_BPDEX; //From the calculation of kinetic ER distribution, we know that the normalization of the distribution
                            //is simply BHER_BPDEX.

  // Find the smallest scale factor that enters the window
  for (n=0; n<num_outputs; n++) if (steps[n].scale > KIN_POWER_A_CONSTRAINT_LOW) break;
  for (; n<num_outputs; n++) 
  {
    // If the scale factor is bigger than the upper limit, stop.
    if (steps[n].scale > KIN_POWER_A_CONSTRAINT_HIGH) break;

    // Only count the kinetic powers in massive halos.
    for (j=(KIN_POWER_M_CONSTRAINT-M_MIN)*BPDEX; j<M_BINS; j++) 
    {
      double L_tmp = 0;
      if (steps[n].t[j] <= 0) continue;

      // mass- and redshift-dependent AGN duty cycle.
      double dc = steps[n].bh_duty[j];
      if (dc < 1e-4) dc = 1e-4;

      for (k=0; k<BHER_BINS; k++) 
      {
          L_tmp += exp10(38.1 + steps[n].log_bh_mass[j] + BHER_MIN + k*BHER_INV_BPDEX + steps[n].bh_eta[j])
                                            * steps[n].bher_dist_kin[j * BHER_BINS + k];
      }
      
      L_tmp /= norm;
      L_kin += L_tmp * dc * steps[n].t[j];
      nd_tot += steps[n].t[j];
    }
  }

  L_kin /= nd_tot;
  return log10(L_kin);
}

// Calculate the average radiative Eddington ratio.
void calc_bh_eta_avg(int n) 
{
  int64_t i, k;
  for (i=0; i<M_BINS; i++) 
  {
    if (steps[n].t[i] <= 0) continue;
    double norm = 0;
    double eta_avg = 0;
    double bher_min = steps[n].ledd_min[i] + steps[n].bh_eta[i];
    double inv_bpdex = 1.0 / steps[n].ledd_bpdex[i];

    // Need to account for the fact that only part of the SMBHs are active, which is parametrized by the duty cycle.
    double dc = steps[n].bh_duty[i];
    if (dc < 1e-4) dc = 1e-4;

    for (k = 0; k < BHER_BINS; k++) 
    {
      eta_avg += exp10(bher_min + k * inv_bpdex) * steps[n].bher_dist[i * BHER_BINS + k];
      norm += steps[n].bher_dist[i * BHER_BINS + k];
    }
    eta_avg /= norm;
    eta_avg *= dc;
    steps[n].bh_eta_avg[i] = log10(eta_avg);
  }
}

// Calculate the average radiative Eddington ratio.
void calc_bh_eta_kin_avg(int n) 
{
  int64_t i, k;
  for (i=0; i<M_BINS; i++) 
  {
    if (steps[n].t[i] <= 0) continue;
    double norm = 0;
    double eta_avg = 0;
    double bher_min = BHER_MIN + steps[n].bh_eta[i];
    double inv_bpdex = 1.0 / BHER_BPDEX;

    // Need to account for the fact that only part of the SMBHs are active, which is parametrized by the duty cycle.
    double dc = steps[n].bh_duty[i];
    if (dc < 1e-4) dc = 1e-4;

    for (k = 0; k < BHER_BINS; k++) 
    {
      eta_avg += exp10(bher_min + k * inv_bpdex) * steps[n].bher_dist_kin[i * BHER_BINS + k];
      norm += steps[n].bher_dist_kin[i * BHER_BINS + k];
    }
    eta_avg /= norm;
    eta_avg *= dc;
    steps[n].bh_eta_kin_avg[i] = log10(eta_avg);
    
  }
}


// Calculate the fraction of kinetically powerful SMBHs in massive halos at low-z.
double recent_kinetic_frac_in_massive_halos(void) 
{
  int64_t n, j, k;
  double nd_tot = 0;
  double frac_tot = 0;
  double norm = BHER_BPDEX;

  // Find the smallest scale factor that enters the window
  for (n=0; n<num_outputs; n++) if (steps[n].scale > KIN_POWER_A_CONSTRAINT_LOW) break;
  for (; n<num_outputs; n++) 
  {
    // If the scale factor is bigger than the upper limit, stop.
    if (steps[n].scale > KIN_POWER_A_CONSTRAINT_HIGH) break;


    // Only count the kinetic powers in massive halos.
    for (j=(KIN_POWER_M_CONSTRAINT-M_MIN)*BPDEX; j<M_BINS; j++) 
    {
      double frac = 0;
      //if (steps[n].t[j]<REAL_ND_CUTOFF) continue;
      if (steps[n].t[j] <= 0) continue;

      // Need to account for the fact that only part of the SMBHs are active, which is parametrized by the duty cycle.
      double dc = steps[n].bh_duty[j];
      if (dc < 1e-4) dc = 1e-4;

      // The critical Eddington ratio is calculated from the kinetic power limit,
      // and further scaled relative to the typical Eddington ratio.
      double eta_crit = KIN_POWER_CONSTRAINT_LOW - (38.1 + steps[n].log_bh_mass[j]) - steps[n].bh_eta[j];
      double f_eta = (eta_crit - BHER_MIN) * BHER_BPDEX;
      int64_t b_eta = f_eta;
      if (f_eta < 0) 
      {
        frac = 1;
      }
      else if (f_eta >= BHER_BINS)
      {
        frac = 0;
      }
      else
      {
        f_eta -= b_eta;
        frac += (1 - f_eta) * steps[n].bher_dist_kin[j * BHER_BINS + b_eta];
        for (k = 0; k < BHER_BINS; k++) 
          {
            if (k < b_eta + 1) continue;
            frac += steps[n].bher_dist_kin[j * BHER_BINS + k];
          }
        frac /= norm;
      }

      frac_tot += frac * dc * steps[n].t[j];
      nd_tot += steps[n].t[j];
    }
  }
  frac_tot /= nd_tot;
  return frac_tot;
}

// Calculate the recent AGN radiative power in massive halos.
double recent_radiative_power_in_massive_halos(void) 
{
  int64_t n, j;
  double L_rad = 0;
  double nd_tot = 0;

  // Find the smallest scale factor that enters the window
  for (n=0; n<num_outputs; n++) if (steps[n].scale > RAD_POWER_A_CONSTRAINT_LOW) break;
  for (; n<num_outputs; n++) 
  {
    // If the scale factor is bigger than the upper limit, stop.
    if (steps[n].scale > RAD_POWER_A_CONSTRAINT_HIGH) break;


    // Only count the kinetic powers in massive halos.
    for (j=(RAD_POWER_M_CONSTRAINT-M_MIN)*BPDEX; j<M_BINS; j++) 
    {
      double L_tmp = 0;
      //if (steps[n].t[j]<REAL_ND_CUTOFF) continue;
      if (steps[n].t[j] <= 0) continue;

      L_tmp = 5.66e46 * steps[n].smhm.bh_efficiency_rad * steps[n].bh_acc_rate[j] * steps[n].bh_f_occ[j];
      L_rad += L_tmp * steps[n].t[j];
      nd_tot += steps[n].t[j];
    
    }

  }

  L_rad /= nd_tot;
  return log10(L_rad);
}

// same as recent_sfh_in_massive_halos, but without the correction
// on the SFR.
double recent_sfh_in_massive_halos_nocorr(void) 
{
  int64_t n, j, count=0;
  double sfr = 0;
  for (n=0; n<num_outputs; n++) if (steps[n].scale > SFR_A_CONSTRAINT) break;
  for (; n<num_outputs; n++) 
  {
    for (j=(SFR_M_CONSTRAINT-M_MIN)*BPDEX; j<M_BINS; j++) 
    {
      if (steps[n].t[j] <= 0) continue;
      sfr += steps[n].sfr[j];
      count++;
    }
  }
  sfr /= count;
  return sfr;
}

double recent_Micl_Mstar_ratio_in_massive_halos(void) 
{
  int64_t n, j;
  double Mstar_tot = 0;
  double Micl_tot = 0;
  double ND_tot = 0;
  for (n=0; n<num_outputs; n++) if (steps[n].scale > ICL_RATIO_A_CONSTRAINT) break;
  for (; n<num_outputs; n++) 
  {
    for (j=(ICL_RATIO_M_CONSTRAINT-M_MIN)*BPDEX; j<M_BINS; j++) 
    {
      Mstar_tot += steps[n].sm_avg[j] * steps[n].t[j];
      Micl_tot += steps[n].sm_icl[j] * steps[n].t[j];
      ND_tot += steps[n].t[j];
    }
  }
  Mstar_tot /= ND_tot;
  Micl_tot /= ND_tot;
  return log10(Micl_tot / Mstar_tot);
}

// Calculate the newly formed stellar masses, SFRs, new BH masses,
// and BH accretion and merger rates.
void calc_new_sm_and_sfr(int n, int i, struct smf_fit *fit) 
{
  double dt;
  int64_t j;
  char buffer[1024];
  dt = steps[n].dt;

  // Ignore the halo mass bins that are too rare.
  if (!steps[n].t[i]) 
  {
    return;
  }

  // New stellar mass is calculated by multiplying SFR by the time interval,
  // with the mass loss due to stellar evolution. Aside from star formation,
  // we also add the merger contribution to it in the lines below.
  steps[n].new_sm[i] = steps[n].sfr[i] * dt * steps[n].smloss[n];

  // // double sm_from_icl = steps[n].smhm.icl_frac * steps[n].sm_icl[i];
  // Instead of having a certain fraction of all the ICL merge into the central galaxy,
  // (the commented line above),
  // we only merge a fraction of ICL that comes from ***the stellar mass*** of incoming
  // satellite galaxies, i.e., steps[n].sm_inc[i]. Here we also fix the small bug that
  // the merged stellar mass does not shrink with time.
  double sm_from_icl = steps[n].smhm.icl_frac * steps[n].sm_inc[i] * steps[n].smloss[n];
  steps[n].mr[i] = sm_from_icl / steps[n].smloss[n] / dt;

  // Add the merged stellar mass.
  steps[n].new_sm[i] += sm_from_icl;

  // The current average stellar mass is simply the sum of old and new stellar masses.
  steps[n].sm_avg[i] = steps[n].old_sm[i] + steps[n].new_sm[i];

  // Assuming a log-normal scatter around the stellar mass--halo mass relation,
  // the median stellar mass is divided by a scatter correction.
  steps[n].sm[i] = steps[n].sm_avg[i] / steps[n].smhm.scatter_corr;
  steps[n].log_sm[i] = (steps[n].sm[i] > 1) ? log10(steps[n].sm[i]) : 0;

  // Calculate the bulge mass with the bulge mass--total stellar mass relation (BMSM).
  // Note that the stellar mass here should be corrected for the systematic offset,
  // since the BMSM is based on ***observed*** data.
  // steps[n].log_bm[i] = bulge_mass(steps[n].log_sm[i]+steps[n].smhm.mu, steps[n].scale);
  steps[n].log_bm[i] = bulge_mass(steps[n].log_sm[i], steps[n].scale);
  steps[n].log_sm_obs[i] = steps[n].log_sm[i]+steps[n].smhm.mu;


  // We are calculating the average and median BH mass at fixed ***halo mass*** bins.
  // So we also need to calculate the scatter in BH mass at fixed halo mass.
  // This scatter is the quadratic sum of the scatter around the black hole mass--
  // bulge mass (BHBM) relation, and that around the stellar mass--halo mass relation,
  // scaled by the slope of the BHBM relation.
  double bh_sm_scatter = steps[n].smhm.scatter*steps[n].smhm.bh_gamma;

  // bh_vdisp_scatter is deprecated for now, since we do not find a good
  // way to use the black hole mass--velocity dispersion (M-sigma) relation.
  double bh_vdisp_scatter = 0.3; //for the time being, use 0.3 dex as the vdisp scatter
  
  // Actual quadratic sum.
  double combined_scatter = (vel_dispersion) ? sqrt(bh_vdisp_scatter  * bh_vdisp_scatter + steps[n].smhm.bh_scatter*steps[n].smhm.bh_scatter) : sqrt(bh_sm_scatter  * bh_sm_scatter + steps[n].smhm.bh_scatter*steps[n].smhm.bh_scatter);
  // The correction between the average and median BH masses assuming a log-normal scatter.
  double bh_scatter_corr = exp(pow(combined_scatter*log(10), 2)/2.0);




  // Calculate median BH mass based on the black hole mass--bulge mass relation.
  // the M-sigma option is deprecated for now.
  steps[n].log_bh_mass[i] = (vel_dispersion) ? 
      calc_bh_at_vdisp(steps[n].vdisp[i], steps[n].smhm)
    : calc_bh_at_bm(steps[n].log_bm[i], steps[n].smhm);

  // the average BH mass is offset from the median by a factor of bh_scatter_corr.
  steps[n].bh_mass_avg[i] = doexp10(steps[n].log_bh_mass[i])*bh_scatter_corr;
  //if (n >= 111 && n <= 152 && i == 9)
  //fprintf(stderr, "n=%d, i=%d, log_sm=%f, log_bm=%f, log_bh_mass=%f, bh_scatter_corr=%e, bh_mass_avg=%e\n", n, i, steps[n].log_sm[i], steps[n].log_bm[i], steps[n].log_bh_mass[i], bh_scatter_corr, steps[n].bh_mass_avg[i]);

  double f_mass = exp((steps[n].med_hm_at_a[i]
                      - steps[n].smhm.dc_mbh) / steps[n].smhm.dc_mbh_w);
  f_mass = f_mass / (1.0 + f_mass);
  steps[n].bh_f_occ[i] = steps[n].smhm.f_occ_min + (1.0 - steps[n].smhm.f_occ_min) * f_mass;
  double dc_pl = doexp10(steps[n].smhm.bh_duty_alpha * (steps[n].med_hm_at_a[i] - steps[n].smhm.bh_duty_m));
  if (dc_pl > 1.0) dc_pl = 1.0;
  //steps[n].bh_duty[i] = steps[n].bh_f_occ[i] * dc_pl;


  steps[n].bh_f_occ_max[i] = 1.0;
  double f_occ_max = 0;
  if (steps[n].scale > 0.142857   && steps[n].med_hm_at_a[i] >= 9)
  {
    double tp = steps[n - 1].t[i];
    double f_occ_p = steps[n - 1].bh_f_occ[i];
    double r = steps[n].merged[i*M_BINS+i];
    f_occ_max = 1.0 / steps[n].t[i] * 
    (steps[n].mmp[i*M_BINS+i] * ((tp - r) * f_occ_p + r * f_occ_p * (2 - f_occ_p)) / (tp - r)
     + steps[n].mmp[(i-1)*M_BINS+i] * steps[n - 1].bh_f_occ[i - 1]);
    steps[n].bh_f_occ[i] = f_occ_max;
    // if (steps[n].bh_f_occ[i] > f_occ_max * (1 + 1e-10))
    // {
    //   fprintf(stderr, 
    //     "Too much increase in SMBH occupation fraction! a: %f, m: %f, tp: %e, r: %e, t: %e, mmp[i,i]: %e, mmp[i-1,i]: %e, f_occ_p: %e, f_occ_i-1: %e, max f_occ: %e, new f_occ: %e\n", 
    //     steps[n].scale, steps[n].med_hm_at_a[i], tp, r, steps[n].t[i], steps[n].mmp[i*M_BINS+i], steps[n].mmp[(i-1)*M_BINS+i], f_occ_p, steps[n - 1].bh_f_occ[i - 1], f_occ_max, steps[n].bh_f_occ[i]);
    //   INVALIDATE(fit, buffer); 
    // }
  }
  steps[n].bh_duty[i] = steps[n].bh_f_occ[i] * dc_pl;

  // steps[n].bh_f_occ[i] = 1.0;
  steps[n].bh_mass_avg[i] /= steps[n].bh_f_occ[i];
  steps[n].log_bh_mass[i] -= log10(steps[n].bh_f_occ[i]);
  // We also have to correct old and wandering BH mass to account for the occupation fraction!
  steps[n].old_bh_mass[i] /= steps[n].bh_f_occ[i];
  steps[n].bh_unmerged[i] /= steps[n].bh_f_occ[i];

  // The new BH mass is the difference between the current and old BH masses.
  float new_bh_mass = steps[n].bh_mass_avg[i] - steps[n].old_bh_mass[i];

  // mfrac is the fraction of merger contribution to the total BH mass growth,
  // which is proportional to the fractional merger contribution to ***galaxy*** growth.
  // The proportionality factor is parametrized as steps[n].smhm.f_merge_bh.
  // See smhm_at_z() for the redshift scaling of steps[n].smhm.f_merge_bh.
  float mfrac = n ? steps[n].smhm.f_merge_bh * sm_from_icl / steps[n].new_sm[i] : 0; 

  // Calculate the accretion and merger rates based on the total BH mass growth
  // and merger contributions. Note that the average BHAR and BHMR are now for
  // SMBH host halos, not for all halos.
  steps[n].new_bh_mass[i] = new_bh_mass;
  steps[n].bh_merge_rate[i] = mfrac * new_bh_mass / dt;
  steps[n].bh_acc_rate[i] = (1 - mfrac) * new_bh_mass /dt;

  if (steps[n].scale > 0.08 && ((steps[n].med_hm_at_a[i] > 11 && steps[n].log_bh_mass[i] >= steps[n].log_sm[i] - 0.5)))
  {
   //      fprintf(stderr, "Too big Mbh compared to Mstar. a: %f, m: %f, mstar: %f, mbh: %f, mbh_i-1: %f, n=%d, i=%d\n", steps[n].scale, steps[n].med_hm_at_a[i], steps[n].log_sm[i], steps[n].log_bh_mass[i], steps[n].log_bh_mass[i-1], n, i);
         INVALIDATE(fit, buffer);
  }


  // Since our parametrization cannot guarantee that the new BH mass is positive,
  // ill-behaved parameter sets are discarded if they produce negative BH mass
  // growth (or accretion rates), and the problematic BHs exist at z < 11.5, with
  // masses above BH_MASS_TO_REQUIRE_GROWTH.
  if ((!(new_bh_mass>=0) || !(steps[n].bh_acc_rate[i]>0)) && (steps[n].scale > 0.08 && i>BPDEX && (steps[n].log_bh_mass[i] > BH_MASS_TO_REQUIRE_GROWTH || steps[n].med_hm_at_a[i] > 11)))
  {
   //fprintf(stderr, "Negative BH accretion rate. a: %f, m: %f, new_bh_mass: %f, old_bh_mass: %f, bh_f_occ: %e, dc_mbh: %f, dc_mbh_w:%f, focc_min: %e, focc_max: %e, mbh_mmp: %e, focc_mmp: %e, nd_mmp: %e, nd_last: %e\n", steps[n].scale, steps[n].med_hm_at_a[i], steps[n].bh_mass_avg[i], steps[n].old_bh_mass[i], steps[n].bh_f_occ[i], steps[n].smhm.dc_mbh, steps[n].smhm.dc_mbh_w, steps[n].smhm.f_occ_min, f_occ_max, steps[n-1].bh_mass_avg[i-1], steps[n-1].bh_f_occ[i-1], steps[n].mmp[(i-1)*M_BINS+i], steps[n-1].t[i-1]);
    INVALIDATE(fit, buffer); 
  }

  // There's a finite amount of unmerged BH masses for central black holes to sonsume
  // during mergers. If the merger fraction is so big that these unmerged BHs cannot
  // cover the merger rate, then we need to invalidate the model.
  if (steps[n].bh_unmerged[i] < mfrac*new_bh_mass) 
  {
    if (steps[n].scale > 0.08 && i>BPDEX && (steps[n].log_bh_mass[i] > BH_MASS_TO_REQUIRE_GROWTH || steps[n].med_hm_at_a[i] > 11))
    {
   //   fprintf(stderr, "Merging rate exceeds available mergers! (sm: %e, m: %e, scale: %f; old_bh: %e; new_bh: %e! DBH: %e; Merge: %e; M_avail: %e \n",
     // steps[n].sm_avg[i], exp10(M_MIN+(i+0.5)/BPDEX), steps[n].scale, steps[n].old_bh_mass[i], steps[n].bh_mass_avg[i], new_bh_mass, new_bh_mass*mfrac, steps[n].bh_unmerged[i]);
       INVALIDATE(fit, buffer); 
    } 
    else 
    {
      steps[n].bh_merge_rate[i] = steps[n].bh_unmerged[i] / dt;
      steps[n].bh_acc_rate[i] = (new_bh_mass - steps[n].bh_unmerged[i]) / dt;
      steps[n].bh_unmerged[i] = 0;
      //for (j=0; j<MBH_BINS; j++) steps[n].bh_unmerged_dist[i*MBH_BINS+j] = 0;
    } 
  }

  // Remove the merged BH masses from the unmerged (or wandering) reservoir of BHs.
  else 
  {
    steps[n].bh_unmerged[i] -= mfrac*new_bh_mass;
    // for (j=0; j<MBH_BINS; j++) steps[n].bh_unmerged_dist[i*MBH_BINS+j] *= steps[n].bh_unmerged[i] / (steps[n].bh_unmerged[i] + mfrac*new_bh_mass);
  }

  // Calculate the ***average*** Eddington ratio for this halo mass bin.
  double bhar_tmp = steps[n].bh_acc_rate[i] > 0 ? steps[n].bh_acc_rate[i] : 1e-8;
  steps[n].bh_eta[i] = log10(bhar_tmp/steps[n].bh_mass_avg[i]*4.5e8*(steps[n].smhm.bh_efficiency_rad));



  steps[n].merged_frac[i] = 0;
  if (steps[n].smhm.icl_frac && steps[n].sm_icl[i]) 
  {
    steps[n].merged_frac[i] = steps[n].smhm.icl_frac; //merged_frac is the fraction of the incoming satellite stellar mass
                                                      //that merge into the central galaxies, which under this parametrization
                                                      //is simply icl_frac.
    // Since galaxy mergers involve the transfer from ICL to central stellar mass,
    // we need to transfer the mass histories as well.
    for (j=0; j<n; j++) 
    {
      steps[n].sm_hist[i*num_outputs + j] += steps[n].smhm.icl_frac*steps[n].sm_inc_hist[i*num_outputs + j];
      steps[n].icl_stars[i*num_outputs + j] -= steps[n].smhm.icl_frac*steps[n].sm_inc_hist[i*num_outputs + j];
    }
    // Give back the merger contributions in the new_sm, so that we can use it
    // to calculate the latest stellar mass history (see below)
    steps[n].new_sm[i] -= (steps[n].smhm.icl_frac * steps[n].sm_inc[i] * steps[n].smloss[n]);
  }

  // Calculate the latest stellar mass history.
  steps[n].new_sm[i] /= steps[n].smloss[n];
  steps[n].sm_hist[i*num_outputs + n] = steps[n].new_sm[i];

  // If the average stellar mass is higher than the overall baryonic fraction in the halo,
  // we should invalidate the model.
  if ((steps[n].t[i] > REAL_ND_CUTOFF) && (steps[n].sm_avg[i] > 0.17*exp10(steps[n].med_hm_at_a[i])) && (!no_z_scaling) && n>2) 
  {
    //fprintf(stderr, "SM exceeds baryon fraction (sm: %e, m: %e, scale: %f!\n",
    //steps[n].sm_avg[i], exp10(M_MIN+(i+0.5)*INV_BPDEX), steps[n].scale);
    INVALIDATE(fit, buffer);
  }

  // If we have a negative SFR, also invalidate the model.
  // However, this is not gonna happen in the current parametrization,
  // since SFRs are calculated from the SFR--Vmax relation and won't be
  // negative. In previous parametrizations (see Behroozi et al. 2013),
  // stellar masses are determined and the SFRs are calculated by differentiation.
  // Only in that case this unphysical situation could occur.
  if ((!no_z_scaling) && (steps[n].sfr[i] < 0 || !isfinite(steps[n].sfr[i]))) 
  {
    char buffer[1024];
    //fprintf(stderr, "Negative SFR at a = %f and m = %f (ND=%g)! (ICL m=%f; ICL frac=%f)\n", steps[n].scale, 
    //i*INV_BPDEX + M_MIN, steps[n].t[i], steps[n].smhm.icl_frac, steps[n].smhm.icl_frac);
    INVALIDATE(fit, buffer);
  }
}

