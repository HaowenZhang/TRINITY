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

static inline double doexp10(double x) {
  double a = exp(M_LN10*x);
  return a;
}
static inline double dolog10(double x) {
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

void calc_sfh(struct smf_fit *f) {
  //fprintf(stderr, "z_start: %f", 1 / steps[0].scale - 1);
   int64_t i,j;
  for (i=0; i<num_outputs; i++) {
    if (!no_z_scaling || ((steps[i].scale < 1.0/(z_min+1.0)) &&
			  (steps[i].scale > 1.0/(z_max+1.0))))
    calc_sm_hist(i, f);
    //for (j=0; j<M_BINS; j++)
    //{
      //fprintf(stdout, "%d %f %f\n", i, steps[i].med_hm_at_a[j], steps[i].log_sm[j]);
    //}
    
  }

#pragma omp for schedule(guided,5)
  for (j=1; j<num_outputs; j++) {
    if (!no_z_scaling || ((steps[j].scale < 1.0/(z_min+1.0)) &&
			  (steps[j].scale > 1.0/(z_max+1.0)))) {
      calc_smf_and_ssfr(j, f);
      calc_uvlf(j, f);
      calc_total_sfr(j);

      calc_bh_acc_rate_distribution(j, f);
      calc_bh_acc_rate_distribution_full(j, f);
      calc_bh_lum_distribution_full(j, f);
      calc_active_bh_fraction(j, f);
      calc_avg_eta_rad(j);
      calc_total_bhar(j);
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
  double ledd_max = steps[n].ledd_max[i];
  double bpdex = steps[n].ledd_bpdex[i];
  double inv_bpdex = 1.0/bpdex;

  
  double scatter_tot = sqrt(steps[n].smhm.bh_scatter * steps[n].smhm.bh_scatter +
                            steps[n].smhm.scatter * steps[n].smhm.bh_gamma * 
                            steps[n].smhm.scatter * steps[n].smhm.bh_gamma);
  double mbh_min = steps[n].log_bh_mass[i] - 8 * scatter_tot;
  double mbh_max = steps[n].log_bh_mass[i] + 8 * scatter_tot;
  double dmbh = 0.05;
  double norm_gauss = 1 / sqrt(2 * M_PI * scatter_tot * scatter_tot);
  double bher_norm = 0;

  for (int j=0; j<BHER_BINS; j++) bher_norm += steps[n].bher_dist[i * BHER_BINS + j];

  double prob_l = 0;
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
  double z = 1 / steps[n].scale - 1;
  double mass_real = 13.5351-0.23712*z+2.0187*exp(-z/4.48394);
}


void remove_bh_merger()
{
  for (int i=1; i<num_outputs; i++)
  {
    for (int j=0; j<M_BINS; j++)
    {
      if (!steps[i].t[j]) continue;
      
      if (isfinite(steps[i].new_bh_mass[j] * steps[i].bh_merge_rate[j] / \
                                (steps[i].bh_acc_rate[j] + steps[i].bh_merge_rate[j])))
      steps[i].bh_unmerged[j] += steps[i].new_bh_mass[j] * steps[i].bh_merge_rate[j] / \
                                (steps[i].bh_acc_rate[j] + steps[i].bh_merge_rate[j]);
     
      if (isfinite(steps[i].bh_acc_rate[j] / \
                                (steps[i].bh_acc_rate[j] + steps[i].bh_merge_rate[j])))
      steps[i].new_bh_mass[j] *= steps[i].bh_acc_rate[j] / \
                                (steps[i].bh_acc_rate[j] + steps[i].bh_merge_rate[j]);

     
      double new_bh_mass_avg = steps[i].old_bh_mass[j] + steps[i].new_bh_mass[j];
      double offset = log10(steps[i].bh_mass_avg[j]) - steps[i].log_bh_mass[j];

      steps[i].log_bh_mass[j] += log10(new_bh_mass_avg / steps[i].bh_mass_avg[j]);

      if (!isfinite(steps[i].log_bh_mass[j])) steps[i].log_bh_mass[j] = -5;
      steps[i].bh_mass_avg[j] = new_bh_mass_avg;
      if (!finite(steps[i].bh_mass_avg[j])) steps[i].bh_mass_avg[j] = exp10(-5 + offset);
      if (!finite(steps[i].bh_unmerged[j])) steps[i].bh_unmerged[j] = 0;

      steps[i].bh_merge_rate[j] = 0;

    }


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

void calc_total_sfr(int n) {
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
 double bhar_minimum = -1;
 double l_min = 41.5 - 2.5 * log10(steps[n].scale); // redshift-dependent lower limit of 2-10KeV QLF
 // double l_min = 42.5;
 double BC = 10.83 * pow(pow(10, l_min) / 3.86e43, 0.28) + 6.08 * pow(pow(10, l_min) / 3.86e43, -0.02);
 l_min += log10(BC); //Convert to bolometric
l_min = 90 - 2.5 * l_min; //Convert it to magnitude
 for (i = 0; i < M_BINS; i++)
 {
   if (steps[n].t[i] == 0) continue;
   //if (steps[n].bh_mass_avg[i] < 1e5) continue;
   // if (steps[n].bh_acc_rate[i] < bhar_minimum) continue;
   // double eta = -(l_min + 5.26) / 2.5 - steps[n].log_bh_mass[i];
  double eta = -(l_min + 5.26) / 2.5 - log10(steps[n].bh_mass_avg[i]);
   double eta_frac = eta - steps[n].bh_eta[i];
   bhar += steps[n].bh_acc_rate[i] * steps[n].t[i];
   if (eta_frac > steps[n].ledd_max[i]) continue;
   if (eta_frac < steps[n].ledd_min[i]) 
    {
      bhar_obs += steps[n].bh_acc_rate[i] * steps[n].t[i];
    }
   else
   {
    double bher_f = (eta_frac - steps[n].ledd_min[i]) * steps[n].ledd_bpdex[i];
    int64_t bher_b = bher_f;
    bher_f -= bher_b;
    double p1 = 0, p2 = 0;
    if (i*BHER_BINS+bher_b+1 < M_BINS * BHER_BINS)
    {
      p1 = steps[n].bher_dist_full[i*BHER_BINS+bher_b];
      p2 = steps[n].bher_dist_full[i*BHER_BINS+bher_b+1];
    }
    else if (i*BHER_BINS+bher_b < M_BINS * BHER_BINS)
    {
      p1 = steps[n].bher_dist_full[i*BHER_BINS+bher_b];
      p2 = p1;
    }
    
    // if (bher_b >= BHER_BINS-1) p2 = p1;
    double prob = p1 + bher_f*(p2-p1);
    for (;bher_b + 1< BHER_BINS; bher_b++) prob += steps[n].bher_dist_full[i*BHER_BINS+bher_b+1];
    double total = 0;
    for (bher_b = 0; bher_b < BHER_BINS; bher_b++) total += steps[n].bher_dist_full[i*BHER_BINS+bher_b];
      //fprintf(stderr, "prob=%f, total=%f\n", prob, total);
    if (prob > 0) prob /= total;
    
    // double f_mass = exp((log10(steps[n].bh_mass_avg[i]) - steps[n].smhm.dc_mbh) / steps[n].smhm.dc_mbh_w);
    // f_mass = f_mass / (1 + f_mass);
    // double dc = steps[n].smhm.bh_duty * f_mass;
    // if (dc < 1e-4) dc = 1e-4;
    // prob *= dc;
    bhar_obs += steps[n].bh_acc_rate[i] * steps[n].t[i] * prob;
   }
 }
 steps[n].cosmic_bhar = bhar;
 steps[n].observed_cosmic_bhar = bhar_obs;
}

// Also calculates observable duty cycle.
void calc_observed_bhar(int n)
{
 int64_t i;
 double bhar = 0;
 double bhar_minimum = -1;
 double l_min = 41.5 - 2.5 * log10(steps[n].scale); // redshift-dependent lower limit of 2-10KeV QLF
 // double l_min = 42.5;
 double BC = 10.83 * pow(pow(10, l_min) / 3.86e43, 0.28) + 6.08 * pow(pow(10, l_min) / 3.86e43, -0.02);
 l_min += log10(BC); //Convert to bolometric
l_min = 90 - 2.5 * l_min; //Convert it to magnitude
 for (i = 0; i < M_BINS; i++)
 {
   if (steps[n].t[i] == 0) continue;
   // if (steps[n].bh_acc_rate[i] < bhar_minimum) continue;
   // double eta = -(l_min + 5.26) / 2.5 - steps[n].log_bh_mass[i];
   double dc = steps[n].smhm.bh_duty;
    double f_mass = exp((log10(steps[n].bh_mass_avg[i]) - steps[n].smhm.dc_mbh) / steps[n].smhm.dc_mbh_w);
    f_mass = f_mass / (1 + f_mass);
    // f_mass = f_mass < 1? f_mass : 1;
    dc *= f_mass;
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

struct smf smhm_at_z(double z, struct smf_fit f) {
  struct smf c;
  double a = 1.0/(1.0+z);
  double a1 = -z/(1.0+z);
  double incompleteness;
  double a2 = a1*a1;

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
  c.delta = DELTA(f); // + (a1*DELTA_A(f) + z8*DELTA_A2(f))*expscale;
  c.beta = BETA(f) + a1*BETA_A(f) + z8*BETA_A2(f);
  c.gamma = 0;
  c.lambda = doexp10(LAMBDA(f) + (a1*LAMBDA_A(f) + z8*LAMBDA_A2(f)));
  c.mu = MU(f) + a1*MU_A(f);
  c.kappa = KAPPA(f) + a1*KAPPA_A(f);
  c.passive_mass = 10.3 + z*0.5 - c.mu;
  c.scatter = SCATTER(f) + z*SCATTER_A(f);
  if (c.scatter > 0.4) c.scatter = 0.4;
  c.scatter_corr = exp(pow(c.scatter*log(10), 2)/2.0);
  c.obs_scatter = (no_obs_scatter) ? 0 : SIGMA_CENTER + SIGMA_Z(f)*(z);
  if (c.obs_scatter > 0.3) c.obs_scatter = 0.3;

  c.icl_frac = exp10(ICL_FRAC(f) + a1 * ICL_FRAC_E(f) + z8 * ICL_FRAC_Z(f));

  incompleteness = BURST_DUST_AMP(f)/(1+exp(BURST_DUST_Z(f)-z));
  if (z < 1.0) incompleteness = 0;
  if (z > 1.0) incompleteness -= BURST_DUST_AMP(f)/(1+exp(BURST_DUST_Z(f)-1.0));
  if (incompleteness < 0) incompleteness = 0;
  c.sm_completeness = 1.0 - incompleteness;
  c.csfr_completeness = 1.0 - (1.0-BURST_DUST(f))*incompleteness;
  c.sfr_sm_corr = 1.0 + (4.0*RHO_05(f)-3.23)*a + (2.46-4.0*RHO_05(f))*a*a;
  if (c.csfr_completeness < 0.01) c.csfr_completeness = 0.01;
  if (c.sm_completeness < 0.01) c.sm_completeness = 0.01;
  c.ssfr_corr = 1.0/(1.0-BURST_DUST(f)*incompleteness);
  c.sm_max = c.sm_min = 0;
  c.f_1 = dolog10(2)-c.delta*pow(dolog10(2), c.gamma)/(1.0+exp(1));
  c.lm_slope = 1.0/c.alpha - 1.0;
  c.mpk = c.mu+c.kappa;
  c.combined_scatter = 0;
  c.valid = 1;

  c.qm = QM_0(f) + QM_1(f) *a1 + QM_2(f) *z8;
  c.qwidth = QWIDTH_0(f) + QWIDTH_1(f) *a1 + QWIDTH_2(f) * z8;

  // // BH part
  c.bh_beta = BH_BETA_0(f) +BH_BETA_1(f)*a1 + BH_BETA_2(f)*z8;
  c.bh_gamma = BH_GAMMA_0(f) + BH_GAMMA_1(f)*a1 + BH_GAMMA_2(f)*z8;
  c.f_merge_bh = exp10(BH_MERGE_F_0(f) + BH_MERGE_F_1(f)*a1);
  c.bh_merge_width = BH_MERGE_W_0(f) + BH_MERGE_W_1(f)*a1;
  c.bh_alpha = BH_ALPHA_0(f) + BH_ALPHA_1(f)*a1;
  c.bh_delta = BH_DELTA_0(f) + BH_DELTA_1(f)*a1;


  c.abhmf_shift = ABHMF_SHIFT(f);
  c.bh_efficiency_rad = exp10(BH_EFFICIENCY_0(f) + a1 * BH_ETA_CRIT_1(f));
  c.bh_eta_crit = BH_ETA_CRIT_0(f)+BH_ETA_CRIT_1(f)*a1;
  c.bh_duty = BH_DUTY_0(f)+BH_DUTY_1(f)*a1;
  if (c.bh_duty < 1e-4) c.bh_duty = 1e-4;
  c.bh_scatter = BH_SCATTER_0(f) + BH_SCATTER_1(f)*a1;

  c.dc_mbh = DC_MBH_0(f);
  c.dc_mbh_w = DC_MBH_W_0(f);
  c.eta_mu = ETA_MU_0(f);

  return c;
}

double calc_bh_at_bm(double bm, struct smf c) {
  return c.bh_beta + c.bh_gamma * (bm - 11.0);
}

double calc_bh_at_vdisp(double vd, struct smf c) {
  return c.bh_beta + c.bh_gamma * (vd - 2.30103);
}

double calc_sm_at_m(double m, struct timestep steps) {
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

double calc_sfr_at_lv(double m, double lv, struct smf c) {
  double vd = lv - c.v_1;
  double vd2 = vd/c.delta;
  
  double sfrac = 1 / (1 + exp((lv - c.qm) / c.qwidth));
  return (sfrac*c.epsilon * (1.0/(doexp10(c.alpha*vd) + doexp10(c.beta*vd)) + c.gamma*exp(-0.5*vd2*vd2)));
}

double calc_sfrac_at_m(double m, struct smf c) {
  double sfrac = 1 / (1 + exp((m - c.qm) / c.qwidth));
  return (sfrac);
}

// double calc_sfrac_at_lv(double lv, struct smf c) {
//   double vd3 = (lv-c.qv)/c.qsig;
//   double qfrac = c.fqmin + (1.0 - c.fqmin) * (0.5 + 0.5 * erf(vd3*M_SQRT1_2));
//   double sfrac = 1.0 - qfrac;
//   return (sfrac);
// }


double calc_smf_at_m(double m, double sm, int64_t n, int64_t n_spline, struct smf_fit *fit, gsl_interp_accel *ga) {
  if (sm < steps[n].smhm.sm_min) 
  {
    return 1e-17;
  }
  if (sm > steps[n].smhm.sm_max) 
  {
    return 1e-17;
  }
  if (!steps[n].smhm.valid) 
  {
    return 1e-17;
  }

  if (m < M_MIN || m > M_MAX) return 1e-17;

  double dlm_dlsm;

  int err; 
  if (n_spline == 1) err = gsl_spline_eval_deriv_e(steps[n].spline, sm, ga, &dlm_dlsm);
  else if (n_spline == 2) err = gsl_spline_eval_deriv_e(steps[n].spline2, sm, ga, &dlm_dlsm);

  if (err || !isfinite(dlm_dlsm))
  {
    return 0;
  }
  dlm_dlsm = fabs(dlm_dlsm);
  double dndlogm = mf_cache(steps[n].scale, m);
  double phi = doexp10(dndlogm)*dlm_dlsm;
  if (!isfinite(phi)) phi = 0;
  return phi;
}

double calc_smf_at_sm(int64_t n, double sm) {
  if (sm < steps[n].smhm.sm_min) return 1e-17;
  if (sm > steps[n].smhm.sm_max) return 1e-17;
  if (!steps[n].smhm.valid) return 1e-17;
  gsl_interp_accel ga = {0};
  double m;
  int err = gsl_spline_eval_e(steps[n].spline, sm, &ga, &m);

  if (err)
    {
      return 0;
    }
  double result = calc_smf_at_m(m, sm, n, 1, NULL, &ga);
  gsl_interp_accel_reset(&ga);
  if (steps[n].alloc2smf)
  {
    err = gsl_spline_eval_e(steps[n].spline2, sm, &ga, &m); //Second segment.
    if (!err) result += calc_smf_at_m(m, sm, n, 2, NULL, &ga);
  }
  
  return result;
}




double calc_uvlf_at_m(double m, double uv, int64_t n, int64_t n_spline, struct smf_fit *fit, gsl_interp_accel *ga) {
  if (uv < steps[n].smhm.uv_min) 
  {
    return 1e-17;
  }
  if (uv > steps[n].smhm.uv_max) 
  {
    return 1e-17;
  }
  if (!steps[n].smhm.valid) 
  {
    return 1e-17;
  }
  if (m < M_MIN || m > M_MAX) return 1e-17;

  double dlm_dMuv;
  int err;
  if (n_spline == 1) err = gsl_spline_eval_deriv_e(steps[n].spline_uv, uv, ga, &dlm_dMuv);
  else if (n_spline == 2) err = gsl_spline_eval_deriv_e(steps[n].spline_uv2, uv, ga, &dlm_dMuv);

  dlm_dMuv = fabs(dlm_dMuv);

  if (err || !isfinite(dlm_dMuv))
  {
    return 0;
  }
  double dndlogm = mf_cache(steps[n].scale, m);
  if (!(dlm_dMuv>0)) {
    dlm_dMuv = 1.0;
  }

  double phi = doexp10(dndlogm)*dlm_dMuv;
  if (!isfinite(phi)) phi = 0;
  return phi;
}

double calc_uvlf_at_uv(int64_t n, double uv) {
  if (uv < steps[n].smhm.uv_min) return 1e-17;
  if (uv > steps[n].smhm.uv_max) return 1e-17;
  if (!steps[n].smhm.valid) return 1e-17;
  gsl_interp_accel ga = {0};

  double m;
  int err = gsl_spline_eval_e(steps[n].spline_uv, uv, &ga, &m);
  if (err)
    {
      return 0;
    }
  double result = calc_uvlf_at_m(m, uv, n, 1, NULL, &ga);
  gsl_interp_accel_reset(&ga);
  if (steps[n].alloc2uvlf)
  {
    err = gsl_spline_eval_e(steps[n].spline_uv2, uv, &ga, &m); //Remember to add the second segment.
    if (!err) result += calc_uvlf_at_m(m, uv, n, 2, NULL, &ga);
  }
  return result;
}

void calc_active_bh_fraction(int n, struct smf_fit *fit) {
  int64_t i, j;
  for (i=0; i<M_BINS; i++) {

    double ledd_min = steps[n].ledd_min[i];
    double ledd_max = steps[n].ledd_max[i];

    double bpdex = steps[n].ledd_bpdex[i];
    double inv_bpdex = 1.0/bpdex;

    double eta_frac = -2.0 - steps[n].bh_eta[i];

    if (eta_frac < ledd_min || eta_frac > ledd_max) continue; 
    
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

    double dc = steps[n].smhm.bh_duty;
    double f_mass = exp((log10(steps[n].bh_mass_avg[i]) - steps[n].smhm.dc_mbh) / steps[n].smhm.dc_mbh_w);
    f_mass = f_mass / (1 + f_mass);
    dc *= f_mass;
    if (dc < 1e-4) dc = 1e-4;

    steps[n].f_active[i] = p_good / p_total * dc;
  }
}

void calc_active_bh_fraction_lim(int n, struct smf_fit *fit, int ledd_or_lum, double lim) {
  int64_t i, j;
  for (i=0; i<M_BINS; i++) {
    double ledd_lim = ledd_or_lum ? lim : lim - 38.1 - steps[n].log_bh_mass[i];
    double ledd_min = steps[n].ledd_min[i];
    double ledd_max = steps[n].ledd_max[i];

    double bpdex = steps[n].ledd_bpdex[i];
    double inv_bpdex = 1.0/bpdex;

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

    double dc = steps[n].smhm.bh_duty;
    double f_mass = exp((log10(steps[n].bh_mass_avg[i]) - steps[n].smhm.dc_mbh) / steps[n].smhm.dc_mbh_w);
    f_mass = f_mass / (1 + f_mass);
    // f_mass = f_mass < 1? f_mass : 1;
    dc *= f_mass;
    if (dc < 1e-4) dc = 1e-4;

    steps[n].f_active[i] = p_good / p_total * dc;
  }
}

void calc_bh_acc_rate_distribution(int n, struct smf_fit *fit) {
  double l_min = 36; double l_max = 50; //The lower/upper limits of the luminosities for
                                        //the calculation of luminosity distribution at
                                        //different BH masses
  int64_t i, j;
  memset(steps[n].bher_dist, 0, sizeof(double)*BHER_BINS);
  steps[n].ledd_min_abs = 100; 
  steps[n].ledd_max_abs = -100;
  double sm_scatter = steps[n].smhm.scatter*steps[n].smhm.bh_gamma;
  double scatter = sqrt(sm_scatter*sm_scatter + steps[n].smhm.bh_scatter*steps[n].smhm.bh_scatter);


  double bher_min = BHER_EFF_MIN;
  //if (nonlinear_luminosity) bher_min = (bher_min - 2) / 2.0;  

  for (i=0; i<M_BINS; i++)
  {
    // mass-dependent modulation of duty cycle
    // double bh_eta_corr = log10(schechter_inv_avg(steps[n].smhm.bh_alpha, bher_min, BHER_EFF_MAX, fit)/(steps[n].smhm.bh_duty * pow(log10(steps[n].bh_mass_avg[i]) / log10(steps[n].bh_mass_avg[M_BINS - 1]), steps[n].smhm.dc_beta)));
    // double bh_eta_corr = log10(schechter_inv_avg(steps[n].smhm.bh_alpha, bher_min, BHER_EFF_MAX, fit)/(steps[n].smhm.bh_duty * pow((M_MIN + (i + 0.5) * INV_BPDEX) - (M_MAX - 0.5 * INV_BPDEX), steps[n].smhm.dc_beta)));
    // double f_mass = pow(exp10((M_MIN + (i + 0.5) * INV_BPDEX) - (M_MAX - 0.5 * INV_BPDEX)), steps[n].smhm.dc_beta);
    // f_mass = f_mass < 1? f_mass : 1;
    double f_mass = exp((log10(steps[n].bh_mass_avg[i]) - steps[n].smhm.dc_mbh) / steps[n].smhm.dc_mbh_w);
    f_mass = f_mass / (1 + f_mass);
    double dc = (steps[n].smhm.bh_duty * f_mass);
    if (dc < 1e-4) dc = 1e-4;
    // printf("the %d-th mass bin.\n", i);
    // double bh_eta_corr = log10(schechter_inv_avg(steps[n].smhm.bh_alpha, bher_min, BHER_EFF_MAX, fit)/dc);
    
    // Note that here we are going to calculate integrals of 1 / (x**a + x**b) dx and 1 / (x**(a+1) + x**(b+1)) dx,
    // while the doublePL_norm is implemented to calculate the norm assuming dP/dlogx = 1 / (x**a + x**b), so there is
    // an additional -1 offset in the following power-law indices!!!!!!
    double nom = doublePL_norm(steps[n].smhm.bh_alpha + 1 - 1, steps[n].smhm.bh_delta + 1 - 1, bher_min, BHER_EFF_MAX, fit);
    double dnom = doublePL_norm(steps[n].smhm.bh_alpha - 1, steps[n].smhm.bh_delta - 1, bher_min, BHER_EFF_MAX, fit);
    double bh_eta_corr = log10(nom/dnom/dc);
    // printf("alpha=%e, delta=%e, nominator=%e, denominator=%e, bh_eta_corr=%e for %d-th bin at %d-th snapshot\n", 
    //             steps[n].smhm.bh_alpha, steps[n].smhm.bh_delta, nom, dnom, log10(nom/dnom), i, n);
    //if (n == 63 && i == 19) fprintf(stderr, "avg bh_eta=%f, bh_eta_corr=%f, bh_eta0 = %f\n", steps[n].bh_eta[i], bh_eta_corr, steps[n].bh_eta[i] + bh_eta_corr);
    //if (n == 63 && i == 19) fprintf(stderr, "avg bh_eta=%f, alpha=%f, delta=%f, nom=%e, dnom=%e, dc=%e, bh_eta_corr=%f, bh_eta0 = %f\n", steps[n].bh_eta[i], steps[n].smhm.bh_alpha, steps[n].smhm.bh_delta, nom, dnom, dc, bh_eta_corr, steps[n].bh_eta[i] + bh_eta_corr);
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

  }
  steps[n].bh_mass_min = 2;
  steps[n].bh_mass_max = 11;

  // if (scatter <= 0) {
  //   for (i=0; i<BHER_BINS; i++) {
  //     double bher = BHER_MIN+i*BHER_INV_BPDEX;
  //     // steps[n].bher_dist[i] = exp10(bher*steps[n].smhm.bh_alpha)*exp(-exp10(bher));
  //     steps[n].bher_dist[i] = 1 / (exp10(bher*steps[n].smhm.bh_alpha) * exp10(bher*steps[n].smhm.bh_delta));
  //   }
  // } else {
  //   const int64_t cmin = -BHER_BINS-3*BHER_BPDEX;
  //   const int64_t cmax = 3*BHER_BPDEX + BHER_BINS;
  //   double *ecache = check_realloc(NULL, sizeof(double)*(cmax-cmin+1), "Allocating ecache");
  //   double inv_scatter = 1.0/scatter;
  //   for (j=cmin; j<=cmax; j++) {
  //     double db = inv_scatter*j*BHER_INV_BPDEX;
  //     ecache[j-cmin] = exp(-0.5*db*db);
  //   }
  //   for (j=-3*BHER_BPDEX; j<BHER_BINS+3*BHER_BPDEX; j++) {
  //     double bher = BHER_MIN+j*BHER_INV_BPDEX;
  //     double prob = 1 / (exp10(bher*steps[n].smhm.bh_alpha) * exp10(bher*steps[n].smhm.bh_delta));
  //     double *ec = ecache + (-j - cmin);
  //     for (i=0; i<BHER_BINS; i++) {
  // //double bher2 = BHER_MIN+i*BHER_INV_BPDEX;
  // //double db = inv_scatter*(bher - bher2);
  // double weight = ec[i]; //exp(-0.5*db*db);
  // //printf("%e %e %"PRId64" %"PRId64" %f %"PRId64"\n", weight, ec[i], i, j, db, i-j-cmin);  
  // steps[n].bher_dist[i] += weight*prob;
  //     }
  //   }
  //   free(ecache);
  // }

  // //Normalize;
  // double total = 0;
  // for (i=0; i<BHER_BINS; i++) total += steps[n].bher_dist[i];
  // total *= BHER_INV_BPDEX;
  // //assert(total > 0);
  // if (total > 0)
  // {
  // total = 1.0 / total;
  // }
  // for (i=0; i<BHER_BINS; i++) steps[n].bher_dist[i] *= total;

  /*  char buffer[1024];
  sprintf(buffer, "bher_dist/dist_%f.dat\n", steps[n].scale);
  FILE *out = check_fopen(buffer, "w");
  fprintf(out, "#ER Prob.\n");
  for (i=0; i<BHER_BINS; i++) fprintf(out, "%f %e\n", BHER_MIN+i*BHER_INV_BPDEX, steps[n].bher_dist[i]);
  fclose(out);
  */
}


void calc_bh_acc_rate_distribution_full(int n, struct smf_fit *fit) {
  int64_t i, j, k;
  // memset(steps[n].bher_dist_full, 0, sizeof(float)*BHER_BINS*M_BINS);
  memset(steps[n].bher_dist, 0, sizeof(double)*BHER_BINS*M_BINS);
  memset(steps[n].bher_dist_norm, 0, sizeof(double)*M_BINS);

  double bh_eta_crit = steps[n].smhm.bh_eta_crit;
  double sm_scatter = steps[n].smhm.scatter*steps[n].smhm.bh_gamma;
  double scatter = sqrt(sm_scatter*sm_scatter + steps[n].smhm.bh_scatter*steps[n].smhm.bh_scatter);

  //  Already corrected in bh_acc_rate_distribution...
  //  double bh_eta_corr = log10(schechter_inv_avg(steps[n].smhm.bh_alpha, BHER_MIN, BHER_EFF_MAX, fit)/steps[n].smhm.bh_duty);
  //  for (i=0; i<M_BINS; i++) steps[n].bh_eta[i] += bh_eta_corr;

  double (*bher_prob)(double, double, double, double, double, double) = &_prob_of_ledd_linear;
  if (nonlinear_luminosity) bher_prob = &_prob_of_ledd_nonlinear;
  //assert(bher_prob == &_prob_of_ledd_linear);

  if (scatter <= 0) {
    for (i=0; i<M_BINS; i++) {
      for (j=0; j<BHER_BINS; j++) {
  double bher = BHER_MIN+j*BHER_INV_BPDEX;
  // steps[n].bher_dist_full[i*BHER_BINS+j] = bher_prob(bher+steps[n].bh_eta[i], steps[n].bh_eta[i], steps[n].smhm.bh_alpha, steps[n].smhm.bh_delta, 1.0, bh_eta_crit)*steps[n].smhm.bh_prob_norm;
  steps[n].bher_dist[i*BHER_BINS+j] = bher_prob(bher+steps[n].bh_eta[i], steps[n].bh_eta[i], steps[n].smhm.bh_alpha, steps[n].smhm.bh_delta, 1.0, bh_eta_crit)*steps[n].smhm.bh_prob_norm;
  steps[n].ledd_min[i] = BHER_MIN;
  steps[n].ledd_max[i] = BHER_MAX;
  steps[n].ledd_bpdex[i] = BHER_BPDEX;
  //exp10(bher*steps[n].smhm.bh_alpha)*exp(-exp10(bher));
      }
    }
  } else {
    // const int64_t cmin = -BHER_BINS-6*BHER_BPDEX;
    // const int64_t cmax = 6*BHER_BPDEX + BHER_BINS;
    // //    int64_t llim=cmin, ulim = cmax;
    // float *ecache = check_realloc(NULL, sizeof(float)*(cmax-cmin+1), "Allocating ecache");


    //fprintf(stderr, "%"PRId64" %"PRId64"\n", llim, ulim);
    for (i=0; i<M_BINS; i++) {
      // float *bher_dist = steps[n].bher_dist_full+i*BHER_BINS;
      double ledd_min = steps[n].bh_eta[i] + BHER_MIN;
      if (nonlinear_luminosity && ledd_min < bh_eta_crit)
  ledd_min = (ledd_min - 0.5 * bh_eta_crit)*2.0;
      ledd_min -= steps[n].bh_eta[i];
      double ledd_eff_min = steps[n].bh_eta[i] + BHER_EFF_MIN;
      // if (nonlinear_luminosity && ledd_min < -2.0) // This is what Peter has written, which is a bug...
      if (nonlinear_luminosity && ledd_eff_min < bh_eta_crit)
  ledd_eff_min = (ledd_eff_min - 0.5 * bh_eta_crit)*2.0;
      ledd_eff_min -= steps[n].bh_eta[i];
      double ledd_max = steps[n].bh_eta[i] + BHER_MAX;
      if (nonlinear_luminosity && ledd_max < bh_eta_crit)
  ledd_max = (ledd_max - 0.5 * bh_eta_crit)*2.0;
      ledd_max -= steps[n].bh_eta[i];


      double bpdex = ((double)BHER_BINS-1)/(ledd_max - ledd_min);
      double inv_bpdex = 1.0/bpdex;
      // int64_t kmax = 6.7*scatter*bpdex;
      // int64_t kmin = -kmax;
      // if (kmax - kmin > cmax - cmin + 1)
      // {
      //   //fprintf(stderr, "Too large scatters in BH mass!\n");
      //   for (k=0; k<BHER_BINS; k++) bher_dist[k] = 0;
      //   //INVALIDATE(fit, "Too large scatter. INVALIDATE.\n");
      //   free(ecache);
      //   return;
      // }
      // double inv_scatter = 1.0/scatter;

      steps[n].ledd_min[i] = ledd_min;
      steps[n].ledd_eff_min[i] = ledd_eff_min;
      steps[n].ledd_max[i] = ledd_max;
      steps[n].ledd_bpdex[i] = bpdex;
      //if (n == 63 && i == 19) fprintf(stderr, "ledd_min=%f, ledd_eff_min=%f, ledd_max=%f, ledd_bpdex=%f\n", ledd_min, ledd_eff_min, ledd_max, bpdex);
  //     for (k=kmin; k<kmax; k++) {
  // double db = inv_scatter*k*inv_bpdex;
  // ecache[k-kmin] = exp(-0.5*db*db);
  //     }

      int64_t jmin = (ledd_eff_min-ledd_min)*bpdex;
      int64_t jmax = (ledd_max-ledd_min)*bpdex;

      for (j=jmin; j<jmax; j++) {
  float bher = ledd_min+j*inv_bpdex;
  // int64_t emin = j+kmin;
  // // if (emin < 0) emin = 0;
  // int64_t emax = j+kmax;
  // // if (emax>BHER_BINS-1) emax = BHER_BINS-1;

  // // if ((j + kmin < 0))
  // // {
  // //   for (k=0; k<BHER_BINS; k++) bher_dist[k] = 0;
  // //   INVALIDATE(fit, "Not good j + kmin. INVALIDATE.\n");
  // //   free(ecache);
  // //   return;
  // // }

  // // float *ec = ecache - (j+kmin);
  // float *ec = ecache - emin;



  float prob = bher_prob(bher+steps[n].bh_eta[i], steps[n].bh_eta[i], steps[n].smhm.bh_alpha, steps[n].smhm.bh_delta, 1.0, bh_eta_crit);
  //if (n == 63 && i == 19)
  //fprintf(stderr, "j=%d, mbh_med=%f, bher=%f, bher_0=%f, bh_alpha=%f, bh_delta=%f, prob=%e\n", j, steps[n].log_bh_mass[i], bher+steps[n].bh_eta[i], steps[n].bh_eta[i], steps[n].smhm.bh_alpha, steps[n].smhm.bh_delta, prob);
  //if (prob < 1e-15) {
  //  if (steps[n].smhm.bh_alpha > 0 && steps[n].smhm.bh_delta > 0) break;
  //  else if (bher+steps[n].bh_eta[i] > 0) break;
  //  else continue;
  //}

  //bher_dist[j] = prob;
  steps[n].bher_dist[i*BHER_BINS+j] = prob;
  // steps[n].bher_dist_full[i*BHER_BINS+j] = prob; //Use this to calculate the active fraction later.
  steps[n].bher_dist_norm[i] += prob;
  
  // for (k=emin; k<emax; k++) {
  //   if (k < 0 || k > BHER_BINS - 1) continue;
  //   float weight = ec[k]; //exp(-0.5*db*db);
  //   bher_dist[k] += weight*prob;
  // } 
  }


  //     //Normalize;
  //     double total = 0;
  //     for (k=0; k<BHER_BINS; k++) total += steps[n].bher_dist_full[i*BHER_BINS+k];
  //     total *= inv_bpdex;
  //     if (!(total>0)) {
  // // fprintf(stderr, "%e %"PRId64" %f %f\n", total, i, steps[n].scale, steps[n].bh_eta[i]);
  //     }
  //     // assert(total > 0);
  //     if (total > 0)
  //       total = 1.0/total;
  //     for (k=0; k<BHER_BINS; k++) steps[n].bher_dist_full[i*BHER_BINS+k] *= total;
      double total = 0;
      total = steps[n].bher_dist_norm[i] * inv_bpdex;
      if (total > 0)
        total = 1.0/total;
      for (k=0; k<BHER_BINS; k++) steps[n].bher_dist[i*BHER_BINS+k] *= total; 
      steps[n].bher_dist_norm[i] = bpdex;

    }
    // free(ecache);
  }

  /*  char buffer[1024];
  sprintf(buffer, "bher_dist/dist_%f.dat\n", steps[n].scale);
  FILE *out = check_fopen(buffer, "w");
  fprintf(out, "#ER Prob.\n");
  for (i=0; i<BHER_BINS; i++) fprintf(out, "%f %e\n", BHER_MIN+i*BHER_INV_BPDEX, steps[n].bher_dist[i]);
  fclose(out);
  */
}


void calc_bh_lum_distribution_full(int n, struct smf_fit *fit)
{
  int64_t i, j, k;
  memset(steps[n].lum_dist_full, 0, sizeof(double) * MBH_BINS * LBOL_BINS);
  memset(steps[n].lum_func_full, 0, sizeof(double) * MBH_BINS * LBOL_BINS);
  double l_min = 36; double l_max = 50;
  double mbh_bpdex = MBH_BINS / (steps[n].bh_mass_max - steps[n].bh_mass_min);
  double mbh_inv_bpdex = 1.0 / mbh_bpdex;
  double sm_scatter = steps[n].smhm.scatter*steps[n].smhm.bh_gamma;
  double scatter = sqrt(sm_scatter*sm_scatter + steps[n].smhm.bh_scatter*steps[n].smhm.bh_scatter);
  double gauss_norm = 1 / sqrt(2 * M_PI) / scatter;
  
  int mbh_match = 0; 
  
  for (i=0; i<MBH_BINS; i++)
  {
    double mbh = steps[n].bh_mass_min + (i + 0.5) * mbh_inv_bpdex;
    //if (n == 134 && fabs(mbh - 8) < 0.3 && mbh_match == 0) fprintf(stderr, "#mbh: %f\n", mbh);
    // Given the current code structure, it is easier to fold in the mass-dependent duty cycle here.
    double f_mass = exp((mbh - steps[n].smhm.dc_mbh) / steps[n].smhm.dc_mbh_w);
    f_mass = f_mass / (1 + f_mass);
    double dc = steps[n].smhm.bh_duty * f_mass;
    if (dc < 1e-4) dc = 1e-4;
    double tnd = 0; //Total ND of BHs with this mass, including dormant and active BHs.
    for (j=0; j < M_BINS; j++)
    {
      double dmbh = (mbh - steps[n].log_bh_mass[j]) / scatter;
      //if (fabs(dmbh) > 5) continue;
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

      // double w_mbh = gauss_norm * exp(-0.5 * dmbh * dmbh); 

      //if (n == 134 && fabs(mbh - 8) < 0.3 && mbh_match == 0) 
      //{  
//	//fprintf(stdout, "%f ", 38.1 + mbh + steps[n].bh_eta[j] + steps[n].ledd_max[j]);
  //      fprintf(stderr, "%f ", steps[n].med_hm_at_a[j]);
    //  }
      tnd += w_mbh * steps[n].t[j];
      //if (n == 134) fprintf(stderr, "n=%d, Mh=%f, mbh=%f, w_mbh=%e\n", n, M_MIN + (j + 0.5) * INV_BPDEX, mbh, w_mbh);
      
      for (k=0; k<LBOL_BINS; k++)
      
      {
        double lbol = LBOL_MIN + (k + 0.5) * LBOL_INV_BPDEX;
        double bh_eta = lbol - 38.1 - mbh;
        double eta_frac = bh_eta - steps[n].bh_eta[j];
        //if (n == 134 && mbh == 7.715) fprintf(stderr, "n=%d, j=%d, mbh=%f, lbol=%f, bh_eta=%f, bh_eta0=%f, eta_frac=%f, ledd_max=%f\n", n, j, mbh, lbol, bh_eta, steps[n].bh_eta[j], eta_frac, steps[n].ledd_max[j]);
        if ((!isfinite(eta_frac)) || eta_frac < steps[n].ledd_min[j] || eta_frac > steps[n].ledd_max[j])
	{
     //           if (n == 134 && fabs(mbh - 8) < 0.3 && mbh_match == 0) fprintf(stderr, "0 ");
                continue;
        } 

        double bher_f = (eta_frac-steps[n].ledd_min[j])*steps[n].ledd_bpdex[j];
        int64_t bher_b = bher_f;
        bher_f -= bher_b;

        double p1 = steps[n].bher_dist[j*BHER_BINS+bher_b];
        double p2 = steps[n].bher_dist[j*BHER_BINS+bher_b+1];
        if (bher_b >= BHER_BINS-1) p2 = p1;
        // double f_CTK = steps[qi->step].f_CTK[mb];
        double nd_l = (p1 + bher_f*(p2-p1)) * w_mbh * steps[n].t[j] * dc; //The number density (Mpc^-3)
                                                                          //of AGN with this luminosity
                                                                          //and this Mbh in this Mh bin.
        //if (n == 63 && j==19 && lbol >= 46.725 && lbol <= 47.075 && mbh >= 7.075 && mbh <= 7.095)
	    //fprintf(stderr, "lbol=%f, j=%d, nd_l=%e\n", lbol, j, nd_l);
	  //  fprintf(stderr, "lbol=%f, j=%d, dmbh=%f-sigma, bher_b=%d, bher_f=%f, p1=%e, p2=%e, w_mbh=%e, dc=%e, nd=%e, nd_l=%e\n", lbol, j, dmbh, bher_b, bher_f, p1, p2, w_mbh, dc, steps[n].t[j], nd_l); 
        //if (n == 134 && fabs(mbh - 8) < 0.3 && mbh_match == 0) fprintf(stderr, "%e ", nd_l);
        //if (n == 134 && mbh == 7.715) fprintf(stderr, "n=%d, j=%d, mbh=%f, lbol=%f, bh_eta=%f, bh_eta0=%f, eta_frac=%f, ledd_max=%f, nd_l=%e, bher_b=%d, bher_f=%f, p1=%e, p2=%e\n", n, j, mbh, lbol, bh_eta, steps[n].bh_eta[j], eta_frac, steps[n].ledd_max[j], nd_l, bher_b, bher_f, p1, p2);
	steps[n].lum_func_full[i*LBOL_BINS+k] += nd_l;
        // steps[n].lum_dist_full[i*MBH_BINS+k] = steps[n].lum_func_full[i*MBH_BINS+k];
        //if (n == 134) fprintf(stderr, "n=%d, Mh=%f, mbh=%f, w_mbh=%e, lbol=%f, eta_frac=%f, bher_f=%f, bher_b=%d, p1=%e, p2=%e, dc=%e\n", n, M_MIN + (j + 0.5) * INV_BPDEX, mbh, w_mbh, lbol, eta_frac, bher_f, bher_b, p1, p2, dc);
      }
     // if (n == 63 && mbh >= 7.075 && mbh <= 7.095)
       //   fprintf(stderr, "\n");
      //if (n == 134 && fabs(mbh - 8) < 0.3 && mbh_match == 0) fprintf(stderr, "\n");
      
    }
    steps[n].bhmf[i] = tnd * mbh_bpdex;
    
    //if (n == 134 && fabs(mbh - 8) < 0.3 && mbh_match == 0) 
    //{
//	mbh_match = 1;
  //      fprintf(stdout, "\n");
    //}
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
void calc_bh_acc_rate_distribution_kinetic(int n, struct smf_fit *fit) {

  if (!nonlinear_luminosity) return; //The linear luminosity ansatz implies that all the
                                      //gravitational energy from accretion gets converted
                                      //into radiative energy, so no kinetic energy.

  int64_t i, j, k;
  memset(steps[n].bher_dist_kin, 0, sizeof(float)*BHER_BINS*M_BINS);

  double bh_eta_crit = steps[n].smhm.bh_eta_crit;
  double sm_scatter = steps[n].smhm.scatter*steps[n].smhm.bh_gamma;
  double scatter = sqrt(sm_scatter*sm_scatter + steps[n].smhm.bh_scatter*steps[n].smhm.bh_scatter);


  // Since this calculation only occurs when the luminosity is nonlinear, there is no
  // need to differentiate between linear and non-linear cases.
  double (*bher_prob)(double, double, double, double, double, double) = &_prob_of_ledd_kinetic;
  // if (nonlinear_luminosity) bher_prob = &_prob_of_ledd_nonlinear;


  if (scatter <= 0) {
    for (i=0; i<M_BINS; i++) {
      for (j=0; j<BHER_BINS; j++) {
  double bher = BHER_MIN+j*BHER_INV_BPDEX;
  steps[n].bher_dist_kin[i*BHER_BINS+j] = bher_prob(bher+steps[n].bh_eta[i], steps[n].bh_eta[i], steps[n].smhm.bh_alpha, steps[n].smhm.bh_delta, 1.0, bh_eta_crit)*steps[n].smhm.bh_prob_norm;
  steps[n].ledd_min[i] = BHER_MIN;
  steps[n].ledd_max[i] = BHER_MAX;
  steps[n].ledd_bpdex[i] = BHER_BPDEX;
  //exp10(bher*steps[n].smhm.bh_alpha)*exp(-exp10(bher));
      }
    }
  } else {
    const int64_t cmin = -BHER_BINS-6*BHER_BPDEX;
    const int64_t cmax = 6*BHER_BPDEX + BHER_BINS;
    //    int64_t llim=cmin, ulim = cmax;
    float *ecache = check_realloc(NULL, sizeof(float)*(cmax-cmin+1), "Allocating ecache");
    /*
    float *prob_dist = check_realloc(NULL, sizeof(float)*(BHER_BINS*2), "Allocating prob dist");
    double inv_scatter = 1.0/scatter;
    for (k=cmin; k<=cmax; k++) {
      double db = inv_scatter*k*BHER_INV_BPDEX;
      ecache[k-cmin] = exp(-0.5*db*db);
      if (ecache[k-cmin] < 1e-7) {
  if (llim==k-1 || llim==cmin) llim = k;
  else if (ulim == cmax) ulim = k;
      }
    }
    */

    //fprintf(stderr, "%"PRId64" %"PRId64"\n", llim, ulim);
    for (i=0; i<M_BINS; i++) {
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
        //fprintf(stderr, "Too large scatters in BH mass!\n");
        for (k=0; k<BHER_BINS; k++) bher_dist[k] = 0;
        //INVALIDATE(fit, "Too large scatter. INVALIDATE.\n");
        free(ecache);
        return;
      }
      double inv_scatter = 1.0/scatter;

      // steps[n].ledd_min[i] = ledd_min;
      // steps[n].ledd_max[i] = ledd_max;
      // steps[n].ledd_bpdex[i] = bpdex;

      for (k=kmin; k<kmax; k++) {
  double db = inv_scatter*k*inv_bpdex;
  ecache[k-kmin] = exp(-0.5*db*db);
      }

      int64_t jmin = (ledd_eff_min-ledd_min)*bpdex;
      int64_t jmax = (ledd_max-ledd_min)*bpdex;

      for (j=jmin; j<jmax; j++) {
  float bher = ledd_min+j*inv_bpdex;
  int64_t emin = j+kmin;
  // if (emin < 0) emin = 0;
  int64_t emax = j+kmax;
  // if (emax>BHER_BINS-1) emax = BHER_BINS-1;

  // if ((j + kmin < 0))
  // {
  //   for (k=0; k<BHER_BINS; k++) bher_dist[k] = 0;
  //   INVALIDATE(fit, "Not good j + kmin. INVALIDATE.\n");
  //   free(ecache);
  //   return;
  // }

  // float *ec = ecache - (j+kmin);
  float *ec = ecache - emin;



  float prob = bher_prob(bher+steps[n].bh_eta[i], steps[n].bh_eta[i], steps[n].smhm.bh_alpha, steps[n].smhm.bh_delta, 1.0, bh_eta_crit);
  //if (steps[n].scale > 0.5 && steps[n].scale < 0.8) 
  //  fprintf(stderr, "n=%d, scale=%f, Mh=%f, bh_eta_0=%f, bh_eta=%f, prob=%e\n", n, steps[n].scale, M_MIN + (i + 0.5) * INV_BPDEX,
  //                    steps[n].bh_eta[i], bher+steps[n].bh_eta[i], prob);
  if (prob < 1e-15) {
    // if (steps[n].smhm.bh_alpha > 0 && steps[n].smhm.bh_delta > 0) break;
    // else if (bher+steps[n].bh_eta[i] > 0) break;
    // else continue;
    continue;
  }

  //bher_dist[j] = prob;

  
  for (k=emin; k<emax; k++) {
    if (k < 0 || k > BHER_BINS - 1) continue;
    float weight = ec[k]; //exp(-0.5*db*db);
    bher_dist[k] += weight*prob;
  } 
      }

      //Normalize;
      double total = 0;
      for (k=0; k<BHER_BINS; k++) total += bher_dist[k];
      total *= inv_bpdex;
      if (!(total>0)) {
  // fprintf(stderr, "%e %"PRId64" %f %f\n", total, i, steps[n].scale, steps[n].bh_eta[i]);
      }
      // assert(total > 0);
      if (total > 0)
  total = 1.0/total;
      for (k=0; k<BHER_BINS; k++) bher_dist[k] *= total;
    }
    free(ecache);
  }

  // for (i=0; i<M_BINS; i++)
  //   for (j=0; j<BHER_BINS; j++)
  //     fprintf(stderr, "n=%d, i=%d, j=%d, bher_dist_kin[%d]=%e\n", n, i, j, j, steps[n].bher_dist_kin[i*BHER_BINS + j]);

  /*  char buffer[1024];
  sprintf(buffer, "bher_dist/dist_%f.dat\n", steps[n].scale);
  FILE *out = check_fopen(buffer, "w");
  fprintf(out, "#ER Prob.\n");
  for (i=0; i<BHER_BINS; i++) fprintf(out, "%f %e\n", BHER_MIN+i*BHER_INV_BPDEX, steps[n].bher_dist[i]);
  fclose(out);
  */
}

void calc_avg_eta_rad(int n)
{
  int64_t i, j;
  for (i = 0; i < M_BINS; i++)
  {
    double eta0 = steps[n].bh_eta[i];
    double ledd_min = steps[n].ledd_min[i];
    double ledd_bpdex = steps[n].ledd_bpdex[i];
    double ledd_inv_bpdex = 1 / ledd_bpdex;
    double norm = 0; 
    double tot = 0;
    double dc = steps[n].smhm.bh_duty;
    //double dc = steps[n].smhm.bh_duty;
    double f_mass = exp((log10(steps[n].bh_mass_avg[i]) - steps[n].smhm.dc_mbh) / steps[n].smhm.dc_mbh_w);
    f_mass = f_mass / (1 + f_mass);
    f_mass = f_mass < 1? f_mass : 1;
    dc *= f_mass;
    if (dc < 1e-4) dc = 1e-4; 
    for (j = 0; j < BHER_BINS; j++)
    {
      double prob = steps[n].bher_dist[i*BHER_BINS + j];
      if isfinite(prob)
      {  
          //norm += prob;
          tot += prob * exp10(eta0 + ledd_min + (j + 0.5) * ledd_inv_bpdex) * ledd_inv_bpdex * dc;
      }
      }
      //if (norm) steps[n].bh_eta_rad_avg[i] = tot / norm;
      if (isfinite(tot)) steps[n].bh_eta_rad_avg[i] = tot;
      //fprintf(stderr, "n=%d, z=%f, i=%d, eta0=%e, eta_avg=%e\n", n, 1 / steps[n].scale - 1, i, eta0, steps[n].bh_eta_rad_avg[i]);
  }
}


// void calc_smf_and_ssfr(int n, struct smf_fit *fit) {
//   int64_t i, j;
//   char buffer[1024];
//   double m, m2, sm, sm_max=0, sm_min=1000;
//   int64_t bin_min=-1, bin_max=-1, bin_peak = 0;
//   gsl_interp_accel ga = {0};
//   steps[n].smhm.sm_max = 0;
//   int64_t count_falling = 0;
//   steps[n].alloc2smf = 0;
//   if (INVALID(*fit)) 
//   {
// 	  // printf("invalid!\n");
// 	  return;
//   }
//   //if (steps[n].scale < 0.1) return;

//   for (i=0; i<M_BINS; i++) {
    
     
//     // INPLEMENTAION 2: bin_max - bin_min + 1 < M_BINS

//     if (steps[n].log_sm[i]>steps[n].smhm.sm_max) {
//       sm_max = steps[n].smhm.sm_max = steps[n].log_sm[i];
//       bin_max = i;
//     }
//     // if (steps[n].log_sm[i]>-1000 && bin_min < 0) {
//     if (steps[n].log_sm[i]>0 && bin_min < 0) {
//       bin_min = i;
//       sm_min = steps[n].smhm.sm_min = steps[n].log_sm[i];
//     }

//     if (i && (steps[n].log_sm[i] <= steps[n].log_sm[i-1])) {
//       // if (steps[n].log_sm[i]==0) { steps[n].log_sm[i]=0; } // I totally don't understand why you would do this.
//       //                                                           // In this case at the very massive end of the SMF,
//       //                                                           // the corresponding halo mass would be very difficult
//       //                                                           // to determine, and cause the crazy behavior there.
//       //                                                           // For now I just put an empty block here.
//       // if (steps[n].log_sm[i]==0 && steps[n].t[i]) 
//       if (steps[n].log_sm[i]==0) 
//         { 
//           // if (bin_min >= 0) steps[n].log_sm[i] = steps[n].log_sm[i-1]+0.001; 
//           // steps[n].log_sm[i] = steps[n].log_sm[i-1]+0.001; 
//           // if (i > bin_max) bin_max = i; //Even if we don't need this artificial massive end, we still
//           //                                                       // need to maintain the bin_max to keep the size consistency between
//           //                                                       // the interpolated data and the gsl_spline object.
//           if (M_MIN+i*INV_BPDEX < 12) 
//           {
//             ;
//             // steps[n].log_sm[i] = steps[n].log_sm[i-1]+0.001;
//           } 
//           else //at the high mass end, we do a constant SM/HM ratio extrapolation.
//           {
//             steps[n].log_sm[i] = steps[n].log_sm[i-1]+INV_BPDEX;
//           }
//           //bin_max = i;
//           if (steps[n].smhm.sm_max < steps[n].log_sm[i]) 
// 	  {
// 		  sm_max = steps[n].smhm.sm_max = steps[n].log_sm[i];
// 		  bin_max = i;
// 	  }
//         }
//       else
//       {
//         if (!count_falling)
//         {
//           count_falling = 1;
//           bin_peak = i - 1;
//         }
//       }
      
//     }

//   }


//     if (bin_min < 0 || bin_max - bin_min + 1 < 10 || (count_falling && bin_peak <= bin_min + 1))
//     {
//             //sprintf(buffer, "All the stellar masses are zero.\n");
//             //fprintf(stderr, "All the stellar masses are zero at z=%f.\n", 1.0 / steps[n].scale - 1);
//             if (n > 0) 
// 		{
// 		fprintf(stderr, "All the stellar masses are zero at z=%f.\n", 1.0 / steps[n].scale - 1);
// 		INVALIDATE(fit, buffer);
//             }
// 	return;
//     }


  
//   double *log_sm_tmp = NULL;
//   double  *hm_tmp = NULL;
//   double  *sfr_tmp = NULL;

//   if (!count_falling)
//   {
//     log_sm_tmp = malloc(sizeof(double)*(bin_max - bin_min + 1));
//     hm_tmp = malloc(sizeof(double)*(bin_max - bin_min + 1));
//     sfr_tmp = malloc(sizeof(double)*(bin_max - bin_min + 1));
//     for (j=bin_min; j <= bin_max; j++)
//     {
//       // printf("log_sm=%f\n", steps[n].log_sm[j]);
//       log_sm_tmp[j - bin_min] = steps[n].log_sm[j];
//       hm_tmp[j - bin_min] = steps[n].med_hm_at_a[j];
//       // if (n == 18) printf("hm_tmp[%d]=%f, med_hm_at_a[%d]=%f\n", j - bin_min, hm_tmp[j - bin_min], j, steps[n].med_hm_at_a[j]);
//       sfr_tmp[j - bin_min] = steps[n].sfr[j];
//     }
//     // The spline does not need to be modified. What should be changed is that the dn/dlogm
//     // should be weighted by sfrac and 1 - sfrac for smf_sf and smf_q, respectively.
//     steps[n].spline = gsl_spline_alloc(gsl_interp_cspline, bin_max - bin_min + 1);
//     steps[n].spline_sfr = gsl_spline_alloc(gsl_interp_cspline, bin_max - bin_min + 1);
//     steps[n].flag_alloc = 1;
//     gsl_spline_init(steps[n].spline, log_sm_tmp, hm_tmp, bin_max - bin_min + 1);
//     gsl_spline_init(steps[n].spline_sfr, log_sm_tmp, sfr_tmp, bin_max - bin_min + 1);
//     free(log_sm_tmp); free(hm_tmp); free(sfr_tmp);
//   }

//   else
//   {
//     log_sm_tmp = malloc(sizeof(double)*(bin_peak - bin_min + 1));
//     hm_tmp = malloc(sizeof(double)*(bin_peak - bin_min + 1));
//     sfr_tmp = malloc(sizeof(double)*(bin_peak - bin_min + 1));
//     for (j=bin_min; j <= bin_peak; j++)
//     {
//       // printf("log_sm=%f\n", steps[n].log_sm[j]);
//       log_sm_tmp[j - bin_min] = steps[n].log_sm[j];
//       hm_tmp[j - bin_min] = steps[n].med_hm_at_a[j];
//       // if (n == 18) printf("hm_tmp[%d]=%f, med_hm_at_a[%d]=%f\n", j - bin_min, hm_tmp[j - bin_min], j, steps[n].med_hm_at_a[j]);
//       sfr_tmp[j - bin_min] = steps[n].sfr[j];

//     }

//     // The spline does not need to be modified. What should be changed is that the dn/dlogm
//     // should be weighted by sfrac and 1 - sfrac for smf_sf and smf_q, respectively.
//     steps[n].spline = gsl_spline_alloc(gsl_interp_cspline, bin_peak - bin_min + 1);
//     steps[n].spline_sfr = gsl_spline_alloc(gsl_interp_cspline, bin_peak - bin_min + 1);
//     gsl_spline_init(steps[n].spline, log_sm_tmp, hm_tmp, bin_peak - bin_min + 1);
//     gsl_spline_init(steps[n].spline_sfr, log_sm_tmp, sfr_tmp, bin_peak - bin_min + 1);
//     steps[n].flag_alloc = 1;
//     free(log_sm_tmp); free(hm_tmp); free(sfr_tmp);


//     int n_seg = M_BINS - bin_peak > 2 ? M_BINS - bin_peak : 3;

//     log_sm_tmp = malloc(sizeof(double)*n_seg);
//     hm_tmp = malloc(sizeof(double)*n_seg);
//     sfr_tmp = malloc(sizeof(double)*n_seg);



//     // log_sm_tmp = (double *)realloc(log_sm_tmp, sizeof(double)*(M_BINS - bin_peak + 1));
//     // hm_tmp = (double *)realloc(hm_tmp, sizeof(double)*(M_BINS - bin_peak + 1));
//     // sfr_tmp = (double *)realloc(sfr_tmp, sizeof(double)*(M_BINS - bin_peak + 1));
//     if (M_BINS - bin_peak > 2)
//     {
//       for (j=bin_peak; j < M_BINS; j++)
//       {
//         // printf("log_sm=%f\n", steps[n].log_sm[j]);
//         log_sm_tmp[j - bin_peak] = steps[n].log_sm[M_BINS - 1 - (j - bin_peak)];
//         hm_tmp[j - bin_peak] = steps[n].med_hm_at_a[M_BINS - 1 - (j - bin_peak)];
//         // if (n == 18) printf("hm_tmp[%d]=%f, med_hm_at_a[%d]=%f\n", j - bin_min, hm_tmp[j - bin_min], j, steps[n].med_hm_at_a[j]);
//         sfr_tmp[j - bin_peak] = steps[n].sfr[M_BINS - 1 - (j - bin_peak)];
//       }
//     }

//     else
//     {
//       log_sm_tmp[2] = steps[n].log_sm[bin_peak];
//       log_sm_tmp[1] = steps[n].log_sm[bin_peak + 1];
//       log_sm_tmp[0] = 2 * log_sm_tmp[1] - log_sm_tmp[2];
//       hm_tmp[0] = 16.3; hm_tmp[1] = 16.1; hm_tmp[2] = 15.9;
//     }
    

//     // if (n_seg > M_BINS - bin_peak)
//     // {
//     //   hm_tmp[2] = 16.3;
//     //   log_sm_tmp[2] = 2 * log_sm_tmp[1] - log_sm_tmp[0];
//     //   sfr_tmp[2] = 2 * sfr_tmp[1] - sfr_tmp[0];
//     // }

//     steps[n].spline2 = gsl_spline_alloc(gsl_interp_cspline, n_seg);
//     steps[n].spline_sfr2 = gsl_spline_alloc(gsl_interp_cspline, n_seg);
//     int err_spline_init1 = gsl_spline_init(steps[n].spline2, log_sm_tmp, hm_tmp, n_seg);
//     int err_spline_init2 = gsl_spline_init(steps[n].spline_sfr2, log_sm_tmp, sfr_tmp, n_seg);
//     free(log_sm_tmp); free(hm_tmp); free(sfr_tmp);
//     if ((err_spline_init1) || (err_spline_init2))
//     {
//       //fprintf(stderr, "More than 1 turning point in the SMHM.\n");
//       if (1.0 / steps[n].scale - 1 < 5) 
// 	{
// 	INVALIDATE(fit, buffer);
//         fprintf(stderr, "More than 1 turning point in the SMHM.\n");
// 	}
// 	gsl_spline_free(steps[n].spline2); gsl_spline_free(steps[n].spline_sfr2);
//       steps[n].flag_alloc = 1;
//       return;
//     }
//     steps[n].flag_alloc = 1;
//     steps[n].alloc2smf = 1;
//     // free(log_sm_tmp); free(hm_tmp); free(sfr_tmp);
//   }

//   // free(log_sm_tmp); free(hm_tmp); free(sfr_tmp);

//   // // The spline does not need to be modified. What should be changed is that the dn/dlogm
//   // // should be weighted by sfrac and 1 - sfrac for smf_sf and smf_q, respectively.
//   // steps[n].spline = gsl_spline_alloc(gsl_interp_cspline, bin_max - bin_min + 1);
//   // steps[n].spline_sfr = gsl_spline_alloc(gsl_interp_cspline, bin_max - bin_min + 1);
//   // steps[n].flag_alloc = 1;
//   // gsl_spline_init(steps[n].spline, log_sm_tmp, hm_tmp, bin_max - bin_min + 1);
//   // gsl_spline_init(steps[n].spline_sfr, log_sm_tmp, sfr_tmp, bin_max - bin_min + 1);

//   i = M_BINS-1;
//   for (j=SM_BINS-1; j>=0; j--) 
//   {
//     sm = SM_MIN + (double)j*SM_INV_BPDEX;
//     if (sm > sm_max || sm < sm_min) 
//     {
//       steps[n].sfr_sm[j+SM_EXTRA] = 0;
//       steps[n].sfr_sm_sf[j+SM_EXTRA] = 0;
//       steps[n].sfr_sm_q[j+SM_EXTRA] = 0;
//       steps[n].smf[j+SM_EXTRA] = 0;
//       //steps[n].smf_sf[j+SM_EXTRA] = 0;
//       //steps[n].smf_q[j+SM_EXTRA] = 0;
//     }
//     else 
//     {

//       int err = gsl_spline_eval_e(steps[n].spline, sm, &ga, &m);
//       // fprintf(stdout, "n=%d, z=%f, sm=%f, interpolated m=%f, err=%d\n", n, 1.0 / steps[n].scale - 1, sm, m, err);
//       if (err || m < 0)
//       //if (err || m < M_MIN - 1.0 || m > M_MAX + 1.0)
//       {
//         // sprintf(buffer, "Error in GSL spline interpolation #2.\n");
//         fprintf(stderr, "Error in GSL spline interpolation #2. m=%f\n", m);
//         INVALIDATE(fit, buffer);
//         // free(log_sm_tmp); free(hm_tmp); free(sfr_tmp);
//         return;
//       }
//       //int64_t b = gsl_interp_accel_find(&ga, steps[n].log_sm+bin_min, (bin_max-bin_min+1), sm)+bin_min;
//       // find out the sfrac to weight the smf_sf and smf_q.
      
//       double f = (m - M_MIN) * BPDEX;
//       int64_t b;
//       if (f >= M_BINS - 1)
//       {
//         b = M_BINS - 2;
//         f = 1;
//       }
//       else
//       {
//         b = f;
//         f -= b;
//       }

      

//       double sfrac_tmp = steps[n].sfrac[b] + f * (steps[n].sfrac[b+1] - steps[n].sfrac[b]);
//       steps[n].sfrac_sm[j+SM_EXTRA] = sfrac_tmp;
//       steps[n].smf[j+SM_EXTRA] = calc_smf_at_m(m,sm,n, 1, fit,&ga);
//       gsl_interp_accel_reset(&ga);
//       err = gsl_spline_eval_e(steps[n].spline_sfr, sm, &ga, &(steps[n].sfr_sm[j + SM_EXTRA]));
      
//       steps[n].sfr_sm_q[j + SM_EXTRA] = exp10(sm) * SSFR_AVG_Q;
//       steps[n].sfr_sm_sf[j + SM_EXTRA] = (steps[n].sfr_sm[j + SM_EXTRA] - (1 - sfrac_tmp) * steps[n].sfr_sm_q[j + SM_EXTRA]) / sfrac_tmp;
//       steps[n].sfr_sm[j + SM_EXTRA] *= steps[n].smf[j+SM_EXTRA];
//       steps[n].sfr_sm_q[j + SM_EXTRA] *= steps[n].smf[j+SM_EXTRA];
//       steps[n].sfr_sm_sf[j + SM_EXTRA] *= steps[n].smf[j+SM_EXTRA];
// 	// steps[n].sfr_sm[j+SM_EXTRA] = gsl_spline_eval(steps[n].spline_sfr, sm, &ga)*steps[n].smf[j+SM_EXTRA];
//       if (err)
// 	//if (err || steps[n].sfr_sm[j + SM_EXTRA] < 0 || !(isfinite(steps[n].sfr_sm[j+SM_EXTRA])))
//       {
//         //sprintf(buffer, "Error in GSL spline interpolation #3.\n");
//         fprintf(stderr, "Error in GSL spline interpolation #3. sfr=%e\n", steps[n].sfr_sm[j+SM_EXTRA]);
//         INVALIDATE(fit, buffer);
//         // free(log_sm_tmp); free(hm_tmp); free(sfr_tmp);
//         return;
//       }

//       if (steps[n].alloc2smf)
//       {
//         gsl_interp_accel_reset(&ga);
//         int err2 = gsl_spline_eval_e(steps[n].spline2, sm, &ga, &m2);
//         if (!err2)
//         {
//           double f2 = (m2 - M_MIN) * BPDEX;
//           int64_t b2;
//           if (f2 >= M_BINS - 1)
//           {
//             b2 = M_BINS - 2;
//             f2 = 1;
//           }
//           else
//           {
//             b2 = f2;
//             f2 -= b2;
//           }

// 	  if (b2 < 0) continue;

//           // Given the second branch of the SMHM, we need to modify the sfrac_sm value
//           // s.t. it is the weighted average of the two HM bins.
//           steps[n].sfrac_sm[j+SM_EXTRA] *= steps[n].t[b] + f * (steps[n].t[b+1] - steps[n].t[b]);
//           steps[n].sfrac_sm[j+SM_EXTRA] += (steps[n].sfrac[b2] + f2 * (steps[n].sfrac[b2+1] - steps[n].sfrac[b2])) *
//                                             (steps[n].t[b2] + f2 * (steps[n].t[b2+1] - steps[n].t[b2]));
//           steps[n].sfrac_sm[j+SM_EXTRA] /= (steps[n].t[b] + f * (steps[n].t[b+1] - steps[n].t[b]) + 
//                                             steps[n].t[b2] + f2 * (steps[n].t[b2+1] - steps[n].t[b2]));
//           sfrac_tmp = steps[n].sfrac[b2] + f2 * (steps[n].sfrac[b2+1] - steps[n].sfrac[b2]);
//           double smf2 = calc_smf_at_m(m2,sm,n,2,fit,&ga);
//           steps[n].smf[j+SM_EXTRA] += smf2;
//     gsl_interp_accel_reset(&ga);
//           double sfr_sm_tmp;
//           double err_sfr = gsl_spline_eval_e(steps[n].spline_sfr2, sm, &ga, &(sfr_sm_tmp));
//           if (!err_sfr)
//           {
//             steps[n].sfr_sm[j + SM_EXTRA] += sfr_sm_tmp * smf2;
//             steps[n].sfr_sm_q[j + SM_EXTRA] += exp10(sm) * SSFR_AVG_Q * smf2;
//             steps[n].sfr_sm_sf[j + SM_EXTRA] += (sfr_sm_tmp - (1 - sfrac_tmp) * exp10(sm) * SSFR_AVG_Q) / sfrac_tmp * smf2;
//           }
//         }
//       }
//     }
//   } 

  

//   //Check status of SMF bins
//   steps[n].smf_ok[SM_BINS+SM_EXTRA-1] = 0;
//   for (j=0; j<SM_BINS-1; j++) {
//     sm = SM_MIN + (double)j*SM_INV_BPDEX;
//     double avg = 0.5*(steps[n].smf[j+SM_EXTRA-1]+steps[n].smf[j+SM_EXTRA+1]);
//     if (fabs(avg-steps[n].smf[j+SM_EXTRA]) > steps[n].smf[j+SM_EXTRA]*5e-4) {
//       steps[n].smf_ok[j+SM_EXTRA] = 0;
//     } else {
//       steps[n].smf_ok[j+SM_EXTRA] = 1;
//     }
//   }
// }


// void calc_uvlf(int n, struct smf_fit *fit) {
//   int64_t i, j;
//   char buffer[1024];
//   double m, m2, uv, uv_max=-1000, uv_min=1000;
//   int64_t bin_min=-1, bin_max=-1, bin_peak=0;
//   gsl_interp_accel ga = {0};
//   steps[n].smhm.uv_max = -1000;
//   int64_t count_falling = 0;
//   steps[n].alloc2uvlf = 0;
//   if (INVALID(*fit)) 
//   {
//     // printf("invalid!\n");
//     return;
//   }
//   //if (steps[n].scale < 0.1) return;

//   for (i=0; i<M_BINS; i++) {
    
     
//     // INPLEMENTAION 2: bin_max - bin_min + 1 < M_BINS

//     if (steps[n].obs_uv[i]>steps[n].smhm.uv_max && steps[n].obs_uv[i] < 1000) {
//       uv_max = steps[n].smhm.uv_max = steps[n].obs_uv[i];
//       bin_max = i;
//     }
//     // if (steps[n].log_sm[i]>-1000 && bin_min < 0) {
//     if (steps[n].obs_uv[i]<steps[n].smhm.uv_min) {
//       bin_min = i;
//       // if (n == 66) fprintf(stderr, "n=66, bin_min=%d\n", bin_min);
//       uv_min = steps[n].smhm.uv_min = steps[n].obs_uv[i];
//     }

//     if (i && (steps[n].obs_uv[i] >= steps[n].obs_uv[i-1])) {
// 	    //fprintf(stderr, "n=%d, z=%f, m1=%f, m2=%f, uv1=%f, uv2=%f\n", n, 1.0/steps[n].scale - 1, M_MIN + (i + 0.5) * INV_BPDEX, M_MIN + (i - 0.5) * INV_BPDEX, steps[n].obs_uv[i], steps[n].obs_uv[i - 1]);
//       // if (steps[n].log_sm[i]==0) { steps[n].log_sm[i]=0; } // I totally don't understand why you would do this.
//       //                                                           // In this case at the very massive end of the SMF,
//       //                                                           // the corresponding halo mass would be very difficult
//       //                                                           // to determine, and cause the crazy behavior there.
//       //                                                           // For now I just put an empty block here.
//       if (steps[n].obs_uv[i]==10000) 
//         { 
//           // if (bin_min >= 0) steps[n].log_sm[i] = steps[n].log_sm[i-1]+0.001; 
//           // steps[n].log_sm[i] = steps[n].log_sm[i-1]+0.001; 
//           // if (i > bin_max) bin_max = i; //Even if we don't need this artificial massive end, we still
//           //                                                       // need to maintain the bin_max to keep the size consistency between
//           //                                                       // the interpolated data and the gsl_spline object.
//           if (M_MIN+i*INV_BPDEX < 12) 
//           {
//             // steps[n].obs_uv[i] = steps[n].obs_uv[i-1]-0.001;
//             ;
//           } 
//           else //at the high mass end, we do a constant SM/HM ratio extrapolation.
//           {
//             steps[n].obs_uv[i] = steps[n].obs_uv[i-1]-INV_BPDEX;
//           }
//           if (steps[n].obs_uv[i] < steps[n].smhm.uv_min) 
// 	  {
//             uv_min = steps[n].smhm.uv_min = steps[n].obs_uv[i];
//             bin_min = i;
// 	  }
//         }
//       else {
//         if (!count_falling)
//         {
//           count_falling = 1;
//           bin_peak = i - 1;
//         }
        
//   // // sprintf(buffer, "Falling SMHM relation: SM(%f) = %f; SM(%f) = %f at scale %f (Mass at z=0: %f)\n", steps[n].med_hm_at_a[i], steps[n].log_sm[i], steps[n].med_hm_at_a[i-1], steps[n].log_sm[i-1], steps[n].scale, i*INV_BPDEX+M_MIN);       
//   // fprintf(stderr, "Falling UVHM relation: UV(%f) = %f; UV(%f) = %f at scale %f (Mass at z=0: %f)\n", steps[n].med_hm_at_a[i], steps[n].obs_uv[i], steps[n].med_hm_at_a[i-1], steps[n].obs_uv[i-1], steps[n].scale, i*INV_BPDEX+M_MIN);       
//   // INVALIDATE(fit, buffer);
//   // steps[n].smhm.valid = 0;
//   // return;
//       }
//     }

//   }
//   // fprintf(stderr, "count_falling=%d, bin_peak=%d, bin_max=%d\n", count_falling, bin_peak, bin_max); 
//  // printf("bin_min=%d, bin_max=%d\n", bin_min, bin_max);
//    //if (n == 1) 
//    //{
// //	   fprintf(stderr, "bin_min=%d, bin_max=%d\n", bin_min,  bin_max);
// //	   for (int j=0; j < M_BINS; j++) fprintf(stderr, "steps[%d].obs_uv[%d]=%f ", n, j, steps[n].obs_uv[j]);
// //	   fprintf(stderr, "\n");
	   
//   // }
//     if (bin_max < 0 || bin_min - bin_max + 1 < 10 || (count_falling && bin_peak <= bin_max + 1))
// {
//         //sprintf(buffer, "All the stellar masses are zero.\n");
//         //fprintf(stderr, "bin_max=%d, bin_min=%d\n", bin_max, bin_min);
// 	//fprintf(stderr, "All the UV magnitudes are zero at z=%f (%d-th snapshot). bin_min=%d, bin_max=%d\n", 1.0 / steps[n].scale - 1, n, bin_min, bin_max);
//         if (1.0 / steps[n].scale - 1 >= 7.5) 
// 	{
// 	fprintf(stderr, "All the UV magnitudes are zero at z=%f (%d-th snapshot). bin_min=%d, bin_max=%d\n", 1.0 / steps[n].scale - 1, n, bin_min, bin_max);
// 	INVALIDATE(fit, buffer);
//         }
// 	return;
// }

  
//   double *obs_uv_tmp = NULL;
//   // double *std_uv_tmp = NULL;
//   double  *hm_tmp = NULL;

//   if (!count_falling)
//   {
//     obs_uv_tmp = malloc(sizeof(double)*(bin_min - bin_max + 1));
//     // std_uv_tmp = malloc(sizeof(double)*(bin_min - bin_max + 1));
//     hm_tmp = malloc(sizeof(double)*(bin_min - bin_max + 1));

//     for (j=bin_max; j <= bin_min; j++)
//     {
//       // printf("log_sm=%f\n", steps[n].log_sm[j]);
//       obs_uv_tmp[j - bin_max] = steps[n].obs_uv[bin_min - (j - bin_max)];
//       // std_uv_tmp[j - bin_max] = steps[n].std_uv[bin_min - (j - bin_max)];
//       hm_tmp[j - bin_max] = steps[n].med_hm_at_a[bin_min - (j - bin_max)];

//     }

    
//     // The spline does not need to be modified. What should be changed is that the dn/dlogm
//     // should be weighted by sfrac and 1 - sfrac for smf_sf and smf_q, respectively.
//     steps[n].spline_uv = gsl_spline_alloc(gsl_interp_cspline, bin_min - bin_max + 1);
//     // steps[n].spline_std_uv = gsl_spline_alloc(gsl_interp_cspline, bin_max - bin_min + 1);



//     steps[n].flag_alloc = 2;
//     gsl_spline_init(steps[n].spline_uv, obs_uv_tmp, hm_tmp, bin_min - bin_max + 1);
//     // gsl_spline_init(steps[n].spline_std_uv, obs_uv_tmp, std_uv_tmp, bin_max - bin_min + 1);
//     free(obs_uv_tmp); free(hm_tmp);
//   }

//   else
//   {
//     obs_uv_tmp = malloc(sizeof(double)*(bin_peak - bin_max + 1));
//     // std_uv_tmp = malloc(sizeof(double)*(bin_min - bin_max + 1));
//     hm_tmp = malloc(sizeof(double)*(bin_peak - bin_max + 1));

//     for (j=bin_max; j <= bin_peak; j++)
//     {
//       // printf("log_sm=%f\n", steps[n].log_sm[j]);
//       obs_uv_tmp[j - bin_max] = steps[n].obs_uv[bin_peak - (j - bin_max)];
//       // std_uv_tmp[j - bin_max] = steps[n].std_uv[bin_peak - (j - bin_max)];
//       hm_tmp[j - bin_max] = steps[n].med_hm_at_a[bin_peak - (j - bin_max)];
//       // fprintf(stderr, "obs_uv[%d]=%f\n", j - bin_max, obs_uv_tmp[j - bin_max]);
//     }

    
//     // The spline does not need to be modified. What should be changed is that the dn/dlogm
//     // should be weighted by sfrac and 1 - sfrac for smf_sf and smf_q, respectively.
//     steps[n].spline_uv = gsl_spline_alloc(gsl_interp_cspline, bin_peak - bin_max + 1);
//     // steps[n].spline_std_uv = gsl_spline_alloc(gsl_interp_cspline, bin_max - bin_min + 1);
//     gsl_spline_init(steps[n].spline_uv, obs_uv_tmp, hm_tmp, bin_peak - bin_max + 1);
//     // gsl_spline_init(steps[n].spline_std_uv, obs_uv_tmp, std_uv_tmp, bin_max - bin_min + 1);
    
//     free(obs_uv_tmp); free(hm_tmp);

//     // Here we need to account for the fact that the segment of UVHM relation after
//     // the turning over point may be fewer than 2 points. In this case
//     // cubic spline interpolation won't work. And my solution is to do a linear
//     // extrapolation out to 10^16.3 Msun to get the third point.
//     int n_seg = M_BINS - bin_peak > 2 ? M_BINS - bin_peak : 3;

//     obs_uv_tmp = malloc(sizeof(double)*n_seg);
//     // std_uv_tmp = malloc(std_uv_tmp, sizeof(double)*(bin_min - bin_max + 1));
//     hm_tmp = malloc(sizeof(double)*n_seg);

//     // obs_uv_tmp = malloc(sizeof(double)*(M_BINS - bin_peak));
//     // // std_uv_tmp = malloc(std_uv_tmp, sizeof(double)*(bin_min - bin_max + 1));
//     // hm_tmp = malloc(sizeof(double)*(M_BINS - bin_peak));

//     // obs_uv_tmp = (double *)realloc(obs_uv_tmp, sizeof(double)*(M_BINS - bin_peak));
//     // // std_uv_tmp = malloc(std_uv_tmp, sizeof(double)*(bin_min - bin_max + 1));
//     // hm_tmp = (double *)realloc(hm_tmp, sizeof(double)*(M_BINS - bin_peak));

//     for (j=bin_peak; j < M_BINS; j++)
//     {
//       // printf("log_sm=%f\n", steps[n].log_sm[j]);
//       obs_uv_tmp[j - bin_peak] = steps[n].obs_uv[j];
//       // std_uv_tmp[j - bin_peak] = steps[n].std_uv[j];
//       hm_tmp[j - bin_peak] = steps[n].med_hm_at_a[j];
//       fprintf(stderr, "bin_peak=%d, bin_min=%d, obs_uv[%d]=%f\n", bin_peak, bin_min, j - bin_peak, obs_uv_tmp[j - bin_peak]);
//     }
    
    

//     // Linear extrapolation.
//     if (n_seg > M_BINS - bin_peak)
//     {
//       // fprintf(stderr, "n_seg=%d, M_BINS - bin_peak=%d\n", n_seg, M_BINS - bin_peak);
// 	    hm_tmp[2] = 16.3; obs_uv_tmp[2] = 2 * obs_uv_tmp[1] - obs_uv_tmp[0];
//     }

//     //for (j=bin_peak; j < M_BINS; j++) fprintf(stderr, "bin_peak=%d, bin_min=%d, obs_uv[%d]=%f\n", bin_peak, bin_min, j - bin_peak, obs_uv_tmp[j - bin_peak]);

//     // The spline does not need to be modified. What should be changed is that the dn/dlogm
//     // should be weighted by sfrac and 1 - sfrac for smf_sf and smf_q, respectively.
//     steps[n].spline_uv2 = gsl_spline_alloc(gsl_interp_cspline, n_seg);
//     // steps[n].spline_std_uv = gsl_spline_alloc(gsl_interp_cspline, bin_max - bin_min + 1);
//     int err_spline_init = gsl_spline_init(steps[n].spline_uv2, obs_uv_tmp, hm_tmp, n_seg);
//     // gsl_spline_init(steps[n].spline_std_uv, obs_uv_tmp, std_uv_tmp, bin_max - bin_min + 1);

    
//     // // The spline does not need to be modified. What should be changed is that the dn/dlogm
//     // // should be weighted by sfrac and 1 - sfrac for smf_sf and smf_q, respectively.
//     // steps[n].spline_uv2 = gsl_spline_alloc(gsl_interp_cspline, M_BINS - bin_peak);
//     // // steps[n].spline_std_uv = gsl_spline_alloc(gsl_interp_cspline, bin_max - bin_min + 1);
//     // int err_spline_init = gsl_spline_init(steps[n].spline_uv2, obs_uv_tmp, hm_tmp, M_BINS - bin_peak);
//     // // gsl_spline_init(steps[n].spline_std_uv, obs_uv_tmp, std_uv_tmp, bin_max - bin_min + 1);
//     //for (j=bin_peak; j<M_BINS; j++)
//       //          fprintf(stderr, "obs_uv_tmp[%d]=%f, hm_tmp[%d]=%f\n", j-bin_peak, obs_uv_tmp[j-bin_peak], j-bin_peak, hm_tmp[j-bin_peak]);
//     free(obs_uv_tmp); free(hm_tmp);
//     if (err_spline_init)
//     {
//       //fprintf(stderr, "More than 1 turning point in the UVHM at %d-th snapshot.\n", n);
//       if (1.0/steps[n].scale - 1 > 7.5) 
// 	{
// 	fprintf(stderr, "More than 1 turning point in the UVHM at %d-th snapshot.\n", n);
// 	//for (j=bin_peak; j<M_BINS; j++)
// 	//	fprintf(stderr, "obs_uv_tmp[%d]=%f, hm_tmp[%d]=%f\n", j-bin_peak, obs_uv_tmp[j-bin_peak], j-bin_peak, hm_tmp[j-bin_peak]); 
// 	INVALIDATE(fit, buffer);
// 	}
//       //free(obs_uv_tmp); free(hm_tmp);
//       gsl_spline_free(steps[n].spline_uv2);
//       steps[n].flag_alloc = 2;
//       return;
//     }
//     steps[n].alloc2uvlf = 1;
//     steps[n].flag_alloc = 2;
    
//   }

//   // free(obs_uv_tmp); free(hm_tmp);

  
//   i = M_BINS-1;
//   for (j=UV_BINS-1; j>=0; j--) {
//     uv = UV_MIN + (double)j*UV_INV_BPMAG;
//     if (uv > uv_max || uv < uv_min) {
//       steps[n].uvlf[j+UV_EXTRA] = 0;
//       steps[n].std_uvlf[j+UV_EXTRA] = 0;
//     }
//     else {

//       int err = gsl_spline_eval_e(steps[n].spline_uv, uv, &ga, &m);
//       // fprintf(stdout, "n=%d, z=%f, uv=%f, interpolated m=%f, err=%d\n", n, 1.0 / steps[n].scale - 1, uv, m, err);
//       if (err || m < 0)
//       //if (err || m < M_MIN - 1.0 || m > M_MAX + 1.0)
//       {
//         // sprintf(buffer, "Error in GSL spline interpolation #2.\n");
//         fprintf(stderr, "Error in GSL spline interpolation #2. m=%f\n", m);
//         INVALIDATE(fit, buffer);
//         // free(obs_uv_tmp); free(hm_tmp);
//         return;
//       }
//       //int64_t b = gsl_interp_accel_find(&ga, steps[n].log_sm+bin_min, (bin_max-bin_min+1), sm)+bin_min;
//       // find out the sfrac to weight the smf_sf and smf_q.
//       double f = (m - M_MIN) * BPDEX;
//       int64_t b;
//       if (f >= M_BINS - 1)
//       {
//         b = M_BINS - 2;
//         f = 1;
//       }
//       else
//       {
//         b = f;
//         f -= b;
//       }
//       steps[n].std_uvlf[j+SM_EXTRA] = steps[n].std_uv[b] + f * (steps[n].std_uv[b+1] - steps[n].std_uv[b]);
//       steps[n].uvlf[j+SM_EXTRA] = calc_uvlf_at_m(m,uv,n,1,fit,&ga);

//       if (steps[n].alloc2uvlf)
//       {
//         gsl_interp_accel_reset(&ga);
//         int err2 = gsl_spline_eval_e(steps[n].spline_uv2, uv, &ga, &m2);
//         if (!err2)
//         {
//           // double f2 = (m2 - M_MIN) * BPDEX;
//           // int64_t b2;
//           // if (f2 >= M_BINS - 1)
//           // {
//           //   b2 = M_BINS - 2;
//           //   f2 = 1;
//           // }
//           // else
//           // {
//           //   b2 = f2;
//           //   f2 -= b2;
//           // }
//           // // I think for the second branch (i.e., massive end) of the UV-HM relation, the scatter there would
//           // // be much smaller than the one from the less massive end. So here I just ignore them.
//           // steps[n].std_uvlf[j+SM_EXTRA] = steps[n].std_uv[b] + f * (steps[n].std_uv[b+1] - steps[n].std_uv[b]);
//           steps[n].uvlf[j+SM_EXTRA] += calc_uvlf_at_m(m2,uv,n,2,fit,&ga);
//         }
//       }
      
//     }
//   }

//   // if (n == 18)
//   // {
//   //   printf("#sm dNdlgsm\n");
//   //   printf("#z=%f\n", 1 / steps[n].scale - 1);
//   //   double sm_tmp = steps[n].smhm.sm_min;
//   //   for (j = 0; sm <= steps[n].smhm.sm_max; j++)
//   //   {
//   //     m = gsl_spline_eval(steps[n].spline, sm, &ga);
//   //     printf("%f, %f\n", sm, log10(calc_smf_at_m(m,sm,n,fit,&ga)));
//   //     sm += 0.1;
//   //   }
//   // }

  

//   //Check status of SMF bins
//   steps[n].uvlf_ok[UV_BINS+UV_EXTRA-1] = 0;
//   for (j=0; j<UV_BINS-1; j++) {
//     uv = UV_MIN + (double)j*UV_INV_BPMAG;
//     double avg = 0.5*(steps[n].uvlf[j+UV_EXTRA-1]+steps[n].uvlf[j+UV_EXTRA+1]);
//     //if (n == 37) fprintf(stderr, "uv=%f, uvlf[%d]=%e\n", uv, j, steps[n].uvlf[j+UV_EXTRA]);
//     if (fabs(avg-steps[n].uvlf[j+UV_EXTRA]) > steps[n].uvlf[j+UV_EXTRA]*5e-4) {
//       steps[n].uvlf_ok[j+UV_EXTRA] = 0;
//     } else {
//       steps[n].uvlf_ok[j+UV_EXTRA] = 1;
//     }
//   }
// }


void calc_smf_and_ssfr(int n, struct smf_fit *fit) {
  int64_t i, j;
  char buffer[1024];
  double m, m2, sm, sm_max=0, sm_min=1000;
  int64_t bin_start=-1, bin_end=-1, bin_peak = 0;
  gsl_interp_accel ga = {0};
  steps[n].smhm.sm_max = 0;
  int64_t count_falling = 0;
  steps[n].alloc2smf = 0;
  if (INVALID(*fit)) 
  {
    // printf("invalid!\n");
    return;
  }
  //if (steps[n].scale < 0.1) return;

  for (i=0; i<M_BINS; i++) {

    // if (steps[n].log_sm[i]>-1000 && bin_min < 0) {
    if (steps[n].log_sm[i]>0) {
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
    
     
    // INPLEMENTAION 2: bin_max - bin_min + 1 < M_BINS

    if (steps[n].log_sm[i]>steps[n].smhm.sm_max) {
      sm_max = steps[n].smhm.sm_max = steps[n].log_sm[i];
      // bin_max = i;
    }

    if (steps[n].log_sm[i]>0 && steps[n].log_sm[i]<sm_min) {
      sm_min = steps[n].smhm.sm_min = steps[n].log_sm[i];
    }

    

    

    if (i && (steps[n].log_sm[i] <= steps[n].log_sm[i-1]) && steps[n].log_sm[i] > 0 && (!count_falling)) 
    {
          count_falling = 1;
          bin_peak = i - 1;
    }

  }


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


  
  double *log_sm_tmp = NULL;
  double  *hm_tmp = NULL;
  double  *sfr_tmp = NULL;

  if (!count_falling)
  {
    log_sm_tmp = malloc(sizeof(double)*(bin_end - bin_start + 1));
    hm_tmp = malloc(sizeof(double)*(bin_end - bin_start + 1));
    sfr_tmp = malloc(sizeof(double)*(bin_end - bin_start + 1));
    for (j=bin_start; j <= bin_end; j++)
    {
      // printf("log_sm=%f\n", steps[n].log_sm[j]);
      log_sm_tmp[j - bin_start] = steps[n].log_sm[j];
      hm_tmp[j - bin_start] = steps[n].med_hm_at_a[j];
      // if (n == 18) printf("hm_tmp[%d]=%f, med_hm_at_a[%d]=%f\n", j - bin_start, hm_tmp[j - bin_start], j, steps[n].med_hm_at_a[j]);
      sfr_tmp[j - bin_start] = steps[n].sfr[j];
    }
    // The spline does not need to be modified. What should be changed is that the dn/dlogm
    // should be weighted by sfrac and 1 - sfrac for smf_sf and smf_q, respectively.
    steps[n].spline = gsl_spline_alloc(gsl_interp_linear, bin_end - bin_start + 1);
    steps[n].spline_sfr = gsl_spline_alloc(gsl_interp_linear, bin_end - bin_start + 1);
    steps[n].flag_alloc = 1;
    gsl_spline_init(steps[n].spline, log_sm_tmp, hm_tmp, bin_end - bin_start + 1);
    gsl_spline_init(steps[n].spline_sfr, log_sm_tmp, sfr_tmp, bin_end - bin_start + 1);
    free(log_sm_tmp); free(hm_tmp); free(sfr_tmp);
  }

  else
  {
    log_sm_tmp = malloc(sizeof(double)*(bin_peak - bin_start + 1));
    hm_tmp = malloc(sizeof(double)*(bin_peak - bin_start + 1));
    sfr_tmp = malloc(sizeof(double)*(bin_peak - bin_start + 1));
    for (j=bin_start; j <= bin_peak; j++)
    {
      // printf("log_sm=%f\n", steps[n].log_sm[j]);
      log_sm_tmp[j - bin_start] = steps[n].log_sm[j];
      hm_tmp[j - bin_start] = steps[n].med_hm_at_a[j];
      // if (n == 18) printf("hm_tmp[%d]=%f, med_hm_at_a[%d]=%f\n", j - bin_start, hm_tmp[j - bin_start], j, steps[n].med_hm_at_a[j]);
      sfr_tmp[j - bin_start] = steps[n].sfr[j];

    }

    // The spline does not need to be modified. What should be changed is that the dn/dlogm
    // should be weighted by sfrac and 1 - sfrac for smf_sf and smf_q, respectively.
    steps[n].spline = gsl_spline_alloc(gsl_interp_linear, bin_peak - bin_start + 1);
    steps[n].spline_sfr = gsl_spline_alloc(gsl_interp_linear, bin_peak - bin_start + 1);
    
    if ((!steps[n].spline) || (!steps[n].spline_sfr))
    {
      //fprintf(stderr, "Too few data points to do akima interpolation for the SMHM, segment 1. scale=%f\n", steps[n].scale);
      INVALIDATE(fit, buffer);
      return;
    }


    gsl_spline_init(steps[n].spline, log_sm_tmp, hm_tmp, bin_peak - bin_start + 1);
    gsl_spline_init(steps[n].spline_sfr, log_sm_tmp, sfr_tmp, bin_peak - bin_start + 1);
    steps[n].flag_alloc = 1;
    free(log_sm_tmp); free(hm_tmp); free(sfr_tmp);


    //int n_seg = bin_end - bin_peak + 1 > 4 ? bin_end - bin_peak + 1 : 5;
    int n_seg = bin_end - bin_peak + 1;

    log_sm_tmp = malloc(sizeof(double)*n_seg);
    hm_tmp = malloc(sizeof(double)*n_seg);
    sfr_tmp = malloc(sizeof(double)*n_seg);



    // log_sm_tmp = (double *)realloc(log_sm_tmp, sizeof(double)*(M_BINS - bin_peak + 1));
    // hm_tmp = (double *)realloc(hm_tmp, sizeof(double)*(M_BINS - bin_peak + 1));
    // sfr_tmp = (double *)realloc(sfr_tmp, sizeof(double)*(M_BINS - bin_peak + 1));
    //if (bin_end - bin_peak + 1 > 4)
    //{
      for (j=bin_peak; j <= bin_end; j++)
      {
        // printf("log_sm=%f\n", steps[n].log_sm[j]);
        log_sm_tmp[j - bin_peak] = steps[n].log_sm[bin_end - (j - bin_peak)];
        hm_tmp[j - bin_peak] = steps[n].med_hm_at_a[bin_end - (j - bin_peak)];
        // if (n == 18) printf("hm_tmp[%d]=%f, med_hm_at_a[%d]=%f\n", j - bin_start, hm_tmp[j - bin_start], j, steps[n].med_hm_at_a[j]);
        sfr_tmp[j - bin_peak] = steps[n].sfr[bin_end - (j - bin_peak)];
      }
    //}

    //else
    //{
      //int npt_real = bin_end - bin_peak + 1;
      //// reverse the order of the last a few points.
      //for (int kk = 0; kk < npt_real; kk++)
      //{
        //hm_tmp[n_seg - kk - 1] = steps[n].med_hm_at_a[bin_peak + kk];
        //log_sm_tmp[n_seg - kk - 1] = steps[n].log_sm[bin_peak + kk];
        //sfr_tmp[n_seg - kk - 1] = steps[n].sfr[bin_peak + kk];
      //}

      // // Linear interpolation beyond the end.
      //for (int kk = npt_real + 1; kk < n_seg; kk++)
      //{
        //hm_tmp[n_seg - kk] = hm_tmp[n_seg - kk + 1] + INV_BPDEX; 
        //log_sm_tmp[n_seg - kk] = 2 * log_sm_tmp[n_seg - kk + 1] - log_sm_tmp[n_seg - kk + 2];
        //sfr_tmp[n_seg - kk] = 2 * sfr_tmp[n_seg - kk + 1] - sfr_tmp[n_seg - kk + 2];
      //}
      // // log_sm_tmp[2] = steps[n].log_sm[bin_peak];
      // // log_sm_tmp[1] = steps[n].log_sm[bin_peak + 1];
      // // log_sm_tmp[0] = 2 * log_sm_tmp[1] - log_sm_tmp[2];
      // // hm_tmp[0] = 16.3; hm_tmp[1] = 16.1; hm_tmp[2] = 15.9;
    //}
    

    // if (n_seg > M_BINS - bin_peak)
    // {
    //   hm_tmp[2] = 16.3;
    //   log_sm_tmp[2] = 2 * log_sm_tmp[1] - log_sm_tmp[0];
    //   sfr_tmp[2] = 2 * sfr_tmp[1] - sfr_tmp[0];
    // }

    steps[n].spline2 = gsl_spline_alloc(gsl_interp_linear, n_seg);
    steps[n].spline_sfr2 = gsl_spline_alloc(gsl_interp_linear, n_seg);
    int err_spline_init1 = gsl_spline_init(steps[n].spline2, log_sm_tmp, hm_tmp, n_seg);
    int err_spline_init2 = gsl_spline_init(steps[n].spline_sfr2, log_sm_tmp, sfr_tmp, n_seg);
    free(log_sm_tmp); free(hm_tmp); free(sfr_tmp);
    if ((err_spline_init1) || (err_spline_init2))
    {
      //fprintf(stderr, "More than 1 turning point in the SMHM.\n");
      if (1.0 / steps[n].scale - 1 < 5) 
  {
         INVALIDATE(fit, buffer);
        //fprintf(stderr, "More than 1 turning point in the SMHM. at the %d-th snapshot.\n", n);
        //for (int ii=0; ii<M_BINS; ii++) fprintf(stderr, "bin_real=%d, log_sm[%d]=%f.\n", steps[n].smhm.bin_real, ii, steps[n].log_sm[ii]);
  }
  gsl_spline_free(steps[n].spline2); gsl_spline_free(steps[n].spline_sfr2);
      steps[n].flag_alloc = 1;
      return;
    }
    steps[n].flag_alloc = 1;
    steps[n].alloc2smf = 1;
    // free(log_sm_tmp); free(hm_tmp); free(sfr_tmp);
  }

  // free(log_sm_tmp); free(hm_tmp); free(sfr_tmp);

  // // The spline does not need to be modified. What should be changed is that the dn/dlogm
  // // should be weighted by sfrac and 1 - sfrac for smf_sf and smf_q, respectively.
  // steps[n].spline = gsl_spline_alloc(gsl_interp_cspline, bin_end - bin_start + 1);
  // steps[n].spline_sfr = gsl_spline_alloc(gsl_interp_cspline, bin_end - bin_start + 1);
  // steps[n].flag_alloc = 1;
  // gsl_spline_init(steps[n].spline, log_sm_tmp, hm_tmp, bin_end - bin_start + 1);
  // gsl_spline_init(steps[n].spline_sfr, log_sm_tmp, sfr_tmp, bin_end - bin_start + 1);

  i = M_BINS-1;
  for (j=SM_BINS-1; j>=0; j--) 
  {
    sm = SM_MIN + (double)j*SM_INV_BPDEX;
    if (sm > sm_max || sm < sm_min) 
    {
      steps[n].sfr_sm[j+SM_EXTRA] = 0;
      steps[n].sfr_sm_sf[j+SM_EXTRA] = 0;
      steps[n].sfr_sm_q[j+SM_EXTRA] = 0;
      steps[n].smf[j+SM_EXTRA] = 0;
      //steps[n].smf_sf[j+SM_EXTRA] = 0;
      //steps[n].smf_q[j+SM_EXTRA] = 0;
    }
    else 
    {
      gsl_interp_accel_reset(&ga);
      int err = gsl_spline_eval_e(steps[n].spline, sm, &ga, &m);
      //fprintf(stderr, "%d %f %f\n", n, sm, m);
      if (err || !isfinite(m))
      //if (err || m < M_MIN - 1.0 || m > M_MAX + 1.0)
      {
        // sprintf(buffer, "Error in GSL spline interpolation #2.\n");
        //fprintf(stderr, "Error in GSL spline interpolation #2. n=%d, sm=%f, sm_min=%f, sm_max=%f, m=%f\n", n, sm, sm_min, sm_max, m);
        //for (int ii=0; ii<M_BINS; ii++) fprintf(stderr, "log_sm[%d]=%f\n", ii, steps[n].log_sm[ii]);
        INVALIDATE(fit, buffer);
        // free(log_sm_tmp); free(hm_tmp); free(sfr_tmp);
        return;
      }
      //int64_t b = gsl_interp_accel_find(&ga, steps[n].log_sm+bin_start, (bin_end-bin_start+1), sm)+bin_start;
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
      steps[n].smf[j+SM_EXTRA] = calc_smf_at_m(m,sm,n, 1, fit,&ga);
      gsl_interp_accel_reset(&ga);
      err = gsl_spline_eval_e(steps[n].spline_sfr, sm, &ga, &(steps[n].sfr_sm[j + SM_EXTRA]));
      
      steps[n].sfr_sm_q[j + SM_EXTRA] = exp10(sm) * SSFR_AVG_Q;
      steps[n].sfr_sm_sf[j + SM_EXTRA] = (steps[n].sfr_sm[j + SM_EXTRA] - (1 - sfrac_tmp) * steps[n].sfr_sm_q[j + SM_EXTRA]) / sfrac_tmp;
      steps[n].sfr_sm[j + SM_EXTRA] *= steps[n].smf[j+SM_EXTRA];
      steps[n].sfr_sm_q[j + SM_EXTRA] *= steps[n].smf[j+SM_EXTRA];
      steps[n].sfr_sm_sf[j + SM_EXTRA] *= steps[n].smf[j+SM_EXTRA];
  // steps[n].sfr_sm[j+SM_EXTRA] = gsl_spline_eval(steps[n].spline_sfr, sm, &ga)*steps[n].smf[j+SM_EXTRA];
      if (err || !isfinite(steps[n].sfr_sm[j + SM_EXTRA]))
  //if (err || steps[n].sfr_sm[j + SM_EXTRA] < 0 || !(isfinite(steps[n].sfr_sm[j+SM_EXTRA])))
      {
        //sprintf(buffer, "Error in GSL spline interpolation #3.\n");
        //fprintf(stderr, "Error in GSL spline interpolation #3. sfr=%e\n", steps[n].sfr_sm[j+SM_EXTRA]);
        INVALIDATE(fit, buffer);
        // free(log_sm_tmp); free(hm_tmp); free(sfr_tmp);
        return;
      }

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

  

  //Check status of SMF bins
  steps[n].smf_ok[SM_BINS+SM_EXTRA-1] = 0;
  for (j=0; j<SM_BINS-1; j++) {
    sm = SM_MIN + (double)j*SM_INV_BPDEX;
    double avg = 0.5*(steps[n].smf[j+SM_EXTRA-1]+steps[n].smf[j+SM_EXTRA+1]);
    if (fabs(avg-steps[n].smf[j+SM_EXTRA]) > steps[n].smf[j+SM_EXTRA]*5e-4) {
      steps[n].smf_ok[j+SM_EXTRA] = 0;
    } else {
      steps[n].smf_ok[j+SM_EXTRA] = 1;
    }
  }
}


void calc_uvlf(int n, struct smf_fit *fit) {
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
    // printf("invalid!\n");
    return;
  }
  //if (steps[n].scale < 0.1) return;

  for (i=0; i<M_BINS; i++) {

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


    // // if (steps[n].log_sm[i]>-1000 && bin_min < 0) {
    // if (steps[n].obs_uv[i]<steps[n].smhm.uv_min) {
    //   bin_min = i;
    //   // if (n == 66) fprintf(stderr, "n=66, bin_min=%d\n", bin_min);
    //   uv_min = steps[n].smhm.uv_min = steps[n].obs_uv[i];
    // }

    if (i && (steps[n].obs_uv[i] >= steps[n].obs_uv[i-1])
        && steps[n].obs_uv[i]<1000 && (!count_falling)) 
      {

          count_falling = 1;
          bin_peak = i - 1;
        
      }
    }
  // fprintf(stderr, "count_falling=%d, bin_peak=%d, bin_max=%d\n", count_falling, bin_peak, bin_max); 
 // printf("bin_min=%d, bin_max=%d\n", bin_min, bin_max);
   //if (n == 1) 
   //{
//     fprintf(stderr, "bin_min=%d, bin_max=%d\n", bin_min,  bin_max);
//     for (int j=0; j < M_BINS; j++) fprintf(stderr, "steps[%d].obs_uv[%d]=%f ", n, j, steps[n].obs_uv[j]);
//     fprintf(stderr, "\n");
     
  // }
    if (bin_start < 0 || bin_end - bin_start + 1 < 10 || (count_falling && bin_peak <= bin_start + 1))
{
        //sprintf(buffer, "All the stellar masses are zero.\n");
        //fprintf(stderr, "bin_max=%d, bin_min=%d\n", bin_max, bin_min);
  //fprintf(stderr, "All the UV magnitudes are zero at z=%f (%d-th snapshot). bin_min=%d, bin_max=%d\n", 1.0 / steps[n].scale - 1, n, bin_min, bin_max);
        if (1.0 / steps[n].scale - 1 >= 7.5) 
  {
  //fprintf(stderr, "All the UV magnitudes are zero at z=%f (%d-th snapshot). bin_end=%d, bin_start=%d\n", 1.0 / steps[n].scale - 1, n, bin_end, bin_start);
  INVALIDATE(fit, buffer);

        }
  return;
  
}

//if (n == 1) fprintf(stderr, "uv_min=%f, uv_max=%f\n", uv_min, uv_max);

  
  double *obs_uv_tmp = NULL;
  // double *std_uv_tmp = NULL;
  double  *hm_tmp = NULL;

  if (!count_falling)
  {
    obs_uv_tmp = malloc(sizeof(double)*(bin_end - bin_start + 1));
    // std_uv_tmp = malloc(sizeof(double)*(bin_min - bin_max + 1));
    hm_tmp = malloc(sizeof(double)*(bin_end - bin_start + 1));

    for (j=bin_start; j <= bin_end; j++)
    {
      // printf("log_sm=%f\n", steps[n].log_sm[j]);
      obs_uv_tmp[j - bin_start] = steps[n].obs_uv[bin_end - (j - bin_start)];
      // std_uv_tmp[j - bin_start] = steps[n].std_uv[bin_end - (j - bin_start)];
      hm_tmp[j - bin_start] = steps[n].med_hm_at_a[bin_end - (j - bin_start)];

    }

    
    // The spline does not need to be modified. What should be changed is that the dn/dlogm
    // should be weighted by sfrac and 1 - sfrac for smf_sf and smf_q, respectively.
    //steps[n].spline_uv = gsl_spline_alloc(gsl_interp_cspline, bin_end - bin_start + 1);
    steps[n].spline_uv = gsl_spline_alloc(gsl_interp_linear, bin_end - bin_start + 1);
    // steps[n].spline_std_uv = gsl_spline_alloc(gsl_interp_cspline, bin_start - bin_end + 1);



    steps[n].flag_alloc = 2;
    gsl_spline_init(steps[n].spline_uv, obs_uv_tmp, hm_tmp, bin_end - bin_start + 1);
    // gsl_spline_init(steps[n].spline_std_uv, obs_uv_tmp, std_uv_tmp, bin_start - bin_end + 1);
    free(obs_uv_tmp); free(hm_tmp);
  }

  else
  {
    obs_uv_tmp = malloc(sizeof(double)*(bin_peak - bin_start + 1));
    // std_uv_tmp = malloc(sizeof(double)*(bin_end - bin_start + 1));
    hm_tmp = malloc(sizeof(double)*(bin_peak - bin_start + 1));

    for (j=bin_start; j <= bin_peak; j++)
    {
      // printf("log_sm=%f\n", steps[n].log_sm[j]);
      obs_uv_tmp[j - bin_start] = steps[n].obs_uv[bin_peak - (j - bin_start)];
      // std_uv_tmp[j - bin_start] = steps[n].std_uv[bin_peak - (j - bin_start)];
      hm_tmp[j - bin_start] = steps[n].med_hm_at_a[bin_peak - (j - bin_start)];
      // fprintf(stderr, "obs_uv[%d]=%f\n", j - bin_start, obs_uv_tmp[j - bin_start]);
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

    // steps[n].spline_std_uv = gsl_spline_alloc(gsl_interp_cspline, bin_start - bin_end + 1);
    gsl_spline_init(steps[n].spline_uv, obs_uv_tmp, hm_tmp, bin_peak - bin_start + 1);
    // gsl_spline_init(steps[n].spline_std_uv, obs_uv_tmp, std_uv_tmp, bin_start - bin_end + 1);
    
    free(obs_uv_tmp); free(hm_tmp);

    // Here we need to account for the fact that the segment of UVHM relation after
    // the turning over point may be fewer than 2 points. In this case
    // cubic spline interpolation won't work. And my solution is to do a linear
    // extrapolation out to 10^16.3 Msun to get the third point.
    // int n_seg = bin_end - bin_peak + 1 > 4 ? bin_end - bin_peak + 1 : 5;
    int n_seg = bin_end - bin_peak + 1;

    obs_uv_tmp = malloc(sizeof(double)*n_seg);
    // std_uv_tmp = malloc(std_uv_tmp, sizeof(double)*(bin_end - bin_start + 1));
    hm_tmp = malloc(sizeof(double)*n_seg);

    // obs_uv_tmp = malloc(sizeof(double)*(M_BINS - bin_peak));
    // // std_uv_tmp = malloc(std_uv_tmp, sizeof(double)*(bin_end - bin_start + 1));
    // hm_tmp = malloc(sizeof(double)*(M_BINS - bin_peak));

    // obs_uv_tmp = (double *)realloc(obs_uv_tmp, sizeof(double)*(M_BINS - bin_peak));
    // // std_uv_tmp = malloc(std_uv_tmp, sizeof(double)*(bin_end - bin_start + 1));
    // hm_tmp = (double *)realloc(hm_tmp, sizeof(double)*(M_BINS - bin_peak));

    for (j=bin_peak; j <= bin_end; j++)
    {
      // printf("log_sm=%f\n", steps[n].log_sm[j]);
      obs_uv_tmp[j - bin_peak] = steps[n].obs_uv[j];
      // std_uv_tmp[j - bin_peak] = steps[n].std_uv[j];
      hm_tmp[j - bin_peak] = steps[n].med_hm_at_a[j];
      //if (n == 1) fprintf(stderr, "bin_peak=%d, bin_end=%d, obs_uv[%d]=%f, sfr[%d]=%e\n", bin_peak, bin_end, j - bin_peak, obs_uv_tmp[j - bin_peak], j - bin_peak, steps[n].sfr[j - bin_peak]);
    }

    //for (j=0;j<M_BINS;j++) if (n == 1) fprintf(stderr, "bin_real=%d, bin_peak=%d, bin_end=%d, obs_uv[%d]=%f, sfr[%d]=%e\n", steps[n].smhm.bin_real, bin_peak, bin_end, j, steps[n].obs_uv[j], j, steps[n].sfr[j]);



    // // Linear extrapolation.
    // if (n_seg > bin_end - bin_peak + 1)
    // {
      // // fprintf(stderr, "n_seg=%d, M_BINS - bin_peak=%d\n", n_seg, M_BINS - bin_peak);
      // int npt_real = bin_end - bin_peak + 1;
      // for (int kk = npt_real; kk < n_seg; kk++)
      // {
        // hm_tmp[kk] = steps[n].med_hm_at_a[bin_end] + (kk - npt_real + 1) * INV_BPDEX; 
        // obs_uv_tmp[kk] = 2 * obs_uv_tmp[kk-1] - obs_uv_tmp[kk-2];
      // }
    // }

    // The spline does not need to be modified. What should be changed is that the dn/dlogm
    // should be weighted by sfrac and 1 - sfrac for smf_sf and smf_q, respectively.
    //steps[n].spline_uv2 = gsl_spline_alloc(gsl_interp_cspline, n_seg);
    steps[n].spline_uv2 = gsl_spline_alloc(gsl_interp_linear, n_seg);
    // steps[n].spline_std_uv = gsl_spline_alloc(gsl_interp_cspline, bin_start - bin_end + 1);
    int err_spline_init = gsl_spline_init(steps[n].spline_uv2, obs_uv_tmp, hm_tmp, n_seg);
    // gsl_spline_init(steps[n].spline_std_uv, obs_uv_tmp, std_uv_tmp, bin_start - bin_end + 1);

    
    // // The spline does not need to be modified. What should be changed is that the dn/dlogm
    // // should be weighted by sfrac and 1 - sfrac for smf_sf and smf_q, respectively.
    // steps[n].spline_uv2 = gsl_spline_alloc(gsl_interp_cspline, M_BINS - bin_peak);
    // // steps[n].spline_std_uv = gsl_spline_alloc(gsl_interp_cspline, bin_start - bin_end + 1);
    // int err_spline_init = gsl_spline_init(steps[n].spline_uv2, obs_uv_tmp, hm_tmp, M_BINS - bin_peak);
    // // gsl_spline_init(steps[n].spline_std_uv, obs_uv_tmp, std_uv_tmp, bin_start - bin_end + 1);
    free(obs_uv_tmp); free(hm_tmp);
    if (err_spline_init)
    {
      //fprintf(stderr, "More than 1 turning point in the UVHM at %d-th snapshot.\n", n);
      if (1.0/steps[n].scale - 1 > 7.5) 
  {
  //fprintf(stderr, "More than 1 turning point in the UVHM at %d-th snapshot.\n", n);
  INVALIDATE(fit, buffer);
  }
      //free(obs_uv_tmp); free(hm_tmp);
      gsl_spline_free(steps[n].spline_uv2);
      steps[n].flag_alloc = 2;
      return;
    }
    steps[n].alloc2uvlf = 1;
    steps[n].flag_alloc = 2;
    
  }

  // free(obs_uv_tmp); free(hm_tmp);

  
  i = M_BINS-1;
  for (j=UV_BINS-1; j>=0; j--) {
    uv = UV_MIN + (double)j*UV_INV_BPMAG;
    if (uv > uv_max || uv < uv_min) {
      steps[n].uvlf[j+UV_EXTRA] = 0;
      steps[n].std_uvlf[j+UV_EXTRA] = 0;
    }
    else {
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
      //int64_t b = gsl_interp_accel_find(&ga, steps[n].log_sm+bin_end, (bin_start-bin_end+1), sm)+bin_end;
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
      //if (n == 32) fprintf(stderr, "uv=%f, m=%f, interpolated uvlf=%e\n", uv, m, steps[n].uvlf[j+SM_EXTRA]);
      if (steps[n].alloc2uvlf)
      {
        gsl_interp_accel_reset(&ga);
        int err2 = gsl_spline_eval_e(steps[n].spline_uv2, uv, &ga, &m2);
        if (!err2 && isfinite(m2))
        {
          // double f2 = (m2 - M_MIN) * BPDEX;
          // int64_t b2;
          // if (f2 >= M_BINS - 1)
          // {
          //   b2 = M_BINS - 2;
          //   f2 = 1;
          // }
          // else
          // {
          //   b2 = f2;
          //   f2 -= b2;
          // }
          // // I think for the second branch (i.e., massive end) of the UV-HM relation, the scatter there would
          // // be much smaller than the one from the less massive end. So here I just ignore them.
          // steps[n].std_uvlf[j+SM_EXTRA] = steps[n].std_uv[b] + f * (steps[n].std_uv[b+1] - steps[n].std_uv[b]);
          steps[n].uvlf[j+SM_EXTRA] += calc_uvlf_at_m(m2,uv,n,2,fit,&ga);
        }
      }
      
    }
  }

  // if (n == 18)
  // {
  //   printf("#sm dNdlgsm\n");
  //   printf("#z=%f\n", 1 / steps[n].scale - 1);
  //   double sm_tmp = steps[n].smhm.sm_min;
  //   for (j = 0; sm <= steps[n].smhm.sm_max; j++)
  //   {
  //     m = gsl_spline_eval(steps[n].spline, sm, &ga);
  //     printf("%f, %f\n", sm, log10(calc_smf_at_m(m,sm,n,fit,&ga)));
  //     sm += 0.1;
  //   }
  // }

  

  //Check status of SMF bins
  steps[n].uvlf_ok[UV_BINS+UV_EXTRA-1] = 0;
  for (j=0; j<UV_BINS-1; j++) {
    uv = UV_MIN + (double)j*UV_INV_BPMAG;
    double avg = 0.5*(steps[n].uvlf[j+UV_EXTRA-1]+steps[n].uvlf[j+UV_EXTRA+1]);
    //if (n == 32) fprintf(stderr, "uv=%f, uvlf[%d]=%e, std_uvlf[%d]=%f\n", uv, j, steps[n].uvlf[j+UV_EXTRA], j, steps[n].std_uvlf[j+UV_EXTRA]);
    if (fabs(avg-steps[n].uvlf[j+UV_EXTRA]) > steps[n].uvlf[j+UV_EXTRA]*5e-4) {
      steps[n].uvlf_ok[j+UV_EXTRA] = 0;
    } else {
      steps[n].uvlf_ok[j+UV_EXTRA] = 1;
    }
  }
}

void create_fake_sfr_hist(int n, int i) {
  double frac_lost = 0, weight = 0;
  double sm, total_time, time, prev_scale;
  int k;
  //if (n == 0)
  //{
//	  steps[n].sm[i] = 0;
//	  steps[n].sm_hist[i*num_outputs] = 0;
//	  return;
  //}
  for (k=0; k<=n; k++) {
    frac_lost += steps[n].smloss[k]*steps[k].dt;
    weight += steps[k].dt;
  }
  frac_lost /= weight; //Average frac lost for constant sfr
  sm = steps[n].sm[i] / frac_lost;
  total_time = scale_to_years(steps[n].scale)-scale_to_years(0);
  for (k=0; k<=n; k++) {
    prev_scale = (n) ? steps[n-1].scale : 0;
    time = scale_to_years(steps[n].scale) - scale_to_years(prev_scale);
    steps[n].sm_hist[i*num_outputs + k] = sm*time/total_time;
  }
  steps[n].new_sm[i] = steps[n].sm_hist[i*num_outputs + n];
  //if (n == 0)
  //{
        //  steps[n].sm[i] = 0;
      //    steps[n].sm_hist[i*num_outputs] = 0;
    //      return;
  //}
}

//inline 
double icl_fract(int64_t i, int64_t j, int64_t n, double ICL_RATIO) {
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


void calc_sm_hist(int n, struct smf_fit *fit) {
  int64_t i, j, k, bin;
  double t; //, ICL_FRACTION_E = ICL_FRAC_E(*fit);
  double *sm_hist_i, *sm_hist_j, sm_inc, sm_mmp;
  double ICL_RATIO = doexp10(ICL_FRAC(*fit));
  double *icl_hist_i, *icl_hist_j;
  double *sm_inc_hist_i, *sm_inc_hist_j;
  double icl_frac, ejec_frac;
  steps[n].smhm = smhm_at_z((1.0/steps[n].scale)-1.0, *fit);
  if (steps[n].smhm.icl_frac < 1e-20 || steps[n].smhm.icl_frac > 1)
  {
        //fprintf(stderr, "Bad ICL fraction: %e at z=%f\n", steps[n].smhm.icl_frac, 1/steps[n].scale-1);
        INVALIDATE(fit, "Bad ICL fraction");
  }

  // Calculate the k_uv and b_uv at the bin "bin_real" before OMP takes effect,
  // otherwise for some mass bins both of these values could be zero and cause 
  // errors.

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
  double scatter_corr = steps[n].smhm.scatter_corr; 
  

  double a200 = steps[n].scale/0.3782;
  double m200 = log10(1.115e12/0.68 / (pow(a200, -0.142441) + pow(a200, -1.78959)));
  double v200 = log10(200);
  
#pragma omp for schedule(dynamic,5) private(j,k,sm_hist_i,sm_hist_j,bin,t)
  for (i=0; i<M_BINS; i++) {
    double lv = v200 + (steps[n].med_hm_at_a[i]-m200)/3.0;
    steps[n].lv[i] = lv;
    // if (n == 150) fprintf(stderr, "Mh=%f, lv=%f\n",steps[n].med_hm_at_a[i], steps[n].lv[i]); 
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
    

    if (n == 0) steps[n].sfr[i] = 0;
     //if (i == 0) fprintf(stderr, "z=%f, sfr=%e, before dilution; ", 1 / steps[n].scale - 1, steps[n].sfr[i]);
    //if (steps[n].t[i]) steps[n].sfr[i] *= steps[n].n[i] / steps[n].t[i]; //account only for the SF from existing halos.
    //fprintf(stderr, "z=%f, sfr=%e, after dilution.\n", 1 / steps[n].scale - 1, steps[n].sfr[i]);
    //fprintf(stderr, "z=%f, Mh=%f, n/t=%e\n", 1 / steps[n].scale - 1, M_MIN + (i + 0.5) * INV_BPDEX, steps[n].n[i] / steps[n].t[i]);

    //steps[n].log_sm[i] = calc_sm_at_m(steps[n].med_hm_at_a[i], steps[n].smhm);
    //steps[n].sm[i] = doexp10(steps[n].log_sm[i]);
    //steps[n].sm_avg[i] = steps[n].sm[i]*scatter_corr;
    
    // Note that we shouldn't calculate the BH mass here, as the real stellar mass will not
    // be calculated until calc_new_sm_and_ssfr().
    steps[n].sm[i] = steps[n].sfr[i]*steps[n].dt;
    steps[n].sm_from_icl[i] = 0;
    steps[n].old_bh_mass[i] = 0;
    steps[n].bh_unmerged[i] = 0;
    // steps[n].bh_merged[i] = 0;
    // for (k=0; k < MBH_BINS; k++) steps[n].bh_unmerged_dist[i*MBH_BINS+k] = 0;
    if (no_z_scaling) continue;

    if (!steps[n].n[i]) {
      if (steps[n].c[i]) {
	create_fake_sfr_hist(n, i);
      }
    }
    else {
      sm_hist_i = &(steps[n].sm_hist[i*num_outputs]);
      memset(sm_hist_i, 0, sizeof(double)*num_outputs);
      icl_hist_i = &(steps[n].icl_stars[i*num_outputs]);
      memset(icl_hist_i, 0, sizeof(double)*num_outputs);
      sm_inc_hist_i = &(steps[n].sm_inc_hist[i*num_outputs]);
      memset(sm_inc_hist_i, 0, sizeof(double)*num_outputs);
      if (n>0) {
	sm_inc = sm_mmp = 0;
	for (j=0; j<M_BINS; j++) {
	  bin = j*M_BINS + i;
    steps[n].old_bh_mass[i] += steps[n-1].bh_mass_avg[j]*steps[n].mmp[bin];
    // if (j == i - 1) steps[n].old_bh_mass_next[i] = steps[n-1].bh_mass_avg[j]*steps[n].mmp[bin];
    // else if (j == i) steps[n].old_bh_mass_self[i] = steps[n-1].bh_mass_avg[j]*steps[n].mmp[bin];
	  steps[n].bh_unmerged[i] += steps[n].merged[bin]*(steps[n-1].bh_mass_avg[j] + 
                              steps[n-1].bh_unmerged[j]) + steps[n].mmp[bin]*steps[n-1].bh_unmerged[j];
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
    


    icl_frac = icl_fract(i,j,n,ICL_RATIO);
	  sm_inc += steps[n].merged[bin]*steps[n-1].sm_avg[j];
	  if (icl_frac > 0.999) continue;
	  sm_mmp += steps[n].mmp[bin]*steps[n-1].sm_avg[j];
	}
	steps[n].sm_inc[i] = sm_inc;
	sm_inc = 0;
	//if (sm_mmp) ejec_frac = (ICL_FRACTION_E-1.0) * sm_inc / sm_mmp;
	//else ejec_frac = 0;
	ejec_frac = 0;
	//if (ICL_FRACTION_E>1.0) steps[n].sm_from_icl[i] -= (ICL_FRACTION_E-1.0)*sm_inc;

	if (ejec_frac > 1) ejec_frac = 1;
	if (ejec_frac < 0) ejec_frac = 0;
	for (j=0; j<M_BINS; j++) {
	  bin = j*M_BINS + i;
	  icl_frac = icl_fract(i,j,n,ICL_RATIO);
	  double incoming_frac = 0; //(1.0-icl_frac)*(1.0-ICL_FRACTION_E);
	  //if (incoming_frac < 0) incoming_frac = 0;
	  t = steps[n].mmp[bin]*(1.0-ejec_frac) + 
	    incoming_frac*steps[n].merged[bin];
	  double iclt1 = (ejec_frac*steps[n].mmp[bin]+
			 (1.0-incoming_frac)*steps[n].merged[bin]);
	  //double iclt2 = (steps[n].mmp[bin]+steps[n].merged[bin]);
	  double iclt2 = (steps[n].mmp[bin]);
	  if (!t && !iclt1 && !iclt2) continue;
	  steps[n].sm_from_icl[i] += steps[n-1].sm_from_icl[j]*steps[n].mmp[bin];
	  sm_hist_j = &(steps[n-1].sm_hist[j*num_outputs]);
	  icl_hist_j = &(steps[n-1].icl_stars[j*num_outputs]);
	  for (k=0; k<n; k++) {
	    sm_hist_i[k] += t*sm_hist_j[k];
	    icl_hist_i[k] += iclt1*sm_hist_j[k];
	    icl_hist_i[k] += iclt2*icl_hist_j[k];
      // sm_inc_hist stores the ICL that come from the stars of satellite galaxies,
      // excluding the ICL from satellites.
      sm_inc_hist_i[k] += steps[n].merged[bin] * sm_hist_j[k];
	  }
	}
      }
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

    calc_old_sm(n,i);
    // correct the SFR by accounting for the contribution from quenched galaxies.
    //double ratio = steps[n].t[i] ? steps[n].n[i] / steps[n].t[i] : 1;
    steps[n].sfr[i] += (1 - steps[n].sfrac[i]) * steps[n].old_sm[i] * SSFR_AVG_Q;
    if (n == 0) steps[n].sfr[i] = 0;


    // Calculate observed UV magnitudes.
    // double z, a1, mh, k_uv, b_uv, k_std_uv, b_std_uv, a_k, b_k, c_k, a_b, b_b, c_b;
    double k_uv, b_uv, std_uv;
    // double a_k_std, b_k_std, a_b_std, b_b_std;
    // z = 1.0 / steps[n].scale - 1;
    // a1 = steps[n].scale - 1;
    // a_k = 0.15367887; b_k = -2.87608013; c_k = 9.4778494 + (-2.37851257) * a1;
    // a_b = -0.34749622; b_b = 6.85270974; c_b = -50.34457998 + 1.99341162 * a1;

    mh = steps[n].med_hm_at_a[i];
    

    k_uv = a_k * mh * mh + b_k * mh;
    b_uv = a_b * mh * mh + b_b * mh;
    k_uv = mh < - 0.5 * b_k / a_k ? -0.25 * b_k * b_k / a_k + c_k : k_uv + c_k;
    b_uv = mh < - 0.5 * b_b / a_b ? -0.25 * b_b * b_b / a_b + c_b : b_uv + c_b;

    std_uv = k_std_uv * mh + b_std_uv;

    // k_std_uv = 0.04724224 * mh + (-0.17632334) + 0.49534345 * a1;
    // b_std_uv = 0.03879428 * mh + (-0.53377296) + (-0.50621441) * a1;

    // steps[n].k_uv[i] = k_uv;
    // steps[n].b_uv[i] = b_uv;
    if (i > steps[n].smhm.bin_real)
    {
      k_uv = k_uv_real;
      b_uv = b_uv_real;
      std_uv = std_uv_real;
      // k_std_uv = k_std_uv_real;
      // b_std_uv = b_std_uv_real;
    }

    steps[n].k_uv[i] = k_uv;
    steps[n].b_uv[i] = b_uv;
    // steps[n].k_std_uv[i] = k_std_uv;
    // steps[n].b_std_uv[i] = b_std_uv;


    //k_std_uv = 0.04724224 * mh + (-0.17632334) + 0.49534345 * a1;
    //b_std_uv = 0.03879428 * mh + (-0.53377296) + (-0.50621441) * a1;


    

    double lgSFR = log10(steps[n].sfr[i]);
    steps[n].obs_uv[i] = k_uv * lgSFR + b_uv;
    //if (n == 1) fprintf(stderr, "k_uv=%f, b_uv=%f, lgSFR=%f, steps[%d].obs_uv[%d]=%f\n", k_uv, b_uv, lgSFR, n, i, steps[n].obs_uv[i]);
    steps[n].std_uv[i] = std_uv;
    // steps[n].std_uv[i] = k_std_uv * lgSFR + b_std_uv;
    if (steps[n].obs_uv[i] > 0) steps[n].obs_uv[i] = 10000;
    if (steps[n].std_uv[i] < 0) steps[n].std_uv[i] = 0.001; 
    if (steps[n].std_uv[i] > 1) steps[n].std_uv[i] = 1.000; 




    //if (i == 0) fprintf(stderr, "z=%f, sfr=%e, after adding QFGs.\n", 1 / steps[n].scale - 1, steps[n].sfr[i]);
    calc_new_sm_and_sfr(n,i,fit);
    //steps[n].sfr[i] -= (1 - steps[n].sfrac[i]) * steps[n].old_sm[i] * SSFR_AVG_Q;
    // // Moved to the calc_new_sm_and_sfr().
    // steps[n].log_bm[i] = bulge_mass(steps[n].log_sm[i]+steps[n].smhm.mu, steps[n].scale);
    // steps[n].log_bh_mass[i] = (vel_dispersion) ? 
    //     calc_bh_at_vdisp(steps[n].vdisp[i], steps[n].smhm)
    //   : calc_bh_at_bm(steps[n].log_bm[i], steps[n].smhm);
    // steps[n].bh_mass_avg[i] = doexp10(steps[n].log_bh_mass[i])*bh_scatter_corr;

    //if (n==32)
  //  {
//	fprintf(stderr, "n=32, Mh=%f, k_uv=%f, b_uv=%f, k_std_uv=%f, b_std_uv=%f, lgSFR=%f, obs_uv=%f\n", steps[n].med_hm_at_a[i], steps[n].k_uv[i], steps[n].b_uv[i], steps[n].k_std_uv[i], steps[n].b_std_uv[i], lgSFR, steps[n].obs_uv[i]);
//	fprintf(stderr, "n=31, Mh=%f, k_uv=%f, b_uv=%f, k_std_uv=%f, b_std_uv=%f, lgSFR=%f, obs_uv=%f\n\n\n", steps[n-1].med_hm_at_a[i], steps[n-1].k_uv[i], steps[n-1].b_uv[i], steps[n-1].k_std_uv[i], steps[n-1].b_std_uv[i], log10(steps[n-1].sfr[i]), steps[n-1].obs_uv[i]);

    //}

  }
}

//inline 
void calc_old_sm(int n, int j) {
  int k;
  steps[n].old_sm[j] = 0;
  for (k=0; k<n; k++)
    steps[n].old_sm[j] += steps[n].smloss[k]*steps[n].sm_hist[j*num_outputs+k];
  steps[n].sm_icl[j] = 0;
  for (k=0; k<n; k++)
    steps[n].sm_icl[j] += steps[n].smloss[k]*steps[n].icl_stars[j*num_outputs+k];
}

// Calculate the number of quasars that exist between (z_low, z_high), and
// above a certain luminosity threshold, lbol, from a survey that covers 
// frac_area of the WHOLE SKY.
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
    //fprintf(stderr, "z_min=%.6f, z_max=%.6f, Vcom=%.6e\n", z_min, z_max, V_comoving);
    // traverse every MBH bin.
    for (i=0; i<MBH_BINS; i++)
    {
      double n_qso_mbh = (1 - lbol_f) * steps[n].lum_func_full[i*LBOL_BINS+lbol_b];
      // fprintf(stderr, "n=%d, a=%.6f, i_MBH=%d, i_lbol=%d, lum_func_full=%e\n", 
      //        n, steps[n].scale,i, lbol_b, steps[n].lum_func_full[i*LBOL_BINS+lbol_b]);
      for (j=lbol_b+1; j<LBOL_BINS; j++)
      {
        n_qso_mbh += steps[n].lum_func_full[i*LBOL_BINS+j];
       // fprintf(stderr, "a=%.6f, i_MBH=%d, i_lbol=%d, lum_func_full=%e\n", 
       //       steps[n].scale,i, j, steps[n].lum_func_full[i*LBOL_BINS+j]);
      }
      // weight every MBH bin by the fraction that BHs in this bin get scattered
      // below 10^8 Msun due to virial estimate. This gives the number ***density***,
      // i.e., the Mpc^-3
      // n_qso_mbh *= frac_below8[i];
      //fprintf(stderr, "n=%d, a=%.6f, i_MBH=%d, n_qso_mbh=%e, frac_below8=%e\n", 
      //       n, steps[n].scale,i, n_qso_mbh, frac_below8[i]);
      n_qso_n += n_qso_mbh * frac_below8[i];
    }
    // Multiply the numeber density with comoving volume.
    // n_qso_n *= V_comoving;
    //fprintf(stderr, "n=%d, n_qso_n * V_comoving = %e\n", n, n_qso_n * V_comoving);
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

// Calculate the number ratio of quasars that exist between (z_low, z_high), and
// above/below a certain luminosity threshold, lbol, from a survey that covers 
// frac_area of the WHOLE SKY.
double ratio_high_z_qso(double z_low, double z_high, double lbol, double frac_area)
{
  int n, i, j;
  // Find out the relevant snapshots
  double a_max = 1.0 / (1 + z_low);
  double a_min = 1.0 / (1 + z_high);
  // double n_qso = 0;
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
    //fprintf(stderr, "z_min=%.6f, z_max=%.6f, Vcom=%.6e\n", z_min, z_max, V_comoving);
    // traverse every MBH bin.
    for (i=0; i<MBH_BINS; i++)
    {
      double n_qso_mbh = (1 - lbol_f) * steps[n].lum_func_full[i*LBOL_BINS+lbol_b];
      //fprintf(stderr, "n=%d, a=%.6f, i_MBH=%d, i_lbol=%d, lum_func_full=%e\n", 
        //      n, steps[n].scale,i, lbol_b, steps[n].lum_func_full[i*LBOL_BINS+lbol_b]);
      for (j=lbol_b+1; j<LBOL_BINS; j++)
      {
        n_qso_mbh += steps[n].lum_func_full[i*LBOL_BINS+j];
       // fprintf(stderr, "a=%.6f, i_MBH=%d, i_lbol=%d, lum_func_full=%e\n", 
         //     steps[n].scale,i, j, steps[n].lum_func_full[i*LBOL_BINS+j]);
      }
      // weight every MBH bin by the fraction that BHs in this bin get scattered
      // below 10^8 Msun due to virial estimate. This gives the number ***density***,
      // i.e., the Mpc^-3
      // n_qso_n += n_qso_mbh * frac_below8[i];
      n_qso_n_below += n_qso_mbh * frac_below8[i];
      n_qso_n_above += n_qso_mbh * (frac_below11[i] - frac_below8[i]);
    }
    // Multiply the numeber density with comoving volume.
    // n_qso_n *= V_comoving;
    // n_qso += n_qso_n * V_comoving;
    n_qso_below += n_qso_n_below * V_comoving;
    n_qso_above += n_qso_n_above * V_comoving;
    // Stop if the redshift is lower than the threshold.
    if (steps[n].scale > a_max) break;
  }
  // The calculated comoving volume is for the 4pi surface area, so we have to
  // multiply it by the fractional survey area relative to the whole sky.
  // Also, the luminosity functions are in the unit of dex^-1, we have to convert
  // them back by multiplying LBOL_INV_BPDEX.
  // n_qso *= frac_area * LBOL_INV_BPDEX;
  // return n_qso;

  double ratio = n_qso_below / n_qso_above;
  return ratio;
}

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
double recent_sfh_in_massive_halos(void) {
  int64_t n, j, count=0;
  double sfr = 0;
  for (n=0; n<num_outputs; n++) if (steps[n].scale > SFR_A_CONSTRAINT) break;
  for (; n<num_outputs; n++) {
    double z = 1.0 / steps[n].scale - 1;
    double corr = doexp10(steps[n].smhm.mu + steps[n].smhm.kappa * exp(-0.5 * (z - 2) * (z - 2)));
    for (j=(SFR_M_CONSTRAINT-M_MIN)*BPDEX; j<M_BINS; j++) {
      //if (steps[n].t[j]<REAL_ND_CUTOFF) continue;
      if (steps[n].t[j] <= 0) continue;
      sfr += steps[n].sfr[j]*corr;
      //sfr += steps[num_outputs-1].sm_hist[j*num_outputs + n]*corr/steps[n].dt;
      count++;
    }
  }
  sfr /= count;
  return sfr;
}

// Calculate the recent AGN kinetic power in massive halos.
double recent_kinetic_power_in_massive_halos(void) {
  int64_t n, j, k, count=0;
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
      //if (steps[n].t[j]<REAL_ND_CUTOFF) continue;
      if (steps[n].t[j] <= 0) continue;

      double dc = steps[n].smhm.bh_duty;
      double f_mass = exp((log10(steps[n].bh_mass_avg[j]) - steps[n].smhm.dc_mbh) / steps[n].smhm.dc_mbh_w);
      f_mass = f_mass / (1 + f_mass);
      // f_mass = f_mass < 1? f_mass : 1;
      dc *= f_mass;
      if (dc < 1e-4) dc = 1e-4;
      // double norm = 0;
      for (k=0; k<BHER_BINS; k++) 
      {
          L_tmp += exp10(38.1 + steps[n].log_bh_mass[j] + BHER_MIN + k*BHER_INV_BPDEX + steps[n].bh_eta[j])
                                            * steps[n].bher_dist_kin[j * BHER_BINS + k];
          // fprintf(stderr, "n=%d, scale=%f, Mh=%f, Mbh=%f, eta=%f, Prob(eta)=%e, L=%f, duty=%e\n", n, steps[n].scale, M_MIN + (j + 0.5) * INV_BPDEX,
          //         steps[n].log_bh_mass[j], steps[n].ledd_min[j] + k * BHER_INV_BPDEX + steps[n].bh_eta[j], 
          //         steps[n].bher_dist_kin[j * BHER_BINS + k],
          //         38.1 + steps[n].log_bh_mass[j] + steps[n].ledd_min[j] + k * BHER_INV_BPDEX + steps[n].bh_eta[j], dc);
          // norm += steps[n].bher_dist_kin[j * BHER_BINS + k];
          // fprintf(stderr, "n=%d, j=%d, k=%d, bher_dist=%f, L_tmp=%e, L=%e\n", n, j, k, 
          //                                                               log10(steps[n].bher_dist_kin[j * BHER_BINS + k]), exp10(38.1 + steps[n].log_bh_mass[j] + BHER_MIN + k*BHER_INV_BPDEX + steps[n].bh_eta[j]) * steps[n].bher_dist_kin[j * BHER_BINS + k],
          //                                                               exp10(38.1 + steps[n].log_bh_mass[j] + BHER_MIN + k*BHER_INV_BPDEX + steps[n].bh_eta[j]));
      }
      // From the calculation of bher_dist, we know that the normalization of bher_dist turns out to be bher_bpdex.
      L_tmp /= norm;
      // fprintf(stderr, "n=%d, scale=%f, Mh=%f, Mbh=%f, L=%e, norm=%e\n", n, steps[n].scale, M_MIN + (j + 0.5) * INV_BPDEX,
      //             steps[n].log_bh_mass[j], L_tmp, norm);
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
    double dc = steps[n].smhm.bh_duty;
    double f_mass = exp((log10(steps[n].bh_mass_avg[i]) - steps[n].smhm.dc_mbh) / steps[n].smhm.dc_mbh_w);
    f_mass = f_mass / (1 + f_mass);
    // f_mass = f_mass < 1? f_mass : 1;
    dc *= f_mass;
    if (dc < 1e-4) dc = 1e-4;

    for (k = 0; k < BHER_BINS; k++) 
    {
      eta_avg += exp10(bher_min + k * inv_bpdex) * steps[n].bher_dist_full[i * BHER_BINS + k];
      norm += steps[n].bher_dist_full[i * BHER_BINS + k];
      //fprintf(stderr, "n=%d, i=%d, k=%d, bher_dist_full=%e\n", n, i, k, steps[n].bher_dist_full[i * BHER_BINS + k]);
    }
    eta_avg /= norm;
    eta_avg *= dc;
    steps[n].bh_eta_avg[i] = log10(eta_avg);
    // fprintf(stderr, "n=%d, i=%d, bh_eta_avg=%f\n", n, i, steps[n].bh_eta_avg[i]); 
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
    double dc = steps[n].smhm.bh_duty;
    double f_mass = exp((log10(steps[n].bh_mass_avg[i]) - steps[n].smhm.dc_mbh) / steps[n].smhm.dc_mbh_w);
    f_mass = f_mass / (1 + f_mass);
    // f_mass = f_mass < 1? f_mass : 1;
    dc *= f_mass;
    if (dc < 1e-4) dc = 1e-4;

    for (k = 0; k < BHER_BINS; k++) 
    {
      eta_avg += exp10(bher_min + k * inv_bpdex) * steps[n].bher_dist_kin[i * BHER_BINS + k];
      norm += steps[n].bher_dist_kin[i * BHER_BINS + k];
      // fprintf(stderr, "n=%d, i=%d, eta_kin=%f, prob=%f\n", n, i, bher_min + k * inv_bpdex, steps[n].bher_dist_kin[i * BHER_BINS + k]);
    }
    // fprintf(stderr, "n=%d, i=%d, bh_eta_avg=%f, norm=%f\n", n, i, eta_avg, norm); 
    eta_avg /= norm;
    eta_avg *= dc;
    steps[n].bh_eta_kin_avg[i] = log10(eta_avg);
    
  }
}


// Calculate the fraction of kinetically powerful SMBHs in massive halos at low-z.
double recent_kinetic_frac_in_massive_halos(void) {
  int64_t n, j, k, count=0;
  double L_kin = 0;
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
      double L_tmp = 0;
      double frac = 0;
      //if (steps[n].t[j]<REAL_ND_CUTOFF) continue;
      if (steps[n].t[j] <= 0) continue;

      // Need to account for the fact that only part of the SMBHs are active, which is parametrized by the duty cycle.
      double dc = steps[n].smhm.bh_duty;
      double f_mass = exp((log10(steps[n].bh_mass_avg[j]) - steps[n].smhm.dc_mbh) / steps[n].smhm.dc_mbh_w);
      f_mass = f_mass / (1 + f_mass);
      // f_mass = f_mass < 1? f_mass : 1;
      dc *= f_mass;
      if (dc < 1e-4) dc = 1e-4;


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
        // for (k = b_eta + 1; k < BHER_BINS; k++)
        for (k = 0; k < BHER_BINS; k++) 
          {
            

            // fprintf(stderr, "n=%d, j=%d, k=%d, eta_frac=%f, bher_dist=%e\n", n, j, k, BHER_MIN + k * BHER_INV_BPDEX, 
            //                                                             steps[n].bher_dist_kin[j * BHER_BINS + k]);
            if (k < b_eta + 1) continue;
            frac += steps[n].bher_dist_kin[j * BHER_BINS + k];
            
          }
        frac /= norm;
      }


      // fprintf(stderr, "n=%d, scale=%f, Mh=%f, Mbh=%f, eta_crit=%f, f_eta=%f, b_eta=%d, frac=%e\n", n, steps[n].scale, M_MIN + (j + 0.5) * INV_BPDEX,
      //             steps[n].log_bh_mass[j], eta_crit, f_eta, b_eta, frac);


      frac_tot += frac * dc * steps[n].t[j];


      nd_tot += steps[n].t[j];
    
    }

  }

  frac_tot /= nd_tot;
  return frac_tot;
}



// // Calculate the recent AGN radiative power in massive halos.
// double recent_radiative_power_in_massive_halos(void) {
//   int64_t n, j, k, count=0;
//   double L_rad = 0;
//   double nd_tot = 0;

//   // Find the smallest scale factor that enters the window
//   for (n=0; n<num_outputs; n++) if (steps[n].scale > RAD_POWER_A_CONSTRAINT_LOW) break;
//   for (; n<num_outputs; n++) 
//   {
//     // If the scale factor is bigger than the upper limit, stop.
//     if (steps[n].scale > RAD_POWER_A_CONSTRAINT_HIGH) break;


//     // Only count the kinetic powers in massive halos.
//     for (j=(RAD_POWER_M_CONSTRAINT-M_MIN)*BPDEX; j<M_BINS; j++) 
//     {
//       double L_tmp = 0;
//       //if (steps[n].t[j]<REAL_ND_CUTOFF) continue;
//       if (steps[n].t[j] <= 0) continue;

//       double dc = steps[n].smhm.bh_duty;
//       double f_mass = exp((log10(steps[n].bh_mass_avg[j]) - steps[n].smhm.dc_mbh) / steps[n].smhm.dc_mbh_w);
//       f_mass = f_mass / (1 + f_mass);
//       // f_mass = f_mass < 1? f_mass : 1;
//       dc *= f_mass;
//       if (dc < 1e-4) dc = 1e-4;
//       double norm = 0;
//       for (k=0; k<BHER_BINS; k++) 
//         {
//           L_tmp += exp10(38.1 + steps[n].log_bh_mass[j] + steps[n].ledd_min[j] + k / steps[n].ledd_bpdex[j] + steps[n].bh_eta[j])
//                                               * steps[n].bher_dist_full[j * BHER_BINS + k];
//           // fprintf(stderr, "n=%d, scale=%f, Mh=%f, Mbh=%f, eta=%f, Prob(eta)=%e, L=%f, duty=%e\n", n, steps[n].scale, M_MIN + (j + 0.5) * INV_BPDEX,
//           //         steps[n].log_bh_mass[j], steps[n].ledd_min[j] + k / steps[n].ledd_bpdex[j] + steps[n].bh_eta[j], 
//           //         steps[n].bher_dist_full[j * BHER_BINS + k],
//           //         38.1 + steps[n].log_bh_mass[j] + steps[n].ledd_min[j] + k / steps[n].ledd_bpdex[j] + steps[n].bh_eta[j], dc);
//           norm += steps[n].bher_dist_full[j * BHER_BINS + k];
//         }
//       // From the calculation of bher_dist, we know that the normalization of bher_dist turns out to be bher_bpdex.
//       // L_rad *= dc * steps[n].t[j] / steps[n].ledd_bpdex[j];
//       L_tmp /= norm;
//       // fprintf(stderr, "n=%d, scale=%f, Mh=%f, Mbh=%f, L=%e, norm=%e\n", n, steps[n].scale, M_MIN + (j + 0.5) * INV_BPDEX,
//       //             steps[n].log_bh_mass[j], L_tmp, norm);
//       // L_tmp *= dc * steps[n].t[j];
//       L_rad += L_tmp * dc * steps[n].t[j];

//       nd_tot += steps[n].t[j];
    
//     }

//   }

//   L_rad /= nd_tot;
//   return log10(L_rad);
// }


// Calculate the recent AGN radiative power in massive halos.
double recent_radiative_power_in_massive_halos(void) {
  int64_t n, j, k, count=0;
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

      L_tmp = 5.66e46 * steps[n].smhm.bh_efficiency_rad * steps[n].bh_acc_rate[j];
      L_rad += L_tmp * steps[n].t[j];

      nd_tot += steps[n].t[j];
    
    }

  }

  L_rad /= nd_tot;
  return log10(L_rad);
}

double recent_sfh_in_massive_halos_nocorr(void) {
  int64_t n, j, count=0;
  double sfr = 0;
  for (n=0; n<num_outputs; n++) if (steps[n].scale > SFR_A_CONSTRAINT) break;
  for (; n<num_outputs; n++) {
    for (j=(SFR_M_CONSTRAINT-M_MIN)*BPDEX; j<M_BINS; j++) {
      if (steps[n].t[j] <= 0) continue;
      //if (steps[n].t[j]<REAL_ND_CUTOFF) continue;
      //sfr += steps[num_outputs-1].sm_hist[j*num_outputs + n]/steps[n].dt;
      sfr += steps[n].sfr[j];
      count++;
    }
  }
  sfr /= count;
  return sfr;
}

double recent_Micl_Mstar_ratio_in_massive_halos(void) {
  int64_t n, j, count=0;
  double Mstar_tot = 0;
  double Micl_tot = 0;
  double ND_tot = 0;
  for (n=0; n<num_outputs; n++) if (steps[n].scale > ICL_RATIO_A_CONSTRAINT) break;
  for (; n<num_outputs; n++) {
    // double corr = doexp10(steps[n].smhm.mu);
    for (j=(ICL_RATIO_M_CONSTRAINT-M_MIN)*BPDEX; j<M_BINS; j++) {
      //if (steps[n].t[j]<REAL_ND_CUTOFF) continue;
      Mstar_tot += steps[n].sm_avg[j] * steps[n].t[j];
      Micl_tot += steps[n].sm_icl[j] * steps[n].t[j];
      ND_tot += steps[n].t[j];
      // count++;
    }
  }
  Mstar_tot /= ND_tot;
  Micl_tot /= ND_tot;
  //fprintf(stderr, "Total Mstar: %e, total Micl: %e\n", Mstar_tot, Micl_tot);
  return log10(Micl_tot / Mstar_tot);
}

void calc_new_sm_and_sfr(int n, int i, struct smf_fit *fit) 
{
  double dt;
  int64_t j, k;
  char buffer[1024];
  dt = steps[n].dt;

  if (!steps[n].t[i]) {
    //steps[n].log_sm[i] = 0;
    //steps[n].new_sm[i] = 0;
    //steps[n].sm[i] = 0;
    return;
  }

  steps[n].new_sm[i] = steps[n].sfr[i] * dt * steps[n].smloss[n];
  // double sm_from_icl = steps[n].smhm.icl_frac * steps[n].sm_icl[i];
  // Instead of having a certain fraction of all the ICL merge into the central galaxy,
  // we only merge a fraction of ICL that comes from ***the stellar mass*** of incoming
  // satellite galaxies, i.e., steps[n].sm_inc[i]. Here we also fix the small bug that
  // the merged stellar mass does not shrink with time.
  double sm_from_icl = steps[n].smhm.icl_frac * steps[n].sm_inc[i] * steps[n].smloss[n];
  steps[n].mr[i] = sm_from_icl / steps[n].smloss[n] / dt;
  steps[n].new_sm[i] += sm_from_icl;
  //if (n == 0) fprintf(stderr, "after ICL addition, new_sm[%d]=%e\n", i, steps[n].new_sm[i]);
  steps[n].sm_avg[i] = steps[n].old_sm[i] + steps[n].new_sm[i];
  steps[n].sm[i] = steps[n].sm_avg[i] / steps[n].smhm.scatter_corr;
  steps[n].log_sm[i] = (steps[n].sm[i] > 1) ? log10(steps[n].sm[i]) : 0;

  steps[n].log_bm[i] = bulge_mass(steps[n].log_sm[i]+steps[n].smhm.mu, steps[n].scale);
  steps[n].log_sm_obs[i] = steps[n].log_sm[i]+steps[n].smhm.mu;



  double bh_sm_scatter = steps[n].smhm.scatter*steps[n].smhm.bh_gamma;
  double bh_vdisp_scatter = 0.3; //for the time being, use 0.3 dex as the vdisp scatter
  double combined_scatter = (vel_dispersion) ? sqrt(bh_vdisp_scatter  * bh_vdisp_scatter + steps[n].smhm.bh_scatter*steps[n].smhm.bh_scatter) : sqrt(bh_sm_scatter  * bh_sm_scatter + steps[n].smhm.bh_scatter*steps[n].smhm.bh_scatter);
  double bh_scatter_corr = exp(pow(combined_scatter*log(10), 2)/2.0);

  steps[n].log_bh_mass[i] = (vel_dispersion) ? 
      calc_bh_at_vdisp(steps[n].vdisp[i], steps[n].smhm)
    : calc_bh_at_bm(steps[n].log_bm[i], steps[n].smhm);
  steps[n].bh_mass_avg[i] = doexp10(steps[n].log_bh_mass[i])*bh_scatter_corr;


  float new_bh_mass = steps[n].bh_mass_avg[i] - steps[n].old_bh_mass[i];
  //float bh_m_dsm = (steps[n].log_sm[i] + steps[n].smhm.mu - steps[n].smhm.bh_merge) / steps[n].smhm.bh_merge_width;
  //float bh_m_dsm = ((M_MIN+(i+0.5)*INV_BPDEX) - steps[n].smhm.bh_merge) / steps[n].smhm.bh_merge_width;
  //float mfrac = 1.0 - 1.0/(1.0+exp(bh_m_dsm));
  float mfrac = n ? steps[n].smhm.f_merge_bh * sm_from_icl / steps[n].new_sm[i] : 0; // Note that the new_sm now includes the merger contribution.
                                                  // here we use this fraction as the BH merger fraction.

  //steps[n].bh_merge_rate[i] = new_bh_mass*mfrac/dt;
  //steps[n].bh_acc_rate[i] = (new_bh_mass*(1.0-mfrac))/dt;

  steps[n].new_bh_mass[i] = new_bh_mass;
  steps[n].bh_merge_rate[i] = mfrac * new_bh_mass / dt;
  steps[n].bh_acc_rate[i] = (1 - mfrac) * new_bh_mass /dt;
  //steps[n].bh_acc_rate[i] = (new_bh_mass)/dt;
  //if (new_bh_mass <= steps[n].bh_merged[i]) steps[n].bh_acc_rate[i] = new_bh_mass / dt;


  
  //  if ((!(new_bh_mass>=0) && steps[n].log_bh_mass[i] > BH_MASS_TO_REQUIRE_GROWTH)){
  // if ((!(new_bh_mass>=0) && steps[n].log_bh_mass[i] > BH_MASS_TO_REQUIRE_GROWTH) || (steps[n].bh_acc_rate[i] < 0 && steps[n].log_bh_mass[i] > BH_MASS_TO_REQUIRE_GROWTH)) {
  //   fprintf(stderr, "Negative BH accretion rate! (sm: %e, m: %e, scale: %f; old_bh: %e; new_bh: %e, unmerged BH mass: %e, new BH mass: %e!\n",
  //       sprintf(buffer, "Negative BH accretion rate! (sm: %e, m: %e, scale: %f; old_bh: %e; new_bh: %e!\n",
  //          steps[n].sm_avg[i], exp10(M_MIN+(i+0.5)/BPDEX), steps[n].scale, steps[n].old_bh_mass[i], steps[n].bh_mass_avg[i], steps[n].bh_unmerged[i], new_bh_mass);
  //              INVALIDATE(fit, buffer); 
  //  }

  //if (!(new_bh_mass>=0) && steps[n].log_bh_mass[i] > BH_MASS_TO_REQUIRE_GROWTH)
  if ((!(new_bh_mass>=0) || !(steps[n].bh_acc_rate[i]>0)) && (steps[n].scale > 0.08 && i>BPDEX && (steps[n].log_bh_mass[i] > BH_MASS_TO_REQUIRE_GROWTH || steps[n].med_hm_at_a[i] > 11)))
  //if (!(new_bh_mass>=0) && steps[n].scale > 0.08 && i < steps[n].smhm.bin_real && steps[n].med_hm_at_a[i] > 11) 
  {
    //fprintf(stderr, "Negative BH accretion rate! (sm: %e, m: %e, scale: %f; old_bh: %e; new_bh: %e!\n",
    //sprintf(buffer, "Negative BH accretion rate! (sm: %e, m: %e, scale: %f; old_bh: %e; new_bh: %e!\n",
    //steps[n].sm_avg[i], exp10(M_MIN+(i+0.5)/BPDEX), steps[n].scale, steps[n].old_bh_mass[i], steps[n].bh_mass_avg[i]);
    INVALIDATE(fit, buffer); 
  }

  if (steps[n].bh_unmerged[i] < mfrac*new_bh_mass) 
  {
    if (steps[n].scale > 0.08 && i>BPDEX && (steps[n].log_bh_mass[i] > BH_MASS_TO_REQUIRE_GROWTH || steps[n].med_hm_at_a[i] > 11))
    //if (steps[n].log_bh_mass[i] > BH_MASS_TO_REQUIRE_GROWTH && steps[n].scale > 0.08 && i>BPDEX && i < steps[n].smhm.bin_real && steps[n].med_hm_at_a[i] > 11)
    //if (steps[n].log_bh_mass[i] > BH_MASS_TO_REQUIRE_GROWTH && steps[n].scale > 0.15 && i>BPDEX && steps[n].t[i] > 1e-8) 
    {
      //fprintf(stderr, "Merging rate exceeds available mergers! (sm: %e, m: %e, scale: %f; old_bh: %e; new_bh: %e! DBH: %e; Merge: %e; M_avail: %e \n",
        // sprintf(buffer, "Merging rate exceeds available mergers! (sm: %e, m: %e, scale: %f; old_bh: %e; new_bh: %e! DBH: %e; Merge: %e; M_avail: %e \n",
      //steps[n].sm_avg[i], exp10(M_MIN+(i+0.5)/BPDEX), steps[n].scale, steps[n].old_bh_mass[i], steps[n].bh_mass_avg[i], new_bh_mass, new_bh_mass*mfrac, steps[n].bh_unmerged[i]);
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
 else 
  {
   steps[n].bh_unmerged[i] -= mfrac*new_bh_mass;
    // for (j=0; j<MBH_BINS; j++) steps[n].bh_unmerged_dist[i*MBH_BINS+j] *= steps[n].bh_unmerged[i] / (steps[n].bh_unmerged[i] + mfrac*new_bh_mass);
 }
  double bhar_tmp = steps[n].bh_acc_rate[i] > 0 ? steps[n].bh_acc_rate[i] : 1e-8;
  //if (steps[n].bh_acc_rate[i] <= 0) steps[n].bh_acc_rate[i] = 1e-8;
  steps[n].bh_eta[i] = 
    log10(//schechter_inv_avg(steps[n].smhm.bh_alpha, BHER_MIN, BHER_EFF_MAX)*
    bhar_tmp/steps[n].bh_mass_avg[i]*4.5e8*(steps[n].smhm.bh_efficiency_rad));

 


  steps[n].merged_frac[i] = 0;


  if (steps[n].smhm.icl_frac && steps[n].sm_icl[i]) 
  {
    // double frac_from_icl = steps[n].smhm.icl_frac; //Since steps[n].icl_frac is the fraction
    //                                           //of the ICL mass that got into the descendant,
    //                                           //***NOT*** the fraction of the new SM that come
    //                                           //from mergers.
    // if (steps[n].sm_inc[i]>0) 
    // {
    //   steps[n].merged_frac[i] = steps[n].sm_icl[i]*frac_from_icl/steps[n].sm_inc[i];
    //   if (steps[n].merged_frac[i]>1) steps[n].merged_frac[i] = 1;
    // }

    steps[n].merged_frac[i] = steps[n].smhm.icl_frac; //merged_frac is the fraction of the incoming satellite stellar mass
                                                      //that merge into the central galaxies, which under this parametrization
                                                      //is simply icl_frac.

    // if (frac_from_icl) 
    // {
    for (j=0; j<n; j++) 
    {
      // steps[n].sm_hist[i*num_outputs + j] += frac_from_icl*steps[n].icl_stars[i*num_outputs + j];
      // steps[n].icl_stars[i*num_outputs + j] *= (1.0-frac_from_icl);
      steps[n].sm_hist[i*num_outputs + j] += steps[n].smhm.icl_frac*steps[n].sm_inc_hist[i*num_outputs + j];
      steps[n].icl_stars[i*num_outputs + j] -= steps[n].smhm.icl_frac*steps[n].sm_inc_hist[i*num_outputs + j];
    }
    // }
    steps[n].new_sm[i] -= (steps[n].smhm.icl_frac * steps[n].sm_inc[i] * steps[n].smloss[n]);
  }


  steps[n].new_sm[i] /= steps[n].smloss[n];
  steps[n].sm_hist[i*num_outputs + n] = steps[n].new_sm[i];
  //steps[n].sfr[i] = steps[n].new_sm[i] / dt;

  if ((steps[n].t[i] > REAL_ND_CUTOFF) && (steps[n].sm_avg[i] > 0.17*exp10(steps[n].med_hm_at_a[i])) && (!no_z_scaling) && n>2) 
  {
    /*sprintf(buffer, "SM exceeds baryon fraction (sm: %e, m: %e, scale: %f!\n",
      steps[n].sm_avg[i], exp10(M_MIN+(i+0.5)*INV_BPDEX), steps[n].scale);
      INVALIDATE(fit, buffer);*/
    // Why did we comment that??? This should be an important constraint.
    //fprintf(stderr, "SM exceeds baryon fraction (sm: %e, m: %e, scale: %f!\n",
      //steps[n].sm_avg[i], exp10(M_MIN+(i+0.5)*INV_BPDEX), steps[n].scale);
    // sprintf(buffer, "SM exceeds baryon fraction (sm: %e, m: %e, scale: %f!\n",
    //   steps[n].sm_avg[i], exp10(M_MIN+(i+0.5)*INV_BPDEX), steps[n].scale);
      INVALIDATE(fit, buffer);
  }

  //if (steps[n].t[i]<REAL_ND_CUTOFF)
  //  steps[n].sfr[i] = 0;
  if ((!no_z_scaling) && (steps[n].sfr[i] < 0 || !isfinite(steps[n].sfr[i]))) 
  {
    char buffer[1024];
    // sprintf(buffer, "Negative SFR at a = %f and m = %f (ND=%g)! (ICL m=%f; ICL frac=%f)\n", steps[n].scale, 
     // i*INV_BPDEX + M_MIN, steps[n].t[i], steps[n].smhm.icl_m, steps[n].icl_frac[i]);
    //fprintf(stderr, "Negative SFR at a = %f and m = %f (ND=%g)! (ICL m=%f; ICL frac=%f)\n", steps[n].scale, 
      //i*INV_BPDEX + M_MIN, steps[n].t[i], steps[n].smhm.icl_frac, steps[n].smhm.icl_frac);
    INVALIDATE(fit, buffer);
  }
}



// void calc_new_sm_and_sfr(int n, int i, struct smf_fit *fit) 
// {
//   double dt;
//   int64_t j, k;
//   char buffer[1024];
//   dt = steps[n].dt;

//   if (!steps[n].t[i]) {
//     //steps[n].log_sm[i] = 0;
//     //steps[n].new_sm[i] = 0;
//     //steps[n].sm[i] = 0;
//     return;
//   }

//   steps[n].new_sm[i] = steps[n].sfr[i] * dt * steps[n].smloss[n];

//   //  steps[n].new_sm[i] = (steps[n].sm_avg[i] - steps[n].old_sm[i]);
//   float dm = (i+0.5)*INV_BPDEX+M_MIN-steps[n].smhm.icl_m;
//   float w = 16.0-steps[n].smhm.icl_m;
//   //if (steps[n].smhm.beta>0) { w = 4.0 / (steps[n].smhm.beta * M_LN10); }
//   if (w < 0.1) w = 0.1;
//   steps[n].icl_frac[i] = 1.0 - 1.0/(1.0+exp(4.0*dm/w));
//   if (steps[n].icl_frac[i] > 0.99) steps[n].icl_frac[i] = 0.99;
//   steps[n].new_sm[i] *= 1.0/(1.0 - steps[n].icl_frac[i]);

//   steps[n].sm_avg[i] = steps[n].old_sm[i] + steps[n].new_sm[i];
//   steps[n].sm[i] = steps[n].sm_avg[i] / steps[n].smhm.scatter_corr;
//   steps[n].log_sm[i] = (steps[n].sm[i] > 1) ? log10(steps[n].sm[i]) : 0;

//   double bh_sm_scatter = steps[n].smhm.scatter*steps[n].smhm.bh_gamma;
//   double bh_vdisp_scatter = 0.3; //for the time being, use 0.3 dex as the vdisp scatter
//   double combined_scatter = (vel_dispersion) ? sqrt(bh_vdisp_scatter  * bh_vdisp_scatter + steps[n].smhm.bh_scatter*steps[n].smhm.bh_scatter) : sqrt(bh_sm_scatter  * bh_sm_scatter + steps[n].smhm.bh_scatter*steps[n].smhm.bh_scatter);
//   double bh_scatter_corr = exp(pow(combined_scatter*log(10), 2)/2.0);

//   steps[n].log_bm[i] = bulge_mass(steps[n].log_sm[i]+steps[n].smhm.mu, steps[n].scale);
  
//   steps[n].log_bh_mass[i] = (vel_dispersion) ? 
//       calc_bh_at_vdisp(steps[n].vdisp[i], steps[n].smhm)
//     : calc_bh_at_bm(steps[n].log_bm[i], steps[n].smhm);
//   steps[n].bh_mass_avg[i] = doexp10(steps[n].log_bh_mass[i])*bh_scatter_corr;

//   steps[n].log_sm_obs[i] = steps[n].log_sm[i]+steps[n].smhm.mu;

//   float new_bh_mass = steps[n].bh_mass_avg[i] - steps[n].old_bh_mass[i];
//   //float bh_m_dsm = (steps[n].log_sm[i] + steps[n].smhm.mu - steps[n].smhm.bh_merge) / steps[n].smhm.bh_merge_width;
//   float bh_m_dsm = ((M_MIN+(i+0.5)*INV_BPDEX) - steps[n].smhm.bh_merge) / steps[n].smhm.bh_merge_width;
//   float mfrac = 1.0 - 1.0/(1.0+exp(bh_m_dsm));

//   steps[n].bh_merge_rate[i] = new_bh_mass*mfrac/dt;
//   steps[n].bh_acc_rate[i] = (new_bh_mass*(1.0-mfrac))/dt;

//   if (i >= 2 && steps[n].bh_acc_rate[i] >= steps[n].bh_acc_rate[i-1] 
//               && steps[n].bh_acc_rate[i-2] >= steps[n].bh_acc_rate[i-1]
//               && log10(steps[n].sm_avg[i]) >= 10)
//   {
//     // double slope = (log10(steps[n].bh_acc_rate[i-1]) - log10(steps[n].bh_acc_rate[i-2]))
//     //                 / (steps[n].log_sm[i-1] - steps[n].log_sm[i-2]);
//      //steps[n].bh_acc_rate[i] = exp10(log10(steps[n].bh_acc_rate[i-1]) + slope * (steps[n].log_sm[i] - steps[n].log_sm[i-1]));
//     //steps[n].bh_acc_rate[i] = steps[n].bh_acc_rate[i-1];
//     if (!isfinite(steps[n].bh_acc_rate[i])) steps[n].bh_acc_rate[i] = 0;
//   }

//   if (!(new_bh_mass>=0) && steps[n].log_bh_mass[i] > BH_MASS_TO_REQUIRE_GROWTH) 
//   {
//     //fprintf(stderr, "Negative BH accretion rate! (sm: %e, m: %e, scale: %f; old_bh: %e; new_bh: %e!\n",
//     //sprintf(buffer, "Negative BH accretion rate! (sm: %e, m: %e, scale: %f; old_bh: %e; new_bh: %e!\n",
//     //steps[n].sm_avg[i], exp10(M_MIN+(i+0.5)/BPDEX), steps[n].scale, steps[n].old_bh_mass[i], steps[n].bh_mass_avg[i]);
//     //INVALIDATE(fit, buffer); 
//   }

//   if (steps[n].bh_unmerged[i] < mfrac*new_bh_mass) 
//   {
//     if (steps[n].log_bh_mass[i] > BH_MASS_TO_REQUIRE_GROWTH && steps[n].scale > 0.15 && i>BPDEX && steps[n].t[i] > 1e-8) {
//       //fprintf(stderr, "Merging rate exceeds available mergers! (sm: %e, m: %e, scale: %f; old_bh: %e; new_bh: %e! DBH: %e; Merge: %e; M_avail: %e \n",
//       // sprintf(buffer, "Merging rate exceeds available mergers! (sm: %e, m: %e, scale: %f; old_bh: %e; new_bh: %e! DBH: %e; Merge: %e; M_avail: %e \n",
//         //steps[n].sm_avg[i], exp10(M_MIN+(i+0.5)/BPDEX), steps[n].scale, steps[n].old_bh_mass[i], steps[n].bh_mass_avg[i], new_bh_mass, new_bh_mass*mfrac, steps[n].bh_unmerged[i]);
//       //INVALIDATE(fit, buffer); 
//     } 
//     else 
//     {
//       steps[n].bh_unmerged[i] = 0;
//       // for (j=0; j<MBH_BINS; j++) steps[n].bh_unmerged_dist[i*MBH_BINS+j] = 0;
//     }
//   } 
//   else 
//   {
//     steps[n].bh_unmerged[i] -= mfrac*new_bh_mass;
//     // for (j=0; j<MBH_BINS; j++) steps[n].bh_unmerged_dist[i*MBH_BINS+j] *= steps[n].bh_unmerged[i] / (steps[n].bh_unmerged[i] + mfrac*new_bh_mass);
//   }

//   // // The fraction of the new BH mass that come from mergers relative to the 
//   // // **remaining** unmerged BH mass. 
//   // double frac_from_unmerged = mfrac * new_bh_mass / steps[n].bh_unmerged[i];

//   // steps[n].n_merge10[i] = 0;
//   // steps[n].n_merge100[i] = 0;

//   // // Count the # of mergers with two different mass ratio thresholds
//   // if (n)
//   // {

//   //   // calculate the threshold BH mass and the corresponding bin numbers f10, b10, f100, b100.
//   //   double mbh_thres10 = log10(steps[n].bh_mass_avg[i]) - 1;
//   //   double mbh_thres100 = mbh_thres10 - 2;
//   //   double f10 = (mbh_thres10 - MBH_MIN) * MBH_BPDEX;
//   //   int64_t b10 = f10;
//   //   f10 -= b10;
//   //   double f100 = (mbh_thres100 - MBH_MIN) * MBH_BPDEX;
//   //   int64_t b100 = f100;
//   //   f100 -= b100;

//   //   // make summations over the unmerged BH mass bins above the threshold bin numbers.
//   //   if (b10 >= MBH_BINS) 
//   //     steps[n].n_merge10[i] = 0;
//   //   else
//   //   {
//   //     steps[n].n_merge10[i] += steps[n].bh_unmerged_dist[i*MBH_BINS+b10] * (1 - f10) * frac_from_unmerged;
//   //     for (j=b10+1; j<MBH_BINS; j++) steps[n].n_merge10[i] += steps[n].bh_unmerged_dist[i*MBH_BINS+b10] * frac_from_unmerged;
//   //   }
    
//   //   if (b100 >= MBH_BINS) 
//   //     steps[n].n_merge100[i] = 0;
//   //   else
//   //   {
//   //     steps[n].n_merge100[i] += steps[n].bh_unmerged_dist[i*MBH_BINS+b100] * (1 - f100) * frac_from_unmerged;
//   //     for (j=b100+1; j<MBH_BINS; j++) steps[n].n_merge100[i] += steps[n].bh_unmerged_dist[i*MBH_BINS+b100] * frac_from_unmerged;
//   //   }

//     //steps[n].n_merge10[i] /= steps[n].t[i];
//     //steps[n].n_merge100[i] /= steps[n].t[i];

//   // }

//   // // Count the # of mergers with two different mass ratio thresholds
//   // if (n)
//   // {
//   //   for (j = 0; j < M_BINS; j++)
//   //   {
//   //     int bin = j * M_BINS + i
//   //     // The new unmerged BH mass that come from the average central SMBH of the merged galaxies
//   //     double m_merger_j = steps[n-1].merged[bin] * steps[n-1].bh_mass_avg[j];
//   //     // The fraction that m_merger_j takes up in the total available unmerged BH mass
//   //     double f_merger_j = m_merger_j / steps[n].bh_unmerged_avail[i];
//   //     // If we think that the mergers happen equally likely to the three parts of the available unmerged BHs,
//   //     // then the total merged SMBH mass that come from the average central SMBH of the merged galaxies should
//   //     // be the total merged mass weighted by f_merger_j. The number of mergers would simply be this total
//   //     // merged mass divided by the average BH mass in this bin at the last snapshot.
//   //     double n_merger_j = mfrac * new_bh_mass * f_merger_j / steps[n-1].bh_mass_avg[j];
//   //     if (steps[n-1].bh_mass_avg[j] / steps[n-1].bh_mass_avg[i] > 0.1)
//   //       steps[n].n_merge10[i] += n_merger_j;
//   //     else if (steps[n-1].bh_mass_avg[j] / steps[n-1].bh_mass_avg[i] > 0.01)
//   //       steps[n].n_merge100[i] += n_merger_j;
//   //   }
//   // }

//   // if (steps[n].bh_acc_rate[i] <= 0) steps[n].bh_acc_rate[i] = 1e-8;
//   // steps[n].bh_eta[i] = 
//   //   log10(//schechter_inv_avg(steps[n].smhm.bh_alpha, BHER_MIN, BHER_EFF_MAX)*
//   //   steps[n].bh_acc_rate[i]/steps[n].bh_mass_avg[i]*4.5e8*(steps[n].smhm.bh_efficiency_rad));
//   //   //steps[n].bh_acc_rate[i]/steps[n].bh_mass_avg[i]*4.5e8*(steps[n].smhm.bh_efficiency/(1.0)));
//   // // steps[n].bh_Lbol[i] = steps[n].smhm.bh_efficiency_rad * steps[n].bh_acc_rate[i] * 1.5e46;
//   // // steps[n].bh_sBHAR[i] = steps[n].bh_Lbol[i] / steps[n].sm_avg[i] / 2.6e35;
  

//   if (steps[n].sm_icl[i] < steps[n].icl_frac[i]*steps[n].new_sm[i]) 
//   {
//     if (steps[n].new_sm[i] > 0 && steps[n].scale > 0.15 && i>BPDEX) 
//     {
//       // sprintf(buffer, "ICL depleted (hm: %f; a: %f; sm: %e, new sm: %e, icl needed: %e, icl available: %e)\n",(i+0.5)*INV_BPDEX+M_MIN, steps[n].scale,  steps[n].sm_avg[i], steps[n].new_sm[i], steps[n].new_sm[i] * steps[n].icl_frac[i], steps[n].sm_icl[i]);
//       //fprintf(stderr, "ICL depleted (hm: %f; a: %f; sm: %e, new sm: %e, icl needed: %e, icl available: %e)\n",(i+0.5)*INV_BPDEX+M_MIN, steps[n].scale,  steps[n].sm_avg[i], steps[n].new_sm[i], steps[n].new_sm[i] * steps[n].icl_frac[i], steps[n].sm_icl[i]);
//       INVALIDATE(fit, buffer);
//       steps[n].icl_frac[i] = steps[n].sm_icl[i]/steps[n].new_sm[i];
//     }
//     else 
//       steps[n].icl_frac[i] = 0;
//   }

//   steps[n].merged_frac[i] = 0;

//   if (steps[n].icl_frac[i] && steps[n].sm_icl[i]) 
//   {
//     double frac_from_icl = steps[n].new_sm[i]*steps[n].icl_frac[i] / steps[n].sm_icl[i];
//     if (steps[n].sm_inc[i]>0) 
//     {
//       steps[n].merged_frac[i] = steps[n].new_sm[i]*steps[n].icl_frac[i]/steps[n].sm_inc[i];
//       if (steps[n].merged_frac[i]>1) steps[n].merged_frac[i] = 1;
//     }
//     if (frac_from_icl) 
//     {
//       for (j=0; j<n; j++) 
//       {
//       	steps[n].sm_hist[i*num_outputs + j] += frac_from_icl*steps[n].icl_stars[i*num_outputs + j];
//       	steps[n].icl_stars[i*num_outputs + j] *= (1.0-frac_from_icl);
//       }
//     }
//     steps[n].new_sm[i] *= (1.0-steps[n].icl_frac[i]);
//   }


//   steps[n].new_sm[i] /= steps[n].smloss[n];
//   steps[n].sm_hist[i*num_outputs + n] = steps[n].new_sm[i];
//   //steps[n].sfr[i] = steps[n].new_sm[i] / dt;

//   if ((steps[n].t[i] > REAL_ND_CUTOFF) && (steps[n].sm_avg[i] > 0.17*exp10(steps[n].med_hm_at_a[i])) && (!no_z_scaling) && n>2) 
//   {
//     /*sprintf(buffer, "SM exceeds baryon fraction (sm: %e, m: %e, scale: %f!\n",
// 	    steps[n].sm_avg[i], exp10(M_MIN+(i+0.5)*INV_BPDEX), steps[n].scale);
// 	    INVALIDATE(fit, buffer);*/
//     // Why did we comment that??? This should be an important constraint.
//     //fprintf(stderr, "SM exceeds baryon fraction (sm: %e, m: %e, scale: %f!\n",
//       //steps[n].sm_avg[i], exp10(M_MIN+(i+0.5)*INV_BPDEX), steps[n].scale);
//     // sprintf(buffer, "SM exceeds baryon fraction (sm: %e, m: %e, scale: %f!\n",
//     //   steps[n].sm_avg[i], exp10(M_MIN+(i+0.5)*INV_BPDEX), steps[n].scale);
//       INVALIDATE(fit, buffer);
//   }

//   //if (steps[n].t[i]<REAL_ND_CUTOFF)
//   //  steps[n].sfr[i] = 0;
//   if ((!no_z_scaling) && (steps[n].sfr[i] < 0 || !isfinite(steps[n].sfr[i]))) 
//   {
//     char buffer[1024];
//     // sprintf(buffer, "Negative SFR at a = %f and m = %f (ND=%g)! (ICL m=%f; ICL frac=%f)\n", steps[n].scale, 
// 	   // i*INV_BPDEX + M_MIN, steps[n].t[i], steps[n].smhm.icl_m, steps[n].icl_frac[i]);
//     //fprintf(stderr, "Negative SFR at a = %f and m = %f (ND=%g)! (ICL m=%f; ICL frac=%f)\n", steps[n].scale, 
//       //i*INV_BPDEX + M_MIN, steps[n].t[i], steps[n].smhm.icl_m, steps[n].icl_frac[i]);
//     INVALIDATE(fit, buffer);
//   }
// }

