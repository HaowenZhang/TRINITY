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
#include <omp.h>
#include <assert.h>

#define LogG -5.36456

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

void calc_sfh(struct smf_fit *f) {
  //fprintf(stderr, "z_start: %f", 1 / steps[0].scale - 1);
   int64_t i,j;
  for (i=0; i<num_outputs; i++) {
    if (!no_z_scaling || ((steps[i].scale < 1.0/(z_min+1.0)) &&
			  (steps[i].scale > 1.0/(z_max+1.0))))
    calc_sm_hist(i, f);
    
  }

#pragma omp for schedule(guided,5)
  for (j=1; j<num_outputs; j++) {
    if (!no_z_scaling || ((steps[j].scale < 1.0/(z_min+1.0)) &&
			  (steps[j].scale > 1.0/(z_max+1.0)))) {
      calc_smf_and_ssfr(j, f);
      calc_uvlf(j, f);
      calc_total_sfr(j);

      // calc_bh_acc_rate_distribution(j, f);
      // calc_bh_acc_rate_distribution_full(j, f);
      // calc_active_bh_fraction(j, f);
      // calc_total_bhar(j);
      // calc_observed_bhar(j);
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
  //for (j=0; j<M_BINS; j++) fprintf("sm_hist[0, %d]=%e\n", j, steps[0].sm_hist[j*num_outputs]);
}

void calc_total_sfr(int n) {
  int64_t i;
  double sfr = 0, obs_sfr = 0;
  double mu = steps[n].smhm.mu;
  double sm_corr = pow(10, mu);
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
   if (steps[n].bh_mass_avg[i] < 1e5) continue;
   // if (steps[n].bh_acc_rate[i] < bhar_minimum) continue;
   // double eta = -(l_min + 5.26) / 2.5 - steps[n].log_bh_mass[i];
  double eta = -(l_min + 5.26) / 2.5 - log10(steps[n].bh_mass_avg[i]);
   double eta_frac = eta - steps[n].bh_eta[i];
   bhar += steps[n].bh_acc_rate[i] * steps[n].t[i];
   //fprintf(stderr, "eta_frac=%e\n", eta_frac);
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
 //printf("cosmic_bhar=%e\n", steps[n].cosmic_bhar);
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
  //  double expscale = 1.0; //doexp10(-EXPSCALE(f)*(a*a)*M_LOG10E);
  double incompleteness;
  double a2 = a1*a1;
  // double a4 = a1*a1*a1*a1;
  double z8 = z;
  //if (z8 > 10) z8=10;
  //fprintf(stderr, "z_ceil: %f\n", Z_CEIL(f));
   if (z8 > 15) z8 = 15;
  a = 1.0/(1.0+z);
  a1 = a - 1.0;
  // c.icl_m = ICL_FRAC(f)+a1*ICL_FRAC_E(f);
  c.icl_frac = ICL_FRAC(f)+a1*ICL_FRAC_E(f);
  c.v_1 = M_1(f) + (a1*M_1_A(f) + log(1 + z) * M_1_A2(f)) + z8 * M_1_A3(f);
  // double z8 = z;
  // if (z8 > 12) z8=12;
  c.sm_0 = 0;
  c.epsilon = doexp10(EFF_0(f) + EFF_0_A(f)*a1 + EFF_0_A2(f) * log(1 + z) + EFF_0_A3(f)*z8);
  c.alpha = ALPHA(f) + (a1*ALPHA_A(f) + log(1 + z) * ALPHA_A2(f)) + z8 * ALPHA_A3(f);
  c.delta = DELTA(f); // + (a1*DELTA_A(f) + z8*DELTA_A2(f))*expscale;
  c.beta = BETA(f) + a1*BETA_A(f) + z8*BETA_A2(f);
  c.gamma = GAMMA(f) + (a1*GAMMA_A(f) + z8*GAMMA_A2(f));
  c.gamma = doexp10(c.gamma);
  c.lambda = doexp10(LAMBDA(f) + (a1*LAMBDA_A(f) + z8*LAMBDA_A2(f)));
  c.mu = MU(f) + a1*MU_A(f);
  c.kappa = KAPPA(f) + a1*KAPPA_A(f);
  c.passive_mass = 10.3 + z*0.5 - c.mu;
  c.scatter = SCATTER(f) + a1*SCATTER_A(f);
  //c.scatter = 0.212;
  c.scatter_corr = exp(pow(c.scatter*log(10), 2)/2.0);
  c.obs_scatter = (no_obs_scatter) ? 0 : SIGMA_CENTER + SIGMA_Z(f)*(z-0.1) + SIGMA_A(f)*a1;
  if (c.obs_scatter > 0.3) c.obs_scatter = 0.3;
  c.icl_frac = exp10(ICL_FRAC(f) + a1 * ICL_FRAC_E(f) + z8 * ICL_FRAC_Z(f));
  //c.icl_frac = exp10(ICL_FRAC(f));
  //c.icl_frac_e = ICL_FRAC_E(f);
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

  //c.qv = 2.232131 + 0.217770*a1 + 0.338483*z8;
  //c.qsig = 0.365818 + 0.284595*a1 +  0.081106*z8;

  c.qm = QM_0(f) + QM_1(f) *a1 + QM_2(f) *z8;
  c.qwidth = QWIDTH_0(f) + QWIDTH_1(f) *a1 + QWIDTH_2(f) * z8;

  //c.qm = 12.214 + 0.0592*a1 + 0.139*z8;
  //c.qwidth = 0.498 + 0.372*a1 + 0.0124 * z8;

  // // BH part
  // c.bh_beta = BH_BETA_0(f) +BH_BETA_1(f)*a1 + BH_BETA_2(f)*z8;
  // c.bh_gamma = BH_GAMMA_0(f) + BH_GAMMA_1(f)*a1 + BH_GAMMA_2(f)*z8;
  // c.bh_merge = BH_MERGE_F_0(f) + BH_MERGE_F_1(f)*a1;
  // c.bh_merge_width = BH_MERGE_W_0(f) + BH_MERGE_W_1(f)*a1;
  // c.bh_alpha = BH_ALPHA_0(f) + BH_ALPHA_1(f)*a1;
  // c.bh_delta = BH_DELTA_0(f) + BH_DELTA_1(f)*a1;

  // if (c.bh_alpha > c.bh_delta)
  // {
  //   double tmp = c.bh_alpha;
  //   c.bh_alpha = c.bh_delta;
  //   c.bh_delta = tmp;
  // }

  // c.bh_efficiency = pow(10, BH_EFFICIENCY_0(f)+BH_EFFICIENCY_1(f)*a1);
  // if (c.bh_efficiency > 0.7) c.bh_efficiency = 0.7;
  // c.bh_efficiency_rad = c.bh_efficiency / (1.0-c.bh_efficiency);
  // c.bh_duty = BH_DUTY_0(f)+BH_DUTY_1(f)*a1;
  // if (c.bh_duty < 1e-4) c.bh_duty = 1e-4;
  // c.bh_scatter = BH_SCATTER_0(f) + BH_SCATTER_1(f)*a1;
  // c.bh_prob_norm = 1.0; //Since in the current model, the scatter is larger than zero and we are renormalizing the bher_dist_full each time, 
  //                                               //we don't really have to use this quantity.
  // c.dc_mbh = DC_MBH_0(f);
  // c.dc_mbh_w = DC_MBH_W_0(f);
  // c.eta_mu = ETA_MU_0(f);
  // c.rho_bh = RHO_BH(f);

  double comb_scatter = sqrt(c.obs_scatter*c.obs_scatter + c.scatter*c.scatter);
  double sm_at_m14 = 11.5; //calc_sm_at_m(14, c)+c.mu;
  double passive_frac = 1.0/(exp10fc(-1.3*(sm_at_m14-c.passive_mass))+1);
  double avg_offset_passive = (1.0-passive_frac)*c.kappa;
  double avg_offset_active = c.kappa - avg_offset_passive;
  c.combined_scatter = sqrt(pow(comb_scatter,2) + passive_frac*pow(avg_offset_passive,2) + (1.0-passive_frac)*pow(avg_offset_active,2)); 
  return c;
}

double calc_bh_at_bm(double bm, struct smf c) {
  return c.bh_beta + c.bh_gamma * (bm - 11.0);
}

double calc_bh_at_vdisp(double vd, struct smf c) {
  return c.bh_beta + c.bh_gamma * (vd - 2.30103);
}

/*double calc_sm_at_m(double m, struct smf c) {
  double dm = m-c.m_1;
  double sm, dma = dm*c.alpha;
  //if (dma > dmb+20) sm = -dma;
  //else if (dmb > dma+20) sm = -dmb;
  //else sm = -dolog10(doexp10(dma) + doexp10(dmb));
  sm = -log10_1p_exp10fc(dma);
  //sm += c.sm_0 + c.delta*(exp10(log10fc(log10_1p_exp10fc(dm))*c.gamma)/(1.0+exp(exp10(-dm)))) + c.f_1;
  sm += c.sm_0 + c.delta*(exp10(dolog10(log10(1.0+exp(dm)))*c.gamma)/(1.0+exp(exp10(-dm)))) + c.f_1;
  return sm;
  }*/

double calc_sfr_at_lv(double m, double lv, struct smf c) {
  double vd = lv - c.v_1;
  double vd2 = vd/c.delta;
  // double vd3 = (lv-c.qv)/c.qsig;
  // double qfrac = c.fqmin + (1.0 - c.fqmin) * (0.5 + 0.5 * erf(vd3*M_SQRT1_2));
  
  double sfrac = 1 / (1 + exp((m - c.qm) / c.qwidth));
  //double sfrac = 0.5+0.5*erf(-vd3*M_SQRT1_2);
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


double calc_smf_at_m(double m, double sm, int64_t n, struct smf_fit *fit, gsl_interp_accel *ga) {
  if (sm < steps[n].smhm.sm_min) 
    {
      // printf("at z=%f, sm=%f < sm_min\n", (1 / steps[n].scale - 1, sm));
      return 1e-17;
    }
  if (sm > steps[n].smhm.sm_max) 
  {
    // printf("at z=%f, sm=%f > sm_max\n", (1 / steps[n].scale - 1, sm));
    return 1e-17;
  }
  if (!steps[n].smhm.valid) 
  {
    // printf("at z=%f, sm=%f, model not valid.\n", (1 / steps[n].scale - 1, sm));
    return 1e-17;
  }
  double dlm_dlsm;
  int err = gsl_spline_eval_deriv_e(steps[n].spline, sm, ga, &dlm_dlsm);

  if (err)
  {
    // sprintf(buffer, "Error in GSL spline interpolation #4.\n");
    //fprintf(stderr, "Error in GSL spline interpolation #4.\n");
    // INVALIDATE(fit, buffer);
    return 0;
  }
  dlm_dlsm = fabs(dlm_dlsm);
  double dndlogm = mf_cache(steps[n].scale, m);
  // double z = 1 / steps[n].scale - 1;
  // if (z > 3.5 && z < 4.0)
	 //  printf("at z=%f, m=%f, sm=%f, sm_max=%f, dlm_dlsm=%f, dndlogm=%f.\n", z, m, sm, steps[n].smhm.sm_max, dlm_dlsm, dndlogm);
  // if (!(dlm_dlsm>0)) { // && fit) {
  //   // printf("at z=%f, sm=%f, Falling SMHM.\n", (1 / steps[n].scale - 1), sm);
  //   // INVALIDATE(fit, "Falling SM(M) relation.");
  //   // return 0;
  //   dlm_dlsm = 1.0;
  // }
  // double dndlogm = mf_cache(steps[n].scale, m);
  double phi = doexp10(dndlogm)*dlm_dlsm;
  if (!isfinite(phi)) phi = 0;
  return phi;
}

double calc_smf_at_sm(int64_t n, double sm) {
  if (sm < steps[n].smhm.sm_min) return 1e-17;
  if (sm > steps[n].smhm.sm_max) return 1e-17;
  if (!steps[n].smhm.valid) return 1e-17;
  gsl_interp_accel ga = {0};
  // printf("sm_min=%f, sm_max=%f\n", steps[n].smhm.sm_min, steps[n].smhm.sm_max);
  double m;
  int err = gsl_spline_eval_e(steps[n].spline, sm, &ga, &m);

  if (err)
    {
      // sprintf(buffer, "Error in GSL spline interpolation #1.\n");
      // fprintf(stderr, "Error in GSL spline interpolation #1.\n");
      // INVALIDATE(fit, buffer);
      return 0;
    }
  double result = calc_smf_at_m(m, sm, n, NULL, &ga);
  gsl_interp_accel_reset(&ga);
  if (steps[n].alloc2smf)
  {
    err = gsl_spline_eval_e(steps[n].spline2, sm, &ga, &m); //Second segment.
    if (!err) result += calc_smf_at_m(m, sm, n, NULL, &ga);
  }
  
  return result;
}




double calc_uvlf_at_m(double m, double uv, int64_t n, struct smf_fit *fit, gsl_interp_accel *ga) {
  if (uv < steps[n].smhm.uv_min) 
    {
      // printf("at z=%f, sm=%f < sm_min\n", (1 / steps[n].scale - 1, sm));
      return 1e-17;
    }
  if (uv > steps[n].smhm.uv_max) 
  {
    // printf("at z=%f, sm=%f > sm_max\n", (1 / steps[n].scale - 1, sm));
    return 1e-17;
  }
  if (!steps[n].smhm.valid) 
  {
    // printf("at z=%f, sm=%f, model not valid.\n", (1 / steps[n].scale - 1, sm));
    return 1e-17;
  }
  double dlm_dMuv;
  int err = gsl_spline_eval_deriv_e(steps[n].spline_uv, uv, ga, &dlm_dMuv);

  dlm_dMuv = fabs(dlm_dMuv);

  if (err)
  {
    // sprintf(buffer, "Error in GSL spline interpolation #4.\n");
    //fprintf(stderr, "Error in GSL spline interpolation #4.\n");
    // INVALIDATE(fit, buffer);
    return 0;
  }
  double dndlogm = mf_cache(steps[n].scale, m);
  // double z = 1 / steps[n].scale - 1;
  // if (z > 3.5 && z < 4.0)
   //  printf("at z=%f, m=%f, sm=%f, sm_max=%f, dlm_dlsm=%f, dndlogm=%f.\n", z, m, sm, steps[n].smhm.sm_max, dlm_dlsm, dndlogm);
  if (!(dlm_dMuv>0)) { // && fit) {
    // printf("at z=%f, sm=%f, Falling SMHM.\n", (1 / steps[n].scale - 1), sm);
    // INVALIDATE(fit, "Falling SM(M) relation.");
    // return 0;
    dlm_dMuv = 1.0;
  }
  // double dndlogm = mf_cache(steps[n].scale, m);
  double phi = doexp10(dndlogm)*dlm_dMuv;
  if (!isfinite(phi)) phi = 0;
  return phi;
}

double calc_uvlf_at_uv(int64_t n, double uv) {
  if (uv < steps[n].smhm.uv_min) return 1e-17;
  if (uv > steps[n].smhm.uv_max) return 1e-17;
  if (!steps[n].smhm.valid) return 1e-17;
  gsl_interp_accel ga = {0};
  // printf("sm_min=%f, sm_max=%f\n", steps[n].smhm.sm_min, steps[n].smhm.sm_max);
  double m;
  int err = gsl_spline_eval_e(steps[n].spline_uv, uv, &ga, &m);
  if (err)
    {
      // sprintf(buffer, "Error in GSL spline interpolation #1.\n");
      // fprintf(stderr, "Error in GSL spline interpolation #1.\n");
      // INVALIDATE(fit, buffer);
      return 0;
    }
  double result = calc_uvlf_at_m(m, uv, n, NULL, &ga);
  gsl_interp_accel_reset(&ga);
  if (steps[n].alloc2uvlf)
  {
    err = gsl_spline_eval_e(steps[n].spline_uv2, uv, &ga, &m); //Remember to add the second segment.
    if (!err) result += calc_uvlf_at_m(m, uv, n, NULL, &ga);
  }
  return result;
}



// void calc_active_bh_fraction(int n, struct smf_fit *fit) {
//   int64_t i;
//   for (i=0; i<M_BINS; i++) {
//     double dc = steps[n].smhm.bh_duty;
//     double alpha = steps[n].smhm.bh_alpha;
//     double delta = steps[n].smhm.bh_delta;
//     // modulate the duty cycle with a power-law dependence on mass.
//     // dc *= pow(log10(steps[n].bh_mass_avg[i]) / log10(steps[n].bh_mass_avg[M_BINS - 1]), steps[n].smhm.dc_beta);
//     // dc *= pow((M_MIN + (i + 0.5) * INV_BPDEX) / (M_MAX - 0.5 * INV_BPDEX), steps[n].smhm.dc_beta);
//     // double f_mass = pow(exp10((M_MIN + (i + 0.5) * INV_BPDEX) - (M_MAX - 0.5 * INV_BPDEX)), steps[n].smhm.dc_beta);
//     double f_mass = exp((log10(steps[n].bh_mass_avg[i]) - steps[n].smhm.dc_mbh) / steps[n].smhm.dc_mbh_w);
//     f_mass = f_mass / (1 + f_mass);
//     // f_mass = f_mass < 1? f_mass : 1;
//     dc *= f_mass;
//     if (dc < 1e-4) dc = 1e-4;
//     double sn = dc/doublePL_frac_above_thresh(BHER_EFF_MIN, alpha, delta, 1.0, fit);
//     double bh_eta = steps[n].bh_eta[i];
//     double thresh = -2.0 - bh_eta;
//     steps[n].f_active[i] = doublePL_frac_above_thresh(thresh, alpha, delta, sn, fit);
//   }
// }

void calc_active_bh_fraction(int n, struct smf_fit *fit) {
  int64_t i, j;
  for (i=0; i<M_BINS; i++) {

    double ledd_min = steps[n].ledd_min[i];
    double ledd_max = steps[n].ledd_max[i];

    double bpdex = steps[n].ledd_bpdex[i];
    double inv_bpdex = 1.0/bpdex;

    double eta_frac = -2.0 - steps[n].bh_eta[i];

    

    //if (n == 177)
    //{
      //fprintf(stderr, "n=%d, i=%ld, eta_frac=%f, ledd_min=%f, ledd_max=%f\n", n, i, eta_frac, steps[n].ledd_min[i], steps[n].ledd_max[i]);
    //}
    if (eta_frac < ledd_min || eta_frac > ledd_max) continue; 
    
    double bher_f = (eta_frac-ledd_min)*bpdex;
    int64_t bher_b = bher_f;
    bher_f -= bher_b;

    double p_good = (1 - bher_f) * steps[n].bher_dist[i*BHER_BINS + bher_b+1]; //The area of the histogram with eta > -2
    double p_total = steps[n].bher_dist_norm[i]; //Total area

    if (!isfinite(p_good) || !isfinite(p_total) || p_total == 0)
    {
      steps[n].f_active[i] = 0;
      continue;
    }

    for (j = bher_b + 2; j < BHER_BINS; j++)
    {
      p_good += steps[n].bher_dist[i*BHER_BINS + j];
      // p_total += steps[n].bher_dist[i*BHER_BINS + j];
    }

    // for (j = 0; j < bher_b + 2; j++)
    // {
    //   p_total += steps[n].bher_dist[i*BHER_BINS + j];
    // }

    double dc = steps[n].smhm.bh_duty;
    double f_mass = exp((log10(steps[n].bh_mass_avg[i]) - steps[n].smhm.dc_mbh) / steps[n].smhm.dc_mbh_w);
    f_mass = f_mass / (1 + f_mass);
    // f_mass = f_mass < 1? f_mass : 1;
    dc *= f_mass;
    if (dc < 1e-4) dc = 1e-4;

    steps[n].f_active[i] = p_good / p_total * dc;
    //if (n == 177)
    //{
      //fprintf(stderr, "n=%d, i=%ld, p_good=%e, p_total=%e, dc=%e\n", n, i, p_good, p_total, dc);
    //}
  }
}

void calc_bh_acc_rate_distribution(int n, struct smf_fit *fit) {
  int64_t i, j;
  memset(steps[n].bher_dist, 0, sizeof(double)*BHER_BINS);
  double sm_scatter = steps[n].smhm.scatter*steps[n].smhm.bh_gamma;
  double scatter_bhar = 0.2; //assuming that the scatter in BHAR at fixed halo mass is 0.2 dex (kind of arbitrary...)
  double scatter = sqrt(sm_scatter*sm_scatter + steps[n].smhm.bh_scatter*steps[n].smhm.bh_scatter);
  // The total scatter that we should use to smooth the ERDF is the joint scatter of these two correlated gaussians,
  // with rho_bh being the correlation coefficient.
  //scatter = sqrt(scatter * scatter + scatter_bhar * scatter_bhar - 2 * steps[n].smhm.rho_bh * scatter * scatter_bhar);

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
    steps[n].bh_eta[i] += bh_eta_corr;
  }

  if (scatter <= 0) {
    for (i=0; i<BHER_BINS; i++) {
      double bher = BHER_MIN+i*BHER_INV_BPDEX;
      // steps[n].bher_dist[i] = exp10(bher*steps[n].smhm.bh_alpha)*exp(-exp10(bher));
      steps[n].bher_dist[i] = 1 / (exp10(bher*steps[n].smhm.bh_alpha) * exp10(bher*steps[n].smhm.bh_delta));
    }
  } else {
    const int64_t cmin = -BHER_BINS-3*BHER_BPDEX;
    const int64_t cmax = 3*BHER_BPDEX + BHER_BINS;
    double *ecache = check_realloc(NULL, sizeof(double)*(cmax-cmin+1), "Allocating ecache");
    double inv_scatter = 1.0/scatter;
    for (j=cmin; j<=cmax; j++) {
      double db = inv_scatter*j*BHER_INV_BPDEX;
      ecache[j-cmin] = exp(-0.5*db*db);
    }
    for (j=-3*BHER_BPDEX; j<BHER_BINS+3*BHER_BPDEX; j++) {
      double bher = BHER_MIN+j*BHER_INV_BPDEX;
      double prob = 1 / (exp10(bher*steps[n].smhm.bh_alpha) * exp10(bher*steps[n].smhm.bh_delta));
      double *ec = ecache + (-j - cmin);
      for (i=0; i<BHER_BINS; i++) {
  //double bher2 = BHER_MIN+i*BHER_INV_BPDEX;
  //double db = inv_scatter*(bher - bher2);
  double weight = ec[i]; //exp(-0.5*db*db);
  //printf("%e %e %"PRId64" %"PRId64" %f %"PRId64"\n", weight, ec[i], i, j, db, i-j-cmin);  
  steps[n].bher_dist[i] += weight*prob;
      }
    }
    free(ecache);
  }

  //Normalize;
  double total = 0;
  for (i=0; i<BHER_BINS; i++) total += steps[n].bher_dist[i];
  total *= BHER_INV_BPDEX;
  //assert(total > 0);
  if (total > 0)
  {
  total = 1.0 / total;
  }
  for (i=0; i<BHER_BINS; i++) steps[n].bher_dist[i] *= total;

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
  memset(steps[n].bher_dist_full, 0, sizeof(float)*BHER_BINS*M_BINS);
  memset(steps[n].bher_dist, 0, sizeof(double)*BHER_BINS*M_BINS);
  memset(steps[n].bher_dist_norm, 0, sizeof(double)*M_BINS);


  double sm_scatter = steps[n].smhm.scatter*steps[n].smhm.bh_gamma;
  double scatter_bhar = 0.2; //assuming that the scatter in BHAR at fixed halo mass is 0.2 dex (kind of arbitrary...)
  double scatter = sqrt(sm_scatter*sm_scatter + steps[n].smhm.bh_scatter*steps[n].smhm.bh_scatter);
  // The total scatter that we should use to smooth the ERDF is the joint scatter of these two correlated gaussians,
  // with rho_bh being the correlation coefficient.
  //scatter = sqrt(scatter * scatter + scatter_bhar * scatter_bhar - 2 * steps[n].smhm.rho_bh * scatter * scatter_bhar);

  //  Already corrected in bh_acc_rate_distribution...
  //  double bh_eta_corr = log10(schechter_inv_avg(steps[n].smhm.bh_alpha, BHER_MIN, BHER_EFF_MAX, fit)/steps[n].smhm.bh_duty);
  //  for (i=0; i<M_BINS; i++) steps[n].bh_eta[i] += bh_eta_corr;

  double (*bher_prob)(double, double, double, double, double) = &_prob_of_ledd_linear;
  if (nonlinear_luminosity) bher_prob = &_prob_of_ledd_nonlinear;
  //assert(bher_prob == &_prob_of_ledd_linear);

  if (scatter <= 0) {
    for (i=0; i<M_BINS; i++) {
      for (j=0; j<BHER_BINS; j++) {
  double bher = BHER_MIN+j*BHER_INV_BPDEX;
  steps[n].bher_dist_full[i*BHER_BINS+j] = bher_prob(bher+steps[n].bh_eta[i], steps[n].bh_eta[i], steps[n].smhm.bh_alpha, steps[n].smhm.bh_delta, 1.0)*steps[n].smhm.bh_prob_norm;
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
      float *bher_dist = steps[n].bher_dist_full+i*BHER_BINS;
      double ledd_min = steps[n].bh_eta[i] + BHER_MIN;
      if (nonlinear_luminosity && ledd_min < -2.0)
  ledd_min = (ledd_min + 1.0)*2.0;
      ledd_min -= steps[n].bh_eta[i];
      double ledd_eff_min = steps[n].bh_eta[i] + BHER_EFF_MIN;
      // if (nonlinear_luminosity && ledd_min < -2.0) // This is what Peter has written, which is a bug...
      if (nonlinear_luminosity && ledd_eff_min < -2.0)
  ledd_eff_min = (ledd_eff_min + 1.0)*2.0;
      ledd_eff_min -= steps[n].bh_eta[i];
      double ledd_max = steps[n].bh_eta[i] + BHER_MAX;
      if (nonlinear_luminosity && ledd_max < -2.0)
  ledd_max = (ledd_max + 1.0)*2.0;
      ledd_max -= steps[n].bh_eta[i];


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

      steps[n].ledd_min[i] = ledd_min;
      steps[n].ledd_max[i] = ledd_max;
      steps[n].ledd_bpdex[i] = bpdex;

      for (k=kmin; k<kmax; k++) {
  double db = inv_scatter*k*inv_bpdex;
  ecache[k-kmin] = exp(-0.5*db*db);
      }

      int64_t jmin = (ledd_eff_min-ledd_min)*bpdex;
      int64_t jmax = (ledd_max-ledd_min)*bpdex;

      for (j=jmin; j<jmax; j++) {
  float bher = ledd_min+j*inv_bpdex;
  int64_t emin = j+kmin;
  if (emin < 0) emin = 0;
  int64_t emax = j+kmax;
  if (emax>BHER_BINS-1) emax = BHER_BINS-1;

  // if ((j + kmin < 0))
  // {
  //   for (k=0; k<BHER_BINS; k++) bher_dist[k] = 0;
  //   INVALIDATE(fit, "Not good j + kmin. INVALIDATE.\n");
  //   free(ecache);
  //   return;
  // }

  // float *ec = ecache - (j+kmin);
  float *ec = ecache - emin;



  float prob = bher_prob(bher+steps[n].bh_eta[i], steps[n].bh_eta[i], steps[n].smhm.bh_alpha, steps[n].smhm.bh_delta, 1.0);
  if (prob < 1e-15) {
    if (steps[n].smhm.bh_alpha > 0 && steps[n].smhm.bh_delta > 0) break;
    else if (bher+steps[n].bh_eta[i] > 0) break;
    else continue;
  }

  //bher_dist[j] = prob;
  steps[n].bher_dist[i*BHER_BINS+j] = prob; //Use this to calculate the active fraction later.
  steps[n].bher_dist_norm[i] += prob;
  
  for (k=emin; k<emax; k++) {
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

  /*  char buffer[1024];
  sprintf(buffer, "bher_dist/dist_%f.dat\n", steps[n].scale);
  FILE *out = check_fopen(buffer, "w");
  fprintf(out, "#ER Prob.\n");
  for (i=0; i<BHER_BINS; i++) fprintf(out, "%f %e\n", BHER_MIN+i*BHER_INV_BPDEX, steps[n].bher_dist[i]);
  fclose(out);
  */
}

void calc_smf_and_ssfr(int n, struct smf_fit *fit) {
  int64_t i, j;
  char buffer[1024];
  double m, m2, sm, sm_max=0, sm_min=1000;
  int64_t bin_min=-1, bin_max=-1, bin_peak = 0;
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
    
     
    // INPLEMENTAION 2: bin_max - bin_min + 1 < M_BINS

    if (steps[n].log_sm[i]>steps[n].smhm.sm_max) {
      sm_max = steps[n].smhm.sm_max = steps[n].log_sm[i];
      bin_max = i;
    }
    // if (steps[n].log_sm[i]>-1000 && bin_min < 0) {
    if (steps[n].log_sm[i]>0 && bin_min < 0) {
      bin_min = i;
      sm_min = steps[n].smhm.sm_min = steps[n].log_sm[i];
    }

    if (i && (steps[n].log_sm[i] <= steps[n].log_sm[i-1])) {
      // if (steps[n].log_sm[i]==0) { steps[n].log_sm[i]=0; } // I totally don't understand why you would do this.
      //                                                           // In this case at the very massive end of the SMF,
      //                                                           // the corresponding halo mass would be very difficult
      //                                                           // to determine, and cause the crazy behavior there.
      //                                                           // For now I just put an empty block here.
      if (steps[n].log_sm[i]==0) 
        { 
          // if (bin_min >= 0) steps[n].log_sm[i] = steps[n].log_sm[i-1]+0.001; 
          // steps[n].log_sm[i] = steps[n].log_sm[i-1]+0.001; 
          // if (i > bin_max) bin_max = i; //Even if we don't need this artificial massive end, we still
          //                                                       // need to maintain the bin_max to keep the size consistency between
          //                                                       // the interpolated data and the gsl_spline object.
          if (M_MIN+i*INV_BPDEX < 12) 
          {
            ;
            // steps[n].log_sm[i] = steps[n].log_sm[i-1]+0.001;
          } 
          else //at the high mass end, we do a constant SM/HM ratio extrapolation.
          {
            steps[n].log_sm[i] = steps[n].log_sm[i-1]+INV_BPDEX;
          }
          //bin_max = i;
          if (steps[n].smhm.sm_max < steps[n].log_sm[i]) 
	  {
		  sm_max = steps[n].smhm.sm_max = steps[n].log_sm[i];
		  bin_max = i;
	  }
        }
      else
      {
        if (!count_falling)
        {
          count_falling = 1;
          bin_peak = i - 1;
        }
      }
      
    }

  }


    if (bin_min < 0 || bin_max - bin_min + 1 < 10 || (count_falling && bin_peak <= bin_min + 1))
    {
            //sprintf(buffer, "All the stellar masses are zero.\n");
            //fprintf(stderr, "All the stellar masses are zero at z=%f.\n", 1.0 / steps[n].scale - 1);
            if (n > 0) INVALIDATE(fit, buffer);
            return;
    }


  
  double *log_sm_tmp = NULL;
  double  *hm_tmp = NULL;
  double  *sfr_tmp = NULL;

  if (!count_falling)
  {
    log_sm_tmp = malloc(sizeof(double)*(bin_max - bin_min + 1));
    hm_tmp = malloc(sizeof(double)*(bin_max - bin_min + 1));
    sfr_tmp = malloc(sizeof(double)*(bin_max - bin_min + 1));
    for (j=bin_min; j <= bin_max; j++)
    {
      // printf("log_sm=%f\n", steps[n].log_sm[j]);
      log_sm_tmp[j - bin_min] = steps[n].log_sm[j];
      hm_tmp[j - bin_min] = steps[n].med_hm_at_a[j];
      // if (n == 18) printf("hm_tmp[%d]=%f, med_hm_at_a[%d]=%f\n", j - bin_min, hm_tmp[j - bin_min], j, steps[n].med_hm_at_a[j]);
      sfr_tmp[j - bin_min] = steps[n].sfr[j];
    }
    // The spline does not need to be modified. What should be changed is that the dn/dlogm
    // should be weighted by sfrac and 1 - sfrac for smf_sf and smf_q, respectively.
    steps[n].spline = gsl_spline_alloc(gsl_interp_cspline, bin_max - bin_min + 1);
    steps[n].spline_sfr = gsl_spline_alloc(gsl_interp_cspline, bin_max - bin_min + 1);
    steps[n].flag_alloc = 1;
    gsl_spline_init(steps[n].spline, log_sm_tmp, hm_tmp, bin_max - bin_min + 1);
    gsl_spline_init(steps[n].spline_sfr, log_sm_tmp, sfr_tmp, bin_max - bin_min + 1);
    free(log_sm_tmp); free(hm_tmp); free(sfr_tmp);
  }

  else
  {
    log_sm_tmp = malloc(sizeof(double)*(bin_peak - bin_min + 1));
    hm_tmp = malloc(sizeof(double)*(bin_peak - bin_min + 1));
    sfr_tmp = malloc(sizeof(double)*(bin_peak - bin_min + 1));
    for (j=bin_min; j <= bin_peak; j++)
    {
      // printf("log_sm=%f\n", steps[n].log_sm[j]);
      log_sm_tmp[j - bin_min] = steps[n].log_sm[j];
      hm_tmp[j - bin_min] = steps[n].med_hm_at_a[j];
      // if (n == 18) printf("hm_tmp[%d]=%f, med_hm_at_a[%d]=%f\n", j - bin_min, hm_tmp[j - bin_min], j, steps[n].med_hm_at_a[j]);
      sfr_tmp[j - bin_min] = steps[n].sfr[j];

    }

    // The spline does not need to be modified. What should be changed is that the dn/dlogm
    // should be weighted by sfrac and 1 - sfrac for smf_sf and smf_q, respectively.
    steps[n].spline = gsl_spline_alloc(gsl_interp_cspline, bin_peak - bin_min + 1);
    steps[n].spline_sfr = gsl_spline_alloc(gsl_interp_cspline, bin_peak - bin_min + 1);
    gsl_spline_init(steps[n].spline, log_sm_tmp, hm_tmp, bin_peak - bin_min + 1);
    gsl_spline_init(steps[n].spline_sfr, log_sm_tmp, sfr_tmp, bin_peak - bin_min + 1);
    steps[n].flag_alloc = 1;
    free(log_sm_tmp); free(hm_tmp); free(sfr_tmp);


    int n_seg = M_BINS - bin_peak > 2 ? M_BINS - bin_peak : 3;

    log_sm_tmp = malloc(sizeof(double)*n_seg);
    hm_tmp = malloc(sizeof(double)*n_seg);
    sfr_tmp = malloc(sizeof(double)*n_seg);



    // log_sm_tmp = (double *)realloc(log_sm_tmp, sizeof(double)*(M_BINS - bin_peak + 1));
    // hm_tmp = (double *)realloc(hm_tmp, sizeof(double)*(M_BINS - bin_peak + 1));
    // sfr_tmp = (double *)realloc(sfr_tmp, sizeof(double)*(M_BINS - bin_peak + 1));
    if (M_BINS - bin_peak > 2)
    {
      for (j=bin_peak; j < M_BINS; j++)
      {
        // printf("log_sm=%f\n", steps[n].log_sm[j]);
        log_sm_tmp[j - bin_peak] = steps[n].log_sm[M_BINS - 1 - (j - bin_peak)];
        hm_tmp[j - bin_peak] = steps[n].med_hm_at_a[M_BINS - 1 - (j - bin_peak)];
        // if (n == 18) printf("hm_tmp[%d]=%f, med_hm_at_a[%d]=%f\n", j - bin_min, hm_tmp[j - bin_min], j, steps[n].med_hm_at_a[j]);
        sfr_tmp[j - bin_peak] = steps[n].sfr[M_BINS - 1 - (j - bin_peak)];
      }
    }

    else
    {
      log_sm_tmp[2] = steps[n].log_sm[bin_peak];
      log_sm_tmp[1] = steps[n].log_sm[bin_peak + 1];
      log_sm_tmp[0] = 2 * log_sm_tmp[1] - log_sm_tmp[2];
      hm_tmp[0] = 16.3; hm_tmp[1] = 16.1; hm_tmp[2] = 15.9;
    }
    

    // if (n_seg > M_BINS - bin_peak)
    // {
    //   hm_tmp[2] = 16.3;
    //   log_sm_tmp[2] = 2 * log_sm_tmp[1] - log_sm_tmp[0];
    //   sfr_tmp[2] = 2 * sfr_tmp[1] - sfr_tmp[0];
    // }

    steps[n].spline2 = gsl_spline_alloc(gsl_interp_cspline, n_seg);
    steps[n].spline_sfr2 = gsl_spline_alloc(gsl_interp_cspline, n_seg);
    int err_spline_init1 = gsl_spline_init(steps[n].spline2, log_sm_tmp, hm_tmp, n_seg);
    int err_spline_init2 = gsl_spline_init(steps[n].spline_sfr2, log_sm_tmp, sfr_tmp, n_seg);
    free(log_sm_tmp); free(hm_tmp); free(sfr_tmp);
    if ((err_spline_init1) || (err_spline_init2))
    {
      //fprintf(stderr, "More than 1 turning point in the SMHM.\n");
      if (1.0 / steps[n].scale - 1 < 5) INVALIDATE(fit, buffer);
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
  // steps[n].spline = gsl_spline_alloc(gsl_interp_cspline, bin_max - bin_min + 1);
  // steps[n].spline_sfr = gsl_spline_alloc(gsl_interp_cspline, bin_max - bin_min + 1);
  // steps[n].flag_alloc = 1;
  // gsl_spline_init(steps[n].spline, log_sm_tmp, hm_tmp, bin_max - bin_min + 1);
  // gsl_spline_init(steps[n].spline_sfr, log_sm_tmp, sfr_tmp, bin_max - bin_min + 1);

  i = M_BINS-1;
  for (j=SM_BINS-1; j>=0; j--) 
  {
    sm = SM_MIN + (double)j*SM_INV_BPDEX;
    if (sm > sm_max || sm < sm_min) 
    {
      steps[n].sfr_sm[j+SM_EXTRA] = 0;
      steps[n].smf[j+SM_EXTRA] = 0;
      //steps[n].smf_sf[j+SM_EXTRA] = 0;
      //steps[n].smf_q[j+SM_EXTRA] = 0;
    }
    else 
    {

      int err = gsl_spline_eval_e(steps[n].spline, sm, &ga, &m);
      // fprintf(stdout, "n=%d, z=%f, sm=%f, interpolated m=%f, err=%d\n", n, 1.0 / steps[n].scale - 1, sm, m, err);
      if (err)
      //if (err || m < M_MIN - 1.0 || m > M_MAX + 1.0)
      {
        // sprintf(buffer, "Error in GSL spline interpolation #2.\n");
        //fprintf(stderr, "Error in GSL spline interpolation #2. m=%f\n", m);
        INVALIDATE(fit, buffer);
        // free(log_sm_tmp); free(hm_tmp); free(sfr_tmp);
        return;
      }
      //int64_t b = gsl_interp_accel_find(&ga, steps[n].log_sm+bin_min, (bin_max-bin_min+1), sm)+bin_min;
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

      

      steps[n].sfrac_sm[j+SM_EXTRA] = steps[n].sfrac[b] + f * (steps[n].sfrac[b+1] - steps[n].sfrac[b]);
      steps[n].smf[j+SM_EXTRA] = calc_smf_at_m(m,sm,n,fit,&ga);
      gsl_interp_accel_reset(&ga);
      err = gsl_spline_eval_e(steps[n].spline_sfr, sm, &ga, &(steps[n].sfr_sm[j + SM_EXTRA]));
      steps[n].sfr_sm[j + SM_EXTRA] *= steps[n].smf[j+SM_EXTRA];
	// steps[n].sfr_sm[j+SM_EXTRA] = gsl_spline_eval(steps[n].spline_sfr, sm, &ga)*steps[n].smf[j+SM_EXTRA];
      if (err)
	//if (err || steps[n].sfr_sm[j + SM_EXTRA] < 0 || !(isfinite(steps[n].sfr_sm[j+SM_EXTRA])))
      {
        //sprintf(buffer, "Error in GSL spline interpolation #3.\n");
        fprintf(stderr, "Error in GSL spline interpolation #3. sfr=%e\n", steps[n].sfr_sm[j+SM_EXTRA]);
        INVALIDATE(fit, buffer);
        // free(log_sm_tmp); free(hm_tmp); free(sfr_tmp);
        return;
      }

      if (steps[n].alloc2smf)
      {
        gsl_interp_accel_reset(&ga);
        int err2 = gsl_spline_eval_e(steps[n].spline2, sm, &ga, &m2);
        if (!err2)
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
          // Given the second branch of the SMHM, we need to modify the sfrac_sm value
          // s.t. it is the weighted average of the two HM bins.
          steps[n].sfrac_sm[j+SM_EXTRA] *= steps[n].t[b] + f * (steps[n].t[b+1] - steps[n].t[b]);
          steps[n].sfrac_sm[j+SM_EXTRA] += (steps[n].sfrac[b2] + f2 * (steps[n].sfrac[b2+1] - steps[n].sfrac[b2])) *
                                            (steps[n].t[b2] + f2 * (steps[n].t[b2+1] - steps[n].t[b2]));
          steps[n].sfrac_sm[j+SM_EXTRA] /= (steps[n].t[b] + f * (steps[n].t[b+1] - steps[n].t[b]) + 
                                            steps[n].t[b2] + f2 * (steps[n].t[b2+1] - steps[n].t[b2]));
          double smf2 = calc_smf_at_m(m2,sm,n,fit,&ga);
          steps[n].smf[j+SM_EXTRA] += smf2;
          double sfr_sm_tmp;
          gsl_interp_accel_reset(&ga);
          double err_sfr = gsl_spline_eval_e(steps[n].spline_sfr2, sm, &ga, &(sfr_sm_tmp));
          if (!err_sfr)
          {
            steps[n].sfr_sm[j + SM_EXTRA] += sfr_sm_tmp * smf2;
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
  int64_t bin_min=-1, bin_max=-1, bin_peak=0;
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
    
     
    // INPLEMENTAION 2: bin_max - bin_min + 1 < M_BINS

    if (steps[n].obs_uv[i]>steps[n].smhm.uv_max && steps[n].obs_uv[i] < 1000) {
      uv_max = steps[n].smhm.uv_max = steps[n].obs_uv[i];
      bin_max = i;
    }
    // if (steps[n].log_sm[i]>-1000 && bin_min < 0) {
    if (steps[n].obs_uv[i]<steps[n].smhm.uv_min) {
      bin_min = i;
      // if (n == 66) fprintf(stderr, "n=66, bin_min=%d\n", bin_min);
      uv_min = steps[n].smhm.uv_min = steps[n].obs_uv[i];
    }

    if (i && (steps[n].obs_uv[i] >= steps[n].obs_uv[i-1])) {
	    //fprintf(stderr, "n=%d, z=%f, m1=%f, m2=%f, uv1=%f, uv2=%f\n", n, 1.0/steps[n].scale - 1, M_MIN + (i + 0.5) * INV_BPDEX, M_MIN + (i - 0.5) * INV_BPDEX, steps[n].obs_uv[i], steps[n].obs_uv[i - 1]);
      // if (steps[n].log_sm[i]==0) { steps[n].log_sm[i]=0; } // I totally don't understand why you would do this.
      //                                                           // In this case at the very massive end of the SMF,
      //                                                           // the corresponding halo mass would be very difficult
      //                                                           // to determine, and cause the crazy behavior there.
      //                                                           // For now I just put an empty block here.
      if (steps[n].obs_uv[i]==10000) 
        { 
          // if (bin_min >= 0) steps[n].log_sm[i] = steps[n].log_sm[i-1]+0.001; 
          // steps[n].log_sm[i] = steps[n].log_sm[i-1]+0.001; 
          // if (i > bin_max) bin_max = i; //Even if we don't need this artificial massive end, we still
          //                                                       // need to maintain the bin_max to keep the size consistency between
          //                                                       // the interpolated data and the gsl_spline object.
          if (M_MIN+i*INV_BPDEX < 12) 
          {
            // steps[n].obs_uv[i] = steps[n].obs_uv[i-1]-0.001;
            ;
          } 
          else //at the high mass end, we do a constant SM/HM ratio extrapolation.
          {
            steps[n].obs_uv[i] = steps[n].obs_uv[i-1]-INV_BPDEX;
          }
          if (steps[n].obs_uv[i] < steps[n].smhm.uv_min) 
	  {
            uv_min = steps[n].smhm.uv_min = steps[n].obs_uv[i];
            bin_min = i;
	  }
        }
      else {
        if (!count_falling)
        {
          count_falling = 1;
          bin_peak = i - 1;
        }
        
  // // sprintf(buffer, "Falling SMHM relation: SM(%f) = %f; SM(%f) = %f at scale %f (Mass at z=0: %f)\n", steps[n].med_hm_at_a[i], steps[n].log_sm[i], steps[n].med_hm_at_a[i-1], steps[n].log_sm[i-1], steps[n].scale, i*INV_BPDEX+M_MIN);       
  // fprintf(stderr, "Falling UVHM relation: UV(%f) = %f; UV(%f) = %f at scale %f (Mass at z=0: %f)\n", steps[n].med_hm_at_a[i], steps[n].obs_uv[i], steps[n].med_hm_at_a[i-1], steps[n].obs_uv[i-1], steps[n].scale, i*INV_BPDEX+M_MIN);       
  // INVALIDATE(fit, buffer);
  // steps[n].smhm.valid = 0;
  // return;
      }
    }

  }
  // fprintf(stderr, "count_falling=%d, bin_peak=%d, bin_max=%d\n", count_falling, bin_peak, bin_max); 
 // printf("bin_min=%d, bin_max=%d\n", bin_min, bin_max);
  //fprintf(stderr, "n=%d, bin_min=%d\n", n, bin_min);
    if (bin_max < 0 || bin_min - bin_max + 1 < 10 || (count_falling && bin_peak <= bin_max + 1))
{
        //sprintf(buffer, "All the stellar masses are zero.\n");
        //fprintf(stderr, "bin_max=%d, bin_min=%d\n", bin_max, bin_min);
	//fprintf(stderr, "All the UV magnitudes are zero at z=%f (%d-th snapshot). bin_min=%d, bin_max=%d\n", 1.0 / steps[n].scale - 1, n, bin_min, bin_max);
        if (1.0 / steps[n].scale - 1 >= 7.5) INVALIDATE(fit, buffer);
        return;
}

  
  double *obs_uv_tmp = NULL;
  // double *std_uv_tmp = NULL;
  double  *hm_tmp = NULL;

  if (!count_falling)
  {
    obs_uv_tmp = malloc(sizeof(double)*(bin_min - bin_max + 1));
    // std_uv_tmp = malloc(sizeof(double)*(bin_min - bin_max + 1));
    hm_tmp = malloc(sizeof(double)*(bin_min - bin_max + 1));

    for (j=bin_max; j <= bin_min; j++)
    {
      // printf("log_sm=%f\n", steps[n].log_sm[j]);
      obs_uv_tmp[j - bin_max] = steps[n].obs_uv[bin_min - (j - bin_max)];
      // std_uv_tmp[j - bin_max] = steps[n].std_uv[bin_min - (j - bin_max)];
      hm_tmp[j - bin_max] = steps[n].med_hm_at_a[bin_min - (j - bin_max)];

    }

    
    // The spline does not need to be modified. What should be changed is that the dn/dlogm
    // should be weighted by sfrac and 1 - sfrac for smf_sf and smf_q, respectively.
    steps[n].spline_uv = gsl_spline_alloc(gsl_interp_cspline, bin_min - bin_max + 1);
    // steps[n].spline_std_uv = gsl_spline_alloc(gsl_interp_cspline, bin_max - bin_min + 1);



    steps[n].flag_alloc = 2;
    gsl_spline_init(steps[n].spline_uv, obs_uv_tmp, hm_tmp, bin_min - bin_max + 1);
    // gsl_spline_init(steps[n].spline_std_uv, obs_uv_tmp, std_uv_tmp, bin_max - bin_min + 1);
    free(obs_uv_tmp); free(hm_tmp);
  }

  else
  {
    obs_uv_tmp = malloc(sizeof(double)*(bin_peak - bin_max + 1));
    // std_uv_tmp = malloc(sizeof(double)*(bin_min - bin_max + 1));
    hm_tmp = malloc(sizeof(double)*(bin_peak - bin_max + 1));

    for (j=bin_max; j <= bin_peak; j++)
    {
      // printf("log_sm=%f\n", steps[n].log_sm[j]);
      obs_uv_tmp[j - bin_max] = steps[n].obs_uv[bin_peak - (j - bin_max)];
      // std_uv_tmp[j - bin_max] = steps[n].std_uv[bin_peak - (j - bin_max)];
      hm_tmp[j - bin_max] = steps[n].med_hm_at_a[bin_peak - (j - bin_max)];
      // fprintf(stderr, "obs_uv[%d]=%f\n", j - bin_max, obs_uv_tmp[j - bin_max]);
    }

    
    // The spline does not need to be modified. What should be changed is that the dn/dlogm
    // should be weighted by sfrac and 1 - sfrac for smf_sf and smf_q, respectively.
    steps[n].spline_uv = gsl_spline_alloc(gsl_interp_cspline, bin_peak - bin_max + 1);
    // steps[n].spline_std_uv = gsl_spline_alloc(gsl_interp_cspline, bin_max - bin_min + 1);
    gsl_spline_init(steps[n].spline_uv, obs_uv_tmp, hm_tmp, bin_peak - bin_max + 1);
    // gsl_spline_init(steps[n].spline_std_uv, obs_uv_tmp, std_uv_tmp, bin_max - bin_min + 1);
    
    free(obs_uv_tmp); free(hm_tmp);

    // Here we need to account for the fact that the segment of UVHM relation after
    // the turning over point may be fewer than 2 points. In this case
    // cubic spline interpolation won't work. And my solution is to do a linear
    // extrapolation out to 10^16.3 Msun to get the third point.
    int n_seg = M_BINS - bin_peak > 2 ? M_BINS - bin_peak : 3;

    obs_uv_tmp = malloc(sizeof(double)*n_seg);
    // std_uv_tmp = malloc(std_uv_tmp, sizeof(double)*(bin_min - bin_max + 1));
    hm_tmp = malloc(sizeof(double)*n_seg);

    // obs_uv_tmp = malloc(sizeof(double)*(M_BINS - bin_peak));
    // // std_uv_tmp = malloc(std_uv_tmp, sizeof(double)*(bin_min - bin_max + 1));
    // hm_tmp = malloc(sizeof(double)*(M_BINS - bin_peak));

    // obs_uv_tmp = (double *)realloc(obs_uv_tmp, sizeof(double)*(M_BINS - bin_peak));
    // // std_uv_tmp = malloc(std_uv_tmp, sizeof(double)*(bin_min - bin_max + 1));
    // hm_tmp = (double *)realloc(hm_tmp, sizeof(double)*(M_BINS - bin_peak));

    for (j=bin_peak; j < M_BINS; j++)
    {
      // printf("log_sm=%f\n", steps[n].log_sm[j]);
      obs_uv_tmp[j - bin_peak] = steps[n].obs_uv[j];
      // std_uv_tmp[j - bin_peak] = steps[n].std_uv[j];
      hm_tmp[j - bin_peak] = steps[n].med_hm_at_a[j];
      // fprintf(stderr, "bin_peak=%d, bin_min=%d, obs_uv[%d]=%f\n", bin_peak, bin_min, j - bin_peak, obs_uv_tmp[j - bin_peak]);
    }

    // Linear extrapolation.
    if (n_seg > M_BINS - bin_peak)
    {
      // fprintf(stderr, "n_seg=%d, M_BINS - bin_peak=%d\n", n_seg, M_BINS - bin_peak);
	    hm_tmp[2] = 16.3; obs_uv_tmp[2] = 2 * obs_uv_tmp[1] - obs_uv_tmp[0];
    }

    // The spline does not need to be modified. What should be changed is that the dn/dlogm
    // should be weighted by sfrac and 1 - sfrac for smf_sf and smf_q, respectively.
    steps[n].spline_uv2 = gsl_spline_alloc(gsl_interp_cspline, n_seg);
    // steps[n].spline_std_uv = gsl_spline_alloc(gsl_interp_cspline, bin_max - bin_min + 1);
    int err_spline_init = gsl_spline_init(steps[n].spline_uv2, obs_uv_tmp, hm_tmp, n_seg);
    // gsl_spline_init(steps[n].spline_std_uv, obs_uv_tmp, std_uv_tmp, bin_max - bin_min + 1);

    
    // // The spline does not need to be modified. What should be changed is that the dn/dlogm
    // // should be weighted by sfrac and 1 - sfrac for smf_sf and smf_q, respectively.
    // steps[n].spline_uv2 = gsl_spline_alloc(gsl_interp_cspline, M_BINS - bin_peak);
    // // steps[n].spline_std_uv = gsl_spline_alloc(gsl_interp_cspline, bin_max - bin_min + 1);
    // int err_spline_init = gsl_spline_init(steps[n].spline_uv2, obs_uv_tmp, hm_tmp, M_BINS - bin_peak);
    // // gsl_spline_init(steps[n].spline_std_uv, obs_uv_tmp, std_uv_tmp, bin_max - bin_min + 1);
    free(obs_uv_tmp); free(hm_tmp);
    if (err_spline_init)
    {
      //fprintf(stderr, "More than 1 turning point in the UVHM at %d-th snapshot.\n", n);
      if (1.0/steps[n].scale - 1 > 7.5) INVALIDATE(fit, buffer);
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

      int err = gsl_spline_eval_e(steps[n].spline_uv, uv, &ga, &m);
      // fprintf(stdout, "n=%d, z=%f, uv=%f, interpolated m=%f, err=%d\n", n, 1.0 / steps[n].scale - 1, uv, m, err);
      if (err)
      //if (err || m < M_MIN - 1.0 || m > M_MAX + 1.0)
      {
        // sprintf(buffer, "Error in GSL spline interpolation #2.\n");
        //fprintf(stderr, "Error in GSL spline interpolation #2. m=%f\n", m);
        INVALIDATE(fit, buffer);
        // free(obs_uv_tmp); free(hm_tmp);
        return;
      }
      //int64_t b = gsl_interp_accel_find(&ga, steps[n].log_sm+bin_min, (bin_max-bin_min+1), sm)+bin_min;
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
      steps[n].uvlf[j+SM_EXTRA] = calc_uvlf_at_m(m,uv,n,fit,&ga);

      if (steps[n].alloc2uvlf)
      {
        gsl_interp_accel_reset(&ga);
        int err2 = gsl_spline_eval_e(steps[n].spline_uv2, uv, &ga, &m2);
        if (!err2)
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
          steps[n].uvlf[j+SM_EXTRA] += calc_uvlf_at_m(m2,uv,n,fit,&ga);
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
    //if (n == 37) fprintf(stderr, "uv=%f, uvlf[%d]=%e\n", uv, j, steps[n].uvlf[j+UV_EXTRA]);
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
  double icl_frac, ejec_frac;
  steps[n].smhm = smhm_at_z((1.0/steps[n].scale)-1.0, *fit);
  steps[n].flag_alloc = 0;
  double scatter_corr = steps[n].smhm.scatter_corr; 
  // double bh_sm_scatter = steps[n].smhm.scatter*steps[n].smhm.bh_gamma;
  // double bh_vdisp_scatter = 0.3; //for the time being, use 0.3 dex as the vdisp scatter
  // double combined_scatter = (vel_dispersion) ? sqrt(bh_vdisp_scatter  * bh_vdisp_scatter + steps[n].smhm.bh_scatter*steps[n].smhm.bh_scatter) : sqrt(bh_sm_scatter  * bh_sm_scatter + steps[n].smhm.bh_scatter*steps[n].smhm.bh_scatter);
  // double bh_scatter_corr = exp(pow(combined_scatter*log(10), 2)/2.0);
  double a200 = steps[n].scale/0.3782;
  double m200 = log10(1.115e12/0.68 / (pow(a200, -0.142441) + pow(a200, -1.78959)));
  double v200 = log10(200);


  // double bh_sm_scatter = steps[n].smhm.scatter*steps[n].smhm.bh_gamma;
  // double bh_vdisp_scatter = 0.3; //for the time being, use 0.3 dex as the vdisp scatter
  // double combined_scatter = (vel_dispersion) ? sqrt(bh_vdisp_scatter  * bh_vdisp_scatter + steps[n].smhm.bh_scatter*steps[n].smhm.bh_scatter) : sqrt(bh_sm_scatter  * bh_sm_scatter + steps[n].smhm.bh_scatter*steps[n].smhm.bh_scatter);

  
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
    // correct for that at this position. Instead we have to do this after calculating the 
    // old_sm.
    // This shouldn't affect anything, because the sm[i] calculated in line 962 is useful
    // only for generating fake SFH, where all the galaxies should be SF-ing.
    // steps[n].sfrac[i] = calc_sfrac_at_lv(lv, steps[n].smhm);
    steps[n].sfrac[i] = calc_sfrac_at_m(steps[n].med_hm_at_a[i], steps[n].smhm);
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
      if (n>0) {
	sm_inc = sm_mmp = 0;
	for (j=0; j<M_BINS; j++) {
	  bin = j*M_BINS + i;
    // steps[n].old_bh_mass[i] += steps[n-1].bh_mass_avg[j]*steps[n].mmp[bin];
    // if (j == i - 1) steps[n].old_bh_mass_next[i] = steps[n-1].bh_mass_avg[j]*steps[n].mmp[bin];
    // else if (j == i) steps[n].old_bh_mass_self[i] = steps[n-1].bh_mass_avg[j]*steps[n].mmp[bin];
	  // steps[n].bh_unmerged[i] += steps[n].merged[bin]*(steps[n-1].bh_mass_avg[j] + 
                                                            // steps[n-1].bh_unmerged[j]) + steps[n].mmp[bin]*steps[n-1].bh_unmerged[j];
    
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
	  }
	}
      }
      for (k=0; k<n; k++) sm_hist_i[k] /= steps[n].t[i];
      for (k=0; k<n; k++) icl_hist_i[k] /= steps[n].t[i];
      steps[n].sm_from_icl[i] /= steps[n].t[i];
      steps[n].sm_inc[i] /= steps[n].t[i];
      // steps[n].old_bh_mass[i] /= steps[n].t[i];
      // steps[n].old_bh_mass_self[i] /= steps[n].t[i];
      // steps[n].old_bh_mass_next[i] /= steps[n].t[i];
      // steps[n].bh_unmerged[i] /= steps[n].t[i];
      // steps[n].bh_unmerged_avail[i] = steps[n].bh_unmerged[i];
      //for (k=0; k<MBH_BINS; k++) steps[n].bh_unmerged_dist /= steps[n].t[i]; // Don't forget the normalization.
    }
    // smooth_bins(steps[n].bh_unmerged_dist + i*MBH_BINS, steps[n].bh_unmerged_dist + i*MBH_BINS,
    //             combined_scatter, MBH_MIN, MBH_BPDEX, MBH_BINS,
    //                               MBH_MIN, MBH_BPDEX, MBH_BINS, 0, 0);

    calc_old_sm(n,i);
    // correct the SFR by accounting for the contribution from quenched galaxies.
    //double ratio = steps[n].t[i] ? steps[n].n[i] / steps[n].t[i] : 1;
    steps[n].sfr[i] += (1 - steps[n].sfrac[i]) * steps[n].old_sm[i] * 1.585E-12 * exp(pow(0.42*log(10), 2)/2.0);
    if (n == 0) steps[n].sfr[i] = 0;


    // Calculate observed UV magnitudes.
    double z, a1, mh, k_uv, b_uv, k_std_uv, b_std_uv, a_k, b_k, c_k, a_b, b_b, c_b;
    double a_k_std, b_k_std, a_b_std, b_b_std;
    z = 1.0 / steps[n].scale - 1;
    a1 = steps[n].scale - 1;
    a_k = 0.15367887; b_k = -2.87608013; c_k = 9.4778494 + (-2.37851257) * a1;
    a_b = -0.34749622; b_b = 6.85270974; c_b = -50.34457998 + 1.99341162 * a1;

    mh = steps[n].med_hm_at_a[i];
    

    k_uv = a_k * mh * mh + b_k * mh;
    b_uv = a_b * mh * mh + b_b * mh;
    k_uv = mh < - 0.5 * b_k / a_k ? -0.25 * b_k * b_k / a_k + c_k : k_uv + c_k;
    b_uv = mh < - 0.5 * b_b / a_b ? -0.25 * b_b * b_b / a_b + c_b : b_uv + c_b;
    k_std_uv = 0.04724224 * mh + (-0.17632334) + 0.49534345 * a1;
    b_std_uv = 0.03879428 * mh + (-0.53377296) + (-0.50621441) * a1;

    double lgSFR = log10(steps[n].sfr[i]);
    steps[n].obs_uv[i] = k_uv * lgSFR + b_uv;
    steps[n].std_uv[i] = k_std_uv * lgSFR + b_std_uv;
    if (steps[n].obs_uv[i] > 0) steps[n].obs_uv[i] = 10000;
    if (steps[n].std_uv[i] < 0) steps[n].std_uv[i] = 0.001; 
    if (steps[n].std_uv[i] > 1) steps[n].std_uv[i] = 1.000; 




    //if (i == 0) fprintf(stderr, "z=%f, sfr=%e, after adding QFGs.\n", 1 / steps[n].scale - 1, steps[n].sfr[i]);
    calc_new_sm_and_sfr(n,i,fit);
    // // Moved to the calc_new_sm_and_sfr().
    // steps[n].log_bm[i] = bulge_mass(steps[n].log_sm[i]+steps[n].smhm.mu, steps[n].scale);
    // steps[n].log_bh_mass[i] = (vel_dispersion) ? 
    //     calc_bh_at_vdisp(steps[n].vdisp[i], steps[n].smhm)
    //   : calc_bh_at_bm(steps[n].log_bm[i], steps[n].smhm);
    // steps[n].bh_mass_avg[i] = doexp10(steps[n].log_bh_mass[i])*bh_scatter_corr;


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

double recent_sfh_in_massive_halos(void) {
  int64_t n, j, count=0;
  double sfr = 0;
  for (n=0; n<num_outputs; n++) if (steps[n].scale > SFR_A_CONSTRAINT) break;
  for (; n<num_outputs; n++) {
    double corr = doexp10(steps[n].smhm.mu);
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
  //if (n == 0) fprintf(stderr, "before ICL addition, new_sm[%d]=%e\n", i, steps[n].new_sm[i]);
  // //  steps[n].new_sm[i] = (steps[n].sm_avg[i] - steps[n].old_sm[i]);
  // float dm = (i+0.5)*INV_BPDEX+M_MIN-steps[n].smhm.icl_m;
  // float w = 16.0-steps[n].smhm.icl_m;
  // //if (steps[n].smhm.beta>0) { w = 4.0 / (steps[n].smhm.beta * M_LN10); }
  // if (w < 0.1) w = 0.1;
  // steps[n].icl_frac[i] = 1.0 - 1.0/(1.0+exp(4.0*dm/w));
  // if (steps[n].icl_frac[i] > 0.99) steps[n].icl_frac[i] = 0.99;
  double sm_from_icl = steps[n].smhm.icl_frac * steps[n].sm_icl[i];

  // steps[n].new_sm[i] *= 1.0/(1.0 - steps[n].icl_frac[i]);
  steps[n].new_sm[i] += sm_from_icl;
  //if (n == 0) fprintf(stderr, "after ICL addition, new_sm[%d]=%e\n", i, steps[n].new_sm[i]);
  steps[n].sm_avg[i] = steps[n].old_sm[i] + steps[n].new_sm[i];
  steps[n].sm[i] = steps[n].sm_avg[i] / steps[n].smhm.scatter_corr;
  steps[n].log_sm[i] = (steps[n].sm[i] > 1) ? log10(steps[n].sm[i]) : 0;

  steps[n].log_bm[i] = bulge_mass(steps[n].log_sm[i]+steps[n].smhm.mu, steps[n].scale);
  steps[n].log_sm_obs[i] = steps[n].log_sm[i]+steps[n].smhm.mu;
 
  // // We no longer need this because we're always taking away a certain fraction (<1)
  // // of ICL mass and feed them to the galaxy. In this way ICL will never get depleted.
  // if (steps[n].sm_icl[i] < steps[n].icl_frac[i]*steps[n].new_sm[i]) 
  // {
  //   if (steps[n].new_sm[i] > 0 && steps[n].scale > 0.15 && i>BPDEX) 
  //   {
  //     // sprintf(buffer, "ICL depleted (hm: %f; a: %f; sm: %e, new sm: %e, icl needed: %e, icl available: %e)\n",(i+0.5)*INV_BPDEX+M_MIN, steps[n].scale,  steps[n].sm_avg[i], steps[n].new_sm[i], steps[n].new_sm[i] * steps[n].icl_frac[i], steps[n].sm_icl[i]);
  //     //fprintf(stderr, "ICL depleted (hm: %f; a: %f; sm: %e, new sm: %e, icl needed: %e, icl available: %e)\n",(i+0.5)*INV_BPDEX+M_MIN, steps[n].scale,  steps[n].sm_avg[i], steps[n].new_sm[i], steps[n].new_sm[i] * steps[n].icl_frac[i], steps[n].sm_icl[i]);
  //     INVALIDATE(fit, buffer);
  //     steps[n].icl_frac[i] = steps[n].sm_icl[i]/steps[n].new_sm[i];
  //   }
  //   else 
  //     steps[n].icl_frac[i] = 0;
  // }

  steps[n].merged_frac[i] = 0;

  // if (steps[n].icl_frac[i] && steps[n].sm_icl[i]) 
  // {
  //   double frac_from_icl = steps[n].new_sm[i]*steps[n].icl_frac[i] / steps[n].sm_icl[i];
  //   if (steps[n].sm_inc[i]>0) 
  //   {
  //     steps[n].merged_frac[i] = steps[n].new_sm[i]*steps[n].icl_frac[i]/steps[n].sm_inc[i];
  //     if (steps[n].merged_frac[i]>1) steps[n].merged_frac[i] = 1;
  //   }
  //   if (frac_from_icl) 
  //   {
  //     for (j=0; j<n; j++) 
  //     {
  //       steps[n].sm_hist[i*num_outputs + j] += frac_from_icl*steps[n].icl_stars[i*num_outputs + j];
  //       steps[n].icl_stars[i*num_outputs + j] *= (1.0-frac_from_icl);
  //     }
  //   }
  //   steps[n].new_sm[i] *= (1.0-steps[n].icl_frac[i]);
  // }


  if (steps[n].smhm.icl_frac && steps[n].sm_icl[i]) 
  {
    double frac_from_icl = steps[n].smhm.icl_frac; //Since steps[n].icl_frac is the fraction
                                              //of the ICL mass that got into the descendant,
                                              //***NOT*** the fraction of the new SM that come
                                              //from mergers.
    if (steps[n].sm_inc[i]>0) 
    {
      steps[n].merged_frac[i] = steps[n].sm_icl[i]*frac_from_icl/steps[n].sm_inc[i];
      if (steps[n].merged_frac[i]>1) steps[n].merged_frac[i] = 1;
    }
    if (frac_from_icl) 
    {
      for (j=0; j<n; j++) 
      {
        steps[n].sm_hist[i*num_outputs + j] += frac_from_icl*steps[n].icl_stars[i*num_outputs + j];
        steps[n].icl_stars[i*num_outputs + j] *= (1.0-frac_from_icl);
      }
    }
    steps[n].new_sm[i] -= (frac_from_icl * steps[n].sm_icl[i]);
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

