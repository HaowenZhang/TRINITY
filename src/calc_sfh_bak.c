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
#include <omp.h>
#include <assert.h>

#define CALC_ICL 

//#undef exp10fc
//#define exp10fc(x) doexp10(x)
//#define log10fc(x) dolog10(x)


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
  int64_t i,j;
  for (i=0; i<num_outputs; i++)
    calc_sm_hist(i, f);

#pragma omp for schedule(guided,5)
  for (j=0; j<num_outputs; j++) {
    calc_smf_and_ssfr(j, f);
    calc_total_sfr(j);
  }
  /*if (omp_get_thread_num() < 1) {
    if (INVALID(*f)) fprintf(stderr, "X");
    else fprintf(stderr, ".");
    }*/
}

void calc_total_sfr(int n) {
  int64_t i;
  double sfr = 0, obs_sfr = 0;
  //double kappa = steps[n].smhm.kappa;
  double mu = steps[n].smhm.mu;
  double sm_corr, sm_meas;
  for (i=0; i<M_BINS; i++) {
    sfr += steps[n].sfr[i]*steps[n].t[i];
    //sm_meas = (log10fc(steps[n].sm[i])*(1.0+kappa) + mu - kappa*SLOPE_SYSTEMATIC_MASS);
    sm_meas = (log10(steps[n].sm[i]) + mu);
    sm_corr = (steps[n].sm[i] > 0) ? doexp10(sm_meas)/steps[n].sm[i] : 1;
    obs_sfr += steps[n].sfr[i]*steps[n].t[i]*sm_corr;
  }
  steps[n].cosmic_sfr = sfr;
  steps[n].observed_cosmic_sfr = obs_sfr*steps[n].smhm.csfr_completeness;
}

struct smf smhm_at_z(double z, struct smf_fit f) {
  struct smf c;
  double a = 1.0/(1.0+z);
  double a1 = -z/(1.0+z);
  double expscale = doexp10(-EXPSCALE(f)*(1+a1)*(1+a1)*M_LOG10E);
  double incompleteness;
  double a4 = z; //a1*a1*a1*a1;
  //double lgz = log(1.0+z);

  /*  c.m_1 = M_1(f) + (a1*M_1_A(f)+z*M_1_A2(f))*expscale;
  c.sm_0 = c.m_1 + EFF_0(f) + (EFF_0_A(f)*a1 + EFF_0_A2(f)*z)*expscale + EFF_0_A3(f)*a1;
  c.alpha = ALPHA(f) + (a1*ALPHA_A(f) + z*ALPHA_A2(f))*expscale;
  c.delta = DELTA(f) + (a1*DELTA_A(f) + z*DELTA_A2(f))*expscale;
  c.beta = BETA(f) + (a1*BETA_A(f) + z*BETA_A2(f))*expscale;
  c.gamma = GAMMA(f) + (a1*GAMMA_A(f) + z*GAMMA_A2(f))*expscale;
  c.lambda = LAMBDA(f) + (a1*LAMBDA_A(f) + z*LAMBDA_A2(f))*expscale;*/

  c.m_1 = M_1(f) + (a1*M_1_A(f)+z*M_1_A2(f))*expscale;
  c.sm_0 = c.m_1 + EFF_0(f) + (EFF_0_A(f)*a1 + EFF_0_A2(f)*z)*expscale + EFF_0_A3(f)*a1;
  c.alpha = ALPHA(f) + (a1*ALPHA_A(f) + a4*ALPHA_A2(f))*expscale;
  c.delta = DELTA(f) + (a1*DELTA_A(f) + z*DELTA_A2(f))*expscale;
  c.beta = BETA(f) + (a1*BETA_A(f) + z*BETA_A2(f))*expscale;
  c.gamma = GAMMA(f) + (a1*GAMMA_A(f) + z*GAMMA_A2(f))*expscale;
  c.lambda = LAMBDA(f) + (a1*LAMBDA_A(f) + z*LAMBDA_A2(f))*expscale;
  //if (c.gamma > 1) c.gamma = 1;
  c.mu = MU(f) + a1*MU_A(f);
  c.kappa = KAPPA(f) + a1*KAPPA_A(f);
  c.passive_mass = 10.3 + z*0.5 - c.mu;
  c.scatter = SCATTER(f) + a1*SCATTER_A(f);
  c.obs_scatter = (no_obs_scatter) ? 0 : SIGMA_CENTER + SIGMA_Z(f)*(z-0.1);
  c.icl_frac = exp10(ICL_FRAC(f));
  c.icl_frac_e = ICL_FRAC_E(f);
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
  return c;
}

/*double calc_sm_at_m(double m, struct smf c) {
  double dm = m-c.m_1;
  double dm2 = dm + c.lambda;
  double sm, dma = dm*c.alpha, dmb = dm*c.beta;
  //if (dma > dmb+20) sm = -dma;
  //else if (dmb > dma+20) sm = -dmb;
  //else sm = -dolog10(doexp10(dma) + doexp10(dmb));
  sm = -dolog10(doexp10(dma) + doexp10(dmb));
  sm += c.sm_0 + c.delta*doexp10(-c.gamma*(dm2*dm2)) + (M_LN2/M_LN10);
  return sm;
  }*/

double calc_sm_at_m(double m, struct smf c) {
  double dm = m-c.m_1;
  double sm, dma = dm*c.alpha;
  //if (dma > dmb+20) sm = -dma;
  //else if (dmb > dma+20) sm = -dmb;
  //else sm = -dolog10(doexp10(dma) + doexp10(dmb));
  sm = -log10_1p_exp10fc(dma);
  //sm += c.sm_0 + c.delta*(exp10(log10fc(log10_1p_exp10fc(dm))*c.gamma)/(1.0+exp(exp10(-dm)))) + c.f_1;
  sm += c.sm_0 + c.delta*(exp10(dolog10(log10(1.0+exp(dm)))*c.gamma)/(1.0+exp(exp10(-dm)))) + c.f_1;
  return sm;
}

/*double calc_smf_at_m(double m, double scale, struct smf c, struct smf_fit *fit){
  double sm1 = calc_sm_at_m(m, c);
  double sm2 = calc_sm_at_m(m+0.001, c);
  double dm = sm2 - sm1;
  if (sm2 <= sm1 && fit) {
    INVALIDATE(fit, "Falling SM(M) relation.");
    //return 0;
    dm = 0.001;
  }
  double dlm_dlsm = 0.001/(sm2-sm1);
  double dndlogm = mf_cache(scale, m);
  double phi = doexp10(dndlogm)*dlm_dlsm;
  if (!isfinite(phi)) phi = 0;
  return phi;
  }

double calc_m_at_sm(double sm, double m, struct smf c) {
  double sm2 = calc_sm_at_m(m, c);
  double dm_max = 5*INV_BPDEX;
  double dm = dm_max/10.0;
  if (m<=0) {
    m = M_MIN;
    dm_max = 10;
  }
  int64_t count = 0;
  while (fabs(sm2-sm)>0.0001 && count<30) {
    double sm3 = calc_sm_at_m(m+dm,c);
    double f = (sm - sm2)/(sm3-sm2);
    double move_m = f*dm;
    if (fabs(move_m) > dm_max)
      return m;
    dm = move_m / 10.0;
    m+=move_m;
    sm2 = calc_sm_at_m(m,c);
    count++;
  }
  if (count>=25) return -1;
  return m;
  }*/

/*double calc_smf_at_sm(int64_t n, int64_t j, double sm) {
  double m = calc_m_at_sm(sm, steps[n].m_at_sm[j], steps[n].smhm);
  return calc_smf_at_m(m, steps[n].scale, steps[n].smhm, NULL);
  }*/

double calc_smf_at_m(double m, double sm, int64_t n, struct smf_fit *fit, gsl_interp_accel *ga) {
  if (sm < steps[n].smhm.sm_min) return 0;
  if (sm > steps[n].smhm.sm_max) return 0;
  double dlm_dlsm = gsl_spline_eval_deriv(steps[n].spline, sm, ga);
  if (!(dlm_dlsm>0) && fit) {
    INVALIDATE(fit, "Falling SM(M) relation.");
    dlm_dlsm = 1.0;
  }
  double dndlogm = mf_cache(steps[n].scale, m);
  double phi = doexp10(dndlogm)*dlm_dlsm;
  if (!isfinite(phi)) phi = 0;
  return phi;
}

double calc_smf_at_sm(int64_t n, double sm) {
  if (sm < steps[n].smhm.sm_min) return 0;
  if (sm > steps[n].smhm.sm_max) return 0;
  gsl_interp_accel ga = {0};
  double m = gsl_spline_eval(steps[n].spline, sm, &ga);
  return calc_smf_at_m(m, sm, n, NULL, &ga);
}

void calc_smf_and_ssfr(int n, struct smf_fit *fit) {
  int64_t i, j;
  char buffer[1024];
  double m, sm, sm_max=0, sm_min=1000;
  int64_t bin_min=-1, bin_max=-1;
  gsl_interp_accel ga = {0};
  steps[n].smhm.sm_max = 0;
  if (INVALID(*fit)) return;
  //if (steps[n].scale < 0.1) return;

  for (i=0; i<M_BINS; i++) {
    if (steps[n].log_sm[i]>steps[n].smhm.sm_max) {
      sm_max = steps[n].smhm.sm_max = steps[n].log_sm[i];
      bin_max = i;
    }
    if (steps[n].log_sm[i]>-1000 && bin_min < 0) {
      bin_min = i;
      sm_min = steps[n].smhm.sm_min = steps[n].log_sm[i];
    }
    if (i && (steps[n].log_sm[i] < steps[n].log_sm[i-1])) {
      sprintf(buffer, "Falling SMHM relation: SM(%f) = %f; SM(%f) = %f at scale %f (Mass at z=0: %f)\n", steps[n].med_hm_at_a[i], steps[n].log_sm[i], steps[n].med_hm_at_a[i-1], steps[n].log_sm[i-1], steps[n].scale, i*INV_BPDEX+M_MIN);	      
      INVALIDATE(fit, buffer);
      return;
    }
  }
  
  if (!((bin_max-bin_min+1) == M_BINS)) {
    sprintf(buffer, "Incorrect number of bins to generate spline at scale %f\n", steps[n].scale);
    INVALIDATE(fit, buffer);
    return;
  }

  gsl_spline_init(steps[n].spline, steps[n].log_sm+bin_min, steps[n].med_hm_at_a+bin_min, bin_max - bin_min + 1);

  i = M_BINS-1;
  for (j=SM_BINS-1; j>=0; j--) {
    sm = SM_MIN + (double)j*SM_INV_BPDEX;
    if (sm > sm_max || sm < sm_min) {
      steps[n].sfr_sm[j+SM_EXTRA] = 0;
      steps[n].smf[j+SM_EXTRA] = 0;
    }
    else {
      m = gsl_spline_eval(steps[n].spline, sm, &ga);
      int64_t b = gsl_interp_accel_find(&ga, steps[n].log_sm+bin_min, (bin_max-bin_min+1), sm)+bin_min;
      steps[n].smf[j+SM_EXTRA] = calc_smf_at_m(m,sm,n,fit,&ga);
      steps[n].sfr_sm[j+SM_EXTRA] = steps[n].sfr[b];
      if (b<bin_max) {
	double f = (sm-steps[n].log_sm[b])/(steps[n].log_sm[b+1]-steps[n].log_sm[b]);
	steps[n].sfr_sm[j+SM_EXTRA] +=  f*(steps[n].sfr[b+1]-steps[n].sfr[b]);
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


/*void calc_smf_and_ssfr(int n, struct smf_fit *fit) {
  int64_t i, j;
  double scale = steps[n].scale;
  double m, sm, lsm1, lsm2;
  struct smf c = steps[n].smhm;
  steps[n].smhm.sm_max = 0;
  
  for (i=0; i<M_BINS; i++) {
    steps[n].sm_nd[i] = calc_smf_at_m(M_MIN + (i+0.5)*INV_BPDEX, scale, c, fit);
    if (dolog10(steps[n].sm[i])>steps[n].smhm.sm_max)
      steps[n].smhm.sm_max = dolog10(steps[n].sm[i]);
  }

  i = M_BINS-1;
  lsm1 = SM_MAX;
  for (j=SM_BINS-1; j>=0; j--) {
    sm = SM_MIN + (double)j*SM_INV_BPDEX;
    while ((steps[n].sm[i]<=0 || ((lsm1 = dolog10(steps[n].sm[i])) > sm)) && i > 0) i--;
    if (i==0 || i==M_BINS-1) {
      steps[n].sfr_sm[j+SM_EXTRA] = 0;
      steps[n].smf[j+SM_EXTRA] = 0;
      steps[n].m_at_sm[j+SM_EXTRA] = M_MAX;
    }
    else {
      lsm2 = dolog10(steps[n].sm[i+1]);
      double f = (sm - lsm1)/(lsm2-lsm1);
      m = M_MIN + (i+0.5+f)*INV_BPDEX;
      m = calc_m_at_sm(sm, M_MIN+(i+0.5)*INV_BPDEX, c);
      if (m<0 || m > M_MAX+INV_BPDEX) {
	if (steps[n].scale > 0.5 && sm > 10 && sm < steps[n].smhm.sm_max) {
	  fprintf(stderr, "#WTF??? sm: %f m: %f i: %f scale: %f; sm_max %f!!!\n", sm, m, M_MIN+(i+0.5+f)*INV_BPDEX, steps[n].scale, steps[n].smhm.sm_max);
	}
	steps[n].sfr_sm[j+SM_EXTRA] = 0;
	steps[n].smf[j+SM_EXTRA] = 0;
	steps[n].m_at_sm[j+SM_EXTRA] = M_MAX;
	continue;
      }
      f = BPDEX*(m - (M_MIN + (i+0.5)*INV_BPDEX));
      double sm2 = calc_sm_at_m(m,c);
      if (fabs(sm2-sm) > 0.002 && steps[n].scale > 0.15) {
	fprintf(stderr, "#WTF??? Sm: %f; Sm2: %f; scale: %f\n", sm, sm2, steps[n].scale);
      }
      steps[n].smf[j+SM_EXTRA] = calc_smf_at_m(m,scale,c,fit);
      steps[n].sfr_sm[j+SM_EXTRA] = steps[n].sfr[i] + f*(steps[n].sfr[i+1]-steps[n].sfr[i]);
      steps[n].m_at_sm[j+SM_EXTRA] = m;
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
  } */

void create_fake_sfr_hist(int n, int i) {
  double frac_lost = 0, weight = 0;
  double sm, total_time, time, prev_scale;
  int k;
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
}

/*inline double icl_fract(int64_t i, int64_t j, int64_t n, double ICL_RATIO) {
  double icl_frac;
  if (j>i) icl_frac = 1;
  else icl_frac = 1.0 - (2.0)/(1.0+doexp10((i-j+0.25)*INV_BPDEX*ICL_RATIO));
  if (n> 0 && steps[n-1].n[i]<0.5*steps[n].n[i]) icl_frac = 1;
  return icl_frac;
  }*/

inline double icl_fract(int64_t i, int64_t j, int64_t n, double ICL_RATIO) {
  double icl_frac;
  /*if (j>i) icl_frac = 1;
    else icl_frac = 1.0 - (2.0)/(1.0+doexp10((i-j+0.25)*INV_BPDEX*ICL_RATIO));*/
  double exponent = steps[n].icl_exp; //0.3*pow(steps[n].scale, -0.7);
  if (j>i) icl_frac = 1.0-0.16;
  else icl_frac = 1.0 - 0.16*doexp10(-1.0*(i-j+0.25)*INV_BPDEX*exponent);
  if (n> 0 && steps[n-1].n[i]<0.5*steps[n].n[i]) icl_frac = 1;
  return icl_frac;
}


void calc_sm_hist(int n, struct smf_fit *fit) {
  int64_t i, j, k, bin;
  double t, ICL_FRACTION_E = ICL_FRAC_E(*fit);
  double *sm_hist_i, *sm_hist_j, sm_inc, sm_mmp;
  double ICL_RATIO = doexp10(ICL_FRAC(*fit));
#ifdef CALC_ICL
  double *icl_hist_i, *icl_hist_j;
#endif /*CALC_ICL*/
  double icl_frac, ejec_frac;
  steps[n].smhm = smhm_at_z((1.0/steps[n].scale)-1.0, *fit);
  double scatter_corr = exp(pow(steps[n].smhm.scatter*log(10), 2)/2.0);
  
#pragma omp for schedule(dynamic,5) private(j,k,sm_hist_i,sm_hist_j,bin,t)
  for (i=0; i<M_BINS; i++) {
    steps[n].log_sm[i] = calc_sm_at_m(steps[n].med_hm_at_a[i], steps[n].smhm);
    steps[n].sm[i] = doexp10(steps[n].log_sm[i]);
    steps[n].sm_avg[i] = steps[n].sm[i]*scatter_corr;
    steps[n].sm_from_icl[i] = 0;

    if (!steps[n].n[i]) {
      if (steps[n].c[i]) {
	create_fake_sfr_hist(n, i);
      }
    }
    else {
      sm_hist_i = &(steps[n].sm_hist[i*num_outputs]);
      memset(sm_hist_i, 0, sizeof(double)*num_outputs);
#ifdef CALC_ICL
      icl_hist_i = &(steps[n].icl_stars[i*num_outputs]);
      memset(icl_hist_i, 0, sizeof(double)*num_outputs);
#endif /* CALC_ICL */
      if (n>0) {
	sm_inc = sm_mmp = 0;
	for (j=0; j<M_BINS; j++) {
	  bin = j*M_BINS + i;
	  icl_frac = icl_fract(i,j,n,ICL_RATIO);
	  if (icl_frac > 0.999) continue;
	  sm_inc += (1.0-icl_frac)*steps[n].merged[bin]*steps[n-1].sm_avg[j];
	  sm_mmp += steps[n].mmp[bin]*steps[n-1].sm_avg[j];
	}
	//	if (sm_mmp) ejec_frac = (ICL_FRACTION_E) * sm_inc / sm_mmp;
	if (sm_mmp) ejec_frac = (ICL_FRACTION_E-1.0) * sm_inc / sm_mmp;
	else ejec_frac = 0;
	if (ICL_FRACTION_E>1.0) steps[n].sm_from_icl[i] -= (ICL_FRACTION_E-1.0)*sm_inc;

	if (ejec_frac > 1) ejec_frac = 1;
	if (ejec_frac < 0) ejec_frac = 0;
	for (j=0; j<M_BINS; j++) {
	  bin = j*M_BINS + i;
	  icl_frac = icl_fract(i,j,n,ICL_RATIO);
	  //float incoming_frac = (1.0-icl_frac);
	  double incoming_frac = (1.0-icl_frac)*(1.0-ICL_FRACTION_E);
	  if (incoming_frac < 0) incoming_frac = 0;
	  t = steps[n].mmp[bin]*(1.0-ejec_frac) + 
	    incoming_frac*steps[n].merged[bin];
#ifdef  CALC_ICL
	  double iclt1 = (ejec_frac*steps[n].mmp[bin]+
			 (1.0-incoming_frac)*steps[n].merged[bin]);
	  double iclt2 = (steps[n].mmp[bin]+steps[n].merged[bin]);
	  if (!t && !iclt1 && !iclt2) continue;
	  steps[n].sm_from_icl[i] += steps[n-1].sm_from_icl[j]*steps[n].mmp[bin];
#else 
	  if (!t) continue;
#endif /*CALC_ICL*/
	  sm_hist_j = &(steps[n-1].sm_hist[j*num_outputs]);
#ifdef CALC_ICL
	  icl_hist_j = &(steps[n-1].icl_stars[j*num_outputs]);
#endif /*CALC_ICL*/
	  for (k=0; k<n; k++) {
	    sm_hist_i[k] += t*sm_hist_j[k];
#ifdef CALC_ICL
	    icl_hist_i[k] += iclt1*sm_hist_j[k];
	    icl_hist_i[k] += iclt2*icl_hist_j[k];
#endif /* CALC_ICL */
	  }
	}
      }
      for (k=0; k<n; k++) sm_hist_i[k] /= steps[n].t[i];
#ifdef CALC_ICL
      for (k=0; k<n; k++) icl_hist_i[k] /= steps[n].t[i];
      steps[n].sm_from_icl[i] /= steps[n].t[i];
#endif /* CALC_ICL */
    }

    calc_old_sm(n,i);
    calc_new_sm_and_sfr(n,i,fit);
  }
}

inline void calc_old_sm(int n, int j) {
  int k;
  steps[n].old_sm[j] = 0;
  for (k=0; k<n; k++)
    steps[n].old_sm[j] += steps[n].smloss[k]*steps[n].sm_hist[j*num_outputs+k];
#ifdef CALC_ICL
  steps[n].sm_icl[j] = 0;
  for (k=0; k<n; k++)
    steps[n].sm_icl[j] += steps[n].smloss[k]*steps[n].icl_stars[j*num_outputs+k];
#endif /* CALC_ICL */
}

double recent_sfh_in_massive_halos(void) {
  int64_t n, j, count=0;
  double sfr = 0;
  for (n=0; n<num_outputs; n++) if (steps[n].scale > SFR_A_CONSTRAINT) break;
  for (; n<num_outputs; n++) {
    double corr = doexp10(steps[n].smhm.mu);
    for (j=(SFR_M_CONSTRAINT-M_MIN)*BPDEX; j<M_BINS; j++) {
      if (steps[n].t[j]<REAL_ND_CUTOFF) continue;
      sfr += steps[num_outputs-1].sm_hist[j*num_outputs + n]*corr/steps[n].dt;
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
      if (steps[n].t[j]<REAL_ND_CUTOFF) continue;
      sfr += steps[num_outputs-1].sm_hist[j*num_outputs + n]/steps[n].dt;
      count++;
    }
  }
  sfr /= count;
  return sfr;
}



inline void calc_new_sm_and_sfr(int n, int i, struct smf_fit *fit) {
  double dt;
  dt = steps[n].dt;
  if (!steps[n].t[i]) return;
  if (steps[n].sm_avg[i] > 0.17*exp10(M_MIN+(i+0.5)*BPDEX)) {
    char buffer[1024];
    sprintf(buffer, "SM exceeds baryon fraction (sm: %e, m: %e, scale: %f!\n",
	    steps[n].sm_avg[i], exp10(M_MIN+(i+0.5)*BPDEX), steps[n].scale);
    INVALIDATE(fit, buffer);
  }

  steps[n].new_sm[i] = (steps[n].sm_avg[i] - steps[n].old_sm[i] - steps[n].sm_from_icl[i])/steps[n].smloss[n];
  if (steps[n].sm_icl[i] > steps[n].smhm.icl_frac*steps[n].sm_avg[i]) {
    double icl_frac = (1.0 - steps[n].smhm.icl_frac*steps[n].sm_avg[i]/steps[n].sm_icl[i]); // /(1.0+exp((1.0/steps[n].scale-1.0))-2.0);
    steps[n].sm_from_icl[i] += icl_frac*steps[n].new_sm[i];
    steps[n].new_sm[i] *= 1.0-icl_frac;
  }
  steps[n].sm_hist[i*num_outputs + n] = steps[n].new_sm[i];
  steps[n].sfr[i] = steps[n].new_sm[i] / dt;

  if (steps[n].t[i]<REAL_ND_CUTOFF)
    steps[n].sfr[i] = 0;

  if (steps[n].sfr[i] < 0 || !isfinite(steps[n].sfr[i])) {
    char buffer[1024];
    sprintf(buffer, "Negative SFR at a = %f and m = %f (ND=%g)! (ICL frac=%f)\n", steps[n].scale, 
	    i*INV_BPDEX + M_MIN, steps[n].t[i], steps[n].smhm.icl_frac);
    INVALIDATE(fit, buffer);
  }
}

