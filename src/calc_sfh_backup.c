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
  int64_t i,j;
  for (i=0; i<num_outputs; i++) {
    if (!no_z_scaling || ((steps[i].scale < 1.0/(z_min+1.0)) &&
			  (steps[i].scale > 1.0/(z_max+1.0))))
    calc_sm_hist(i, f);
    
  }

#pragma omp for schedule(guided,5)
  for (j=0; j<num_outputs; j++) {
    if (!no_z_scaling || ((steps[j].scale < 1.0/(z_min+1.0)) &&
			  (steps[j].scale > 1.0/(z_max+1.0)))) {
      calc_smf_and_ssfr(j, f);
      calc_total_sfr(j);
      gsl_spline_free(steps[j].spline); //free the spline object since we are allocating it every time.
      gsl_spline_free(steps[j].spline_sfr);
    }
  }
}

void calc_total_sfr(int n) {
  int64_t i;
  double sfr = 0, obs_sfr = 0;
  double mu = steps[n].smhm.mu;
  double sm_corr = pow(10, mu);
  double sfr_minimum = BOUWENS_SFR_LIMIT;
  if (steps[n].scale > 1.0/(1.0+4.0))
    sfr_minimum *= (1.0/steps[n].scale - 0.9)/4.1;
  for (i=0; i<M_BINS; i++) {
    sfr += steps[n].sfr[i]*steps[n].t[i];
    if (steps[n].sfr[i]<sfr_minimum) continue;
    obs_sfr += steps[n].sfr[i]*steps[n].t[i]*sm_corr;
  }
  steps[n].cosmic_sfr = sfr;
  steps[n].observed_cosmic_sfr = obs_sfr*steps[n].smhm.csfr_completeness;
}

struct smf smhm_at_z(double z, struct smf_fit f) {
  struct smf c;
  double a = 1.0/(1.0+z);
  double a1 = -z/(1.0+z);
  //  double expscale = 1.0; //doexp10(-EXPSCALE(f)*(a*a)*M_LOG10E);
  double incompleteness;
  double a2 = a1*a1;

  c.icl_m = ICL_FRAC(f)+a1*ICL_FRAC_E(f);
  c.v_1 = M_1(f) + (a1*M_1_A(f)+z*M_1_A2(f)) + a2*M_1_A3(f);
  double z8 = z;
  if (z8 > 12) z8=12;
  c.sm_0 = 0;
  c.epsilon = doexp10(EFF_0(f) + EFF_0_A(f)*a1 + EFF_0_A2(f)*z8 + EFF_0_A3(f)*a2);
  c.alpha = ALPHA(f) + (a1*ALPHA_A(f) + z8*ALPHA_A2(f)) + a2*ALPHA_A3(f);
  c.delta = DELTA(f); // + (a1*DELTA_A(f) + z8*DELTA_A2(f))*expscale;
  c.beta = BETA(f) + a1*BETA_A(f) + z8*BETA_A2(f);
  c.gamma = GAMMA(f) + (a1*GAMMA_A(f) + z8*GAMMA_A2(f));
  if (c.gamma < 0) c.gamma = 0;
  c.lambda = doexp10(LAMBDA(f) + (a1*LAMBDA_A(f) + z8*LAMBDA_A2(f)));
  c.mu = MU(f) + a1*MU_A(f);
  c.kappa = KAPPA(f) + a1*KAPPA_A(f);
  c.passive_mass = 10.3 + z*0.5 - c.mu;
  c.scatter = SCATTER(f) + a1*SCATTER_A(f);
  //c.scatter = 0.212;
  c.scatter_corr = exp(pow(c.scatter*log(10), 2)/2.0);
  c.obs_scatter = (no_obs_scatter) ? 0 : SIGMA_CENTER + SIGMA_Z(f)*(z-0.1) + SIGMA_A(f)*a1;
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
  c.lm_slope = 1.0/c.alpha - 1.0;
  c.mpk = c.mu+c.kappa;
  c.combined_scatter = 0;
  c.valid = 1;

  c.qv = 2.232131 + 0.217770*a1 + 0.338483*z8;
  c.qsig = 0.365818 + 0.284595*a1 +  0.081106*z8;

  double comb_scatter = sqrt(c.obs_scatter*c.obs_scatter + c.scatter*c.scatter);
  double sm_at_m14 = 11.5; //calc_sm_at_m(14, c)+c.mu;
  double passive_frac = 1.0/(exp10fc(-1.3*(sm_at_m14-c.passive_mass))+1);
  double avg_offset_passive = (1.0-passive_frac)*c.kappa;
  double avg_offset_active = c.kappa - avg_offset_passive;
  c.combined_scatter = sqrt(pow(comb_scatter,2) + passive_frac*pow(avg_offset_passive,2) + (1.0-passive_frac)*pow(avg_offset_active,2)); 
  return c;
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

double calc_sfr_at_lv(double lv, struct smf c) {
  double vd = lv - c.v_1;
  double vd2 = vd/c.delta;
  double vd3 = (lv-c.qv)/c.qsig;
  double sfrac = 0.5+0.5*erf(-vd3*M_SQRT1_2);
  return (sfrac*c.epsilon * (1.0/(doexp10(c.alpha*vd) + doexp10(c.beta*vd)) + c.gamma*exp(-0.5*vd2*vd2)));
}

double calc_sfrac_at_lv(double lv, struct smf c) {
  double vd3 = (lv-c.qv)/c.qsig;
  double sfrac = 0.5+0.5*erf(-vd3*M_SQRT1_2);
  return (sfrac);
}


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
  double dlm_dlsm = gsl_spline_eval_deriv(steps[n].spline, sm, ga);
  double dndlogm = mf_cache(steps[n].scale, m);
  double z = 1 / steps[n].scale - 1;
  // if (z > 3.5 && z < 4.0)
	 //  printf("at z=%f, m=%f, sm=%f, sm_max=%f, dlm_dlsm=%f, dndlogm=%f.\n", z, m, sm, steps[n].smhm.sm_max, dlm_dlsm, dndlogm);
  if (!(dlm_dlsm>0)) { // && fit) {
    // printf("at z=%f, sm=%f, Falling SMHM.\n", (1 / steps[n].scale - 1), sm);
    // INVALIDATE(fit, "Falling SM(M) relation.");
    // return 0;
    dlm_dlsm = 1.0;
  }
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
  if (INVALID(*fit)) 
  {
	  // printf("invalid!\n");
	  return;
  }
  //if (steps[n].scale < 0.1) return;

  for (i=0; i<M_BINS; i++) {
    if (i && (steps[n].log_sm[i] <= steps[n].log_sm[i-1])) {
      if (steps[n].log_sm[i]==0) { steps[n].log_sm[i]=0; } // I totally don't understand why you would do this.
                                                                // In this case at the very massive end of the SMF,
                                                                // the corresponding halo mass would be very difficult
                                                                // to determine, and cause the crazy behavior there.
                                                                // For now I just put an empty block here.
      // if (steps[n].log_sm[i]==0) { steps[n].log_sm[i] = steps[n].log_sm[i-1]+0.001; }
      else {
	// sprintf(buffer, "Falling SMHM relation: SM(%f) = %f; SM(%f) = %f at scale %f (Mass at z=0: %f)\n", steps[n].med_hm_at_a[i], steps[n].log_sm[i], steps[n].med_hm_at_a[i-1], steps[n].log_sm[i-1], steps[n].scale, i*INV_BPDEX+M_MIN);	      
	//fprintf(stderr, "Falling SMHM relation: SM(%f) = %f; SM(%f) = %f at scale %f (Mass at z=0: %f)\n", steps[n].med_hm_at_a[i], steps[n].log_sm[i], steps[n].med_hm_at_a[i-1], steps[n].log_sm[i-1], steps[n].scale, i*INV_BPDEX+M_MIN);       
  INVALIDATE(fit, buffer);
	steps[n].smhm.valid = 0;
	return;
      }
    }

    if (steps[n].log_sm[i]>steps[n].smhm.sm_max) {
      sm_max = steps[n].smhm.sm_max = steps[n].log_sm[i];
      bin_max = i;
    }
    if (steps[n].log_sm[i]>-1000 && bin_min < 0 && i < M_BINS - 1 && steps[n].log_sm[i + 1] > 0) {
      bin_min = i;
      sm_min = steps[n].smhm.sm_min = steps[n].log_sm[i];
    }
  }

    if (bin_min < 0)
{
        sprintf(buffer, "All the stellar masses are zero.\n");
        //fprintf(stderr, "All the stellar masses are zero.\n");
        INVALIDATE(fit, buffer);
        return;
}
  
  // if (!((bin_max-bin_min+1) == M_BINS)) {
  //   // sprintf(buffer, "Incorrect number of bins to generate spline at scale %f\n", steps[n].scale);
  //   fprintf(stderr, "Incorrect number %d vs. %d of bins to generate spline at scale %f\n", (bin_max - bin_min + 1), (M_BINS), steps[n].scale);
  //   INVALIDATE(fit, buffer);
  //   return;
  // }

  // double z = 1 / steps[n].scale - 1;
  //   if (z > 3.5 && z < 4.0)
  //   {
  //     for (j=bin_min; j < bin_max; j++)
  //     {
  //       printf("z=%f, m=%f, log_sm=%f\n", z, steps[n].med_hm_at_a[j], steps[n].log_sm[j]);
  //     }
  //   }

  double *log_sm_tmp = NULL;
  double  *hm_tmp = NULL;
  double  *sfr_tmp = NULL;
  log_sm_tmp = malloc(sizeof(double)*(bin_max - bin_min + 1));
  hm_tmp = malloc(sizeof(double)*(bin_max - bin_min + 1));
  sfr_tmp = malloc(sizeof(double)*(bin_max - bin_min + 1));
  for (j=bin_min; j <= bin_max; j++)
  {
    // printf("log_sm=%f\n", steps[n].log_sm[j]);
    log_sm_tmp[j - bin_min] = steps[n].log_sm[j];
    hm_tmp[j - bin_min] = steps[n].med_hm_at_a[j];
    sfr_tmp[j - bin_min] = steps[n].sfr[j];
  }

  // printf("\n");

  // for (j=0; j < bin_max - bin_min + 1; j++)
  // {
  //   printf("log_sm=%f\n", log_sm_tmp[j]);
  // }

  steps[n].spline = gsl_spline_alloc(gsl_interp_cspline, bin_max - bin_min + 1);
  steps[n].spline_sfr = gsl_spline_alloc(gsl_interp_cspline, bin_max - bin_min + 1);

  gsl_spline_init(steps[n].spline, log_sm_tmp, hm_tmp, bin_max - bin_min + 1);
  gsl_spline_init(steps[n].spline_sfr, log_sm_tmp, sfr_tmp, bin_max - bin_min + 1);

  free(log_sm_tmp); free(hm_tmp); free(sfr_tmp);

  // gsl_spline_init(steps[n].spline, steps[n].log_sm+bin_min, steps[n].med_hm_at_a+bin_min, bin_max - bin_min + 1);
  // gsl_spline_init(steps[n].spline_sfr, steps[n].log_sm+bin_min, steps[n].sfr+bin_min, bin_max - bin_min + 1);

  i = M_BINS-1;
  for (j=SM_BINS-1; j>=0; j--) {
    sm = SM_MIN + (double)j*SM_INV_BPDEX;
    if (sm > sm_max || sm < sm_min) {
      steps[n].sfr_sm[j+SM_EXTRA] = 0;
      steps[n].smf[j+SM_EXTRA] = 0;
    }
    else {

      m = gsl_spline_eval(steps[n].spline, sm, &ga);
      //int64_t b = gsl_interp_accel_find(&ga, steps[n].log_sm+bin_min, (bin_max-bin_min+1), sm)+bin_min;
      steps[n].smf[j+SM_EXTRA] = calc_smf_at_m(m,sm,n,fit,&ga);
      steps[n].sfr_sm[j+SM_EXTRA] = gsl_spline_eval(steps[n].spline_sfr, sm, &ga)*steps[n].smf[j+SM_EXTRA];

      /*steps[n].sfr[b];
      if (b<bin_max) {
	double f = (sm-steps[n].log_sm[b])/(steps[n].log_sm[b+1]-steps[n].log_sm[b]);
	steps[n].sfr_sm[j+SM_EXTRA] +=  f*(steps[n].sfr[b+1]-steps[n].sfr[b]);
	}*/
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
  double scatter_corr = steps[n].smhm.scatter_corr; 
  double a200 = steps[n].scale/0.3782;
  double m200 = log10(1.115e12/0.68 / (pow(a200, -0.142441) + pow(a200, -1.78959)));
  double v200 = log10(200);
  
#pragma omp for schedule(dynamic,5) private(j,k,sm_hist_i,sm_hist_j,bin,t)
  for (i=0; i<M_BINS; i++) {
    double lv = v200 + (steps[n].med_hm_at_a[i]-m200)/3.0;
    steps[n].lv[i] = lv;
    steps[n].sfrac[i] = calc_sfrac_at_lv(lv, steps[n].smhm);
    steps[n].sfr[i] = calc_sfr_at_lv(lv, steps[n].smhm)*scatter_corr;
    //steps[n].log_sm[i] = calc_sm_at_m(steps[n].med_hm_at_a[i], steps[n].smhm);
    //steps[n].sm[i] = doexp10(steps[n].log_sm[i]);
    //steps[n].sm_avg[i] = steps[n].sm[i]*scatter_corr;
    steps[n].sm[i] = steps[n].sfr[i]*steps[n].dt;
    steps[n].sm_from_icl[i] = 0;
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
    }

    calc_old_sm(n,i);
    calc_new_sm_and_sfr(n,i,fit);
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



void calc_new_sm_and_sfr(int n, int i, struct smf_fit *fit) {
  double dt;
  int64_t j;
  char buffer[1024];
  dt = steps[n].dt;

  if (!steps[n].t[i]) {
    //steps[n].log_sm[i] = 0;
    //steps[n].new_sm[i] = 0;
    //steps[n].sm[i] = 0;
    return;
  }

  steps[n].new_sm[i] = steps[n].sfr[i] * dt * steps[n].smloss[n];

  //  steps[n].new_sm[i] = (steps[n].sm_avg[i] - steps[n].old_sm[i]);
  float dm = (i+0.5)*INV_BPDEX+M_MIN-steps[n].smhm.icl_m;
  float w = 16.0-steps[n].smhm.icl_m;
  //if (steps[n].smhm.beta>0) { w = 4.0 / (steps[n].smhm.beta * M_LN10); }
  if (w < 0.1) w = 0.1;
  steps[n].icl_frac[i] = 1.0 - 1.0/(1.0+exp(4.0*dm/w));
  if (steps[n].icl_frac[i] > 0.95) steps[n].icl_frac[i] = 0.95;
  steps[n].new_sm[i] *= 1.0/(1.0 - steps[n].icl_frac[i]);

  steps[n].sm_avg[i] = steps[n].old_sm[i] + steps[n].new_sm[i];
  steps[n].sm[i] = steps[n].sm_avg[i] / steps[n].smhm.scatter_corr;
  steps[n].log_sm[i] = (steps[n].sm[i] > 1) ? log10(steps[n].sm[i]) : 0;

  if (steps[n].sm_icl[i] < steps[n].icl_frac[i]*steps[n].new_sm[i]) {
    if (steps[n].new_sm[i] > 0 && steps[n].scale > 0.15 && i>BPDEX) {
      // sprintf(buffer, "ICL depleted (hm: %f; a: %f; sm: %e, new sm: %e, icl: %e)\n",(i+0.5)*INV_BPDEX+M_MIN, steps[n].scale,  steps[n].sm_avg[i], steps[n].new_sm[i], steps[n].icl_frac[i]);
      //fprintf(stderr, "ICL depleted (hm: %f; a: %f; sm: %e, new sm: %e, icl: %e)\n",(i+0.5)*INV_BPDEX+M_MIN, steps[n].scale,  steps[n].sm_avg[i], steps[n].new_sm[i], steps[n].icl_frac[i]);
      INVALIDATE(fit, buffer);
      steps[n].icl_frac[i] = steps[n].sm_icl[i]/steps[n].new_sm[i];
    }
    else steps[n].icl_frac[i] = 0;
  }
  steps[n].merged_frac[i] = 0;
  if (steps[n].icl_frac[i] && steps[n].sm_icl[i]) {
    double frac_from_icl = steps[n].new_sm[i]*steps[n].icl_frac[i] / steps[n].sm_icl[i];
    if (steps[n].sm_inc[i]>0) {
      steps[n].merged_frac[i] = steps[n].new_sm[i]*steps[n].icl_frac[i]/steps[n].sm_inc[i];
      if (steps[n].merged_frac[i]>1) steps[n].merged_frac[i] = 1;
    }
    if (frac_from_icl) {
      for (j=0; j<n; j++) {
	steps[n].sm_hist[i*num_outputs + j] += frac_from_icl*steps[n].icl_stars[i*num_outputs + j];
	steps[n].icl_stars[i*num_outputs + j] *= (1.0-frac_from_icl);
      }
    }
    steps[n].new_sm[i] *= (1.0-steps[n].icl_frac[i]);
  }


  steps[n].new_sm[i] /= steps[n].smloss[n];
  steps[n].sm_hist[i*num_outputs + n] = steps[n].new_sm[i];
  //steps[n].sfr[i] = steps[n].new_sm[i] / dt;

  if ((steps[n].t[i] > REAL_ND_CUTOFF) && (steps[n].sm_avg[i] > 0.17*exp10(steps[n].med_hm_at_a[i])) && (!no_z_scaling) && n>2) {
    /*sprintf(buffer, "SM exceeds baryon fraction (sm: %e, m: %e, scale: %f!\n",
	    steps[n].sm_avg[i], exp10(M_MIN+(i+0.5)*INV_BPDEX), steps[n].scale);
	    INVALIDATE(fit, buffer);*/
  }

  //if (steps[n].t[i]<REAL_ND_CUTOFF)
  //  steps[n].sfr[i] = 0;

  if ((!no_z_scaling) && (steps[n].sfr[i] < 0 || !isfinite(steps[n].sfr[i]))) {
    char buffer[1024];
    // sprintf(buffer, "Negative SFR at a = %f and m = %f (ND=%g)! (ICL m=%f; ICL frac=%f)\n", steps[n].scale, 
	   // i*INV_BPDEX + M_MIN, steps[n].t[i], steps[n].smhm.icl_m, steps[n].icl_frac[i]);
    //fprintf(stderr, "Negative SFR at a = %f and m = %f (ND=%g)! (ICL m=%f; ICL frac=%f)\n", steps[n].scale, 
     // i*INV_BPDEX + M_MIN, steps[n].t[i], steps[n].smhm.icl_m, steps[n].icl_frac[i]);
    INVALIDATE(fit, buffer);
  }
}

