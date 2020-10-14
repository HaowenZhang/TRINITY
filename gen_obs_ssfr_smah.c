#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_spline.h>
#include "observations.h"
#include "smf.h"
#include "all_smf.h"
#include "distance.h"
#include "integrate.h"
#include "mlist.h"
#include "calc_sfh.h"
#include "mah.h"
#include "check_syscalls.h"
#include "expcache2.h"
#include "sm_limits.h"

#define NUM_ZS 16
extern int64_t num_outputs;
extern struct timestep *steps;

double _total_sm(int64_t n, int64_t j) {
  int64_t i;
  double sm=0;
  for (i=0; i<=n; i++) sm += steps[n].sm_hist[j*num_outputs + i];
  return sm;
}

double nd_to_mass(double scale, double nd) {
  double mass = 17;
  double last_tot = 0;
  double tot = 0;
  for (; tot < nd; mass -= 0.01) {
    last_tot = tot;
    tot += pow(10, mf_cache(scale, mass))/100.0;
  }
  mass += 0.015;
  return (mass - 0.01*(nd-last_tot)/(tot-last_tot)) ;
}


float _mar_from_mbins(int64_t n, int64_t j) {
  int64_t i;
  if (!n) return pow(10, M_MIN+(j+0.5)*INV_BPDEX)/steps[n].dt;
  if (n>=num_outputs-1) return _mar_from_mbins(num_outputs-2, j);
  if (j>=M_BINS-1) return 0;
  //Find out average progenitor mass:
  double sum = 0;
  double count = 0;
  for (i=0; i<M_BINS; i++) {
    if (!steps[n].mmp[i*M_BINS + j]) continue;
    sum += pow(10, M_MIN+(i+0.5)*INV_BPDEX)*steps[n].mmp[i*M_BINS+j];
    count += steps[n].mmp[i*M_BINS+j];
  }
  if (!count) return 0;
  sum /= count;
  return ((pow(10, M_MIN+(j+0.5)*INV_BPDEX) - sum)/steps[n].dt);
}

float mar_from_mbins(int64_t n, int64_t j) {
  float mar1 = _mar_from_mbins(n,j);
  float mar2 = _mar_from_mbins(n+1,j);
  return (0.5*(mar1+mar2));
}

float biterp (float a, float b, float c, float d, float f1, float f2) {
  float al = log10(a);
  float bl = log10(b);
  float cl = log10(c);
  float dl = log10(d);
  float e = al+f1*(bl-al);
  float f = cl+f1*(dl-cl);
  return (e+f2*(f-e));
}

int main(int argc, char **argv)
{
  struct smf_fit the_smf;
  int64_t i, j;//, k, zi;
  float fm;
  double ft;
  FILE *direct_ssfr_f;
  char buffer[1024];
  double smah_ratio[M_BINS];
  double smah_norm[M_BINS];
  double smah_mass[M_BINS];
  double smah_started[M_BINS];
  double ssfrs[M_BINS];
  double sms[M_BINS];
  double mars[M_BINS];
  gsl_interp_accel acc = {0};
  gsl_spline *sp = gsl_spline_alloc(gsl_interp_linear, M_BINS);
  gsl_spline *sp2 = gsl_spline_alloc(gsl_interp_linear, M_BINS);

  if (argc<4+NUM_PARAMS) {
    fprintf(stderr, "Usage: %s mass_cache starting_z (mcmc output)\n", argv[0]);
    exit(1);
  }

  double starting_z = atof(argv[2]);
  for (i=0; i<NUM_PARAMS; i++)
    the_smf.params[i] = atof(argv[i+3]);

  gen_exp10cache();
  setup_psf(1);
  load_mf_cache(argv[1]);
  init_timesteps();
  INVALID(the_smf) = 0;
  calc_sfh(&the_smf);

  sprintf(buffer, "results/obs_ssfr_z%.1f.dat", starting_z);
  direct_ssfr_f = check_fopen(buffer, "w");
  fprintf(direct_ssfr_f, "#Z SM SSFR O_SSFR SHMAR\n");
  calc_step_at_z(starting_z, &i, &ft);
  for (j=0; j<M_BINS; j++) smah_started[j] = 0;

  for (; i>0; i--) {
    float z = 1.0/steps[i].scale - 1.0;
    if (z < starting_z) continue;
    float max_m = nd_to_mass(steps[i].scale, 1e-9);
    //float min_sm = sm_min_mass(1.0/steps[i].scale-1.0);
    //float max_sm = sm_max_mass(1.0/steps[i].scale-1.0);
    //fprintf(smhm_smah_f, "#Min: %f; Max: %f\n", min_sm, max_sm);
    double csfr = 0;
    double max_pred = 0;
    for (j=0; j<M_BINS; j++) {

      float mnow = m_now(1.0/(1.0+starting_z), M_MIN+j*INV_BPDEX);
      float mnow_u = m_now(1.0/(1.0+starting_z), M_MIN+(j+0.5)*INV_BPDEX);
      float mnow_d = m_now(1.0/(1.0+starting_z), M_MIN+(j-0.5)*INV_BPDEX);
      float m = m_at_a(mnow, steps[i].scale);
      float dm = m_at_a(mnow_u, steps[i].scale) - m_at_a(mnow_d, steps[i].scale);
      m = m_evolution_avg(M_MIN+j*INV_BPDEX, 1.0/(1.0+starting_z), steps[i].scale);
      dm = m_evolution_avg(M_MIN+(j+0.5)*INV_BPDEX, 1.0/(1.0+starting_z), steps[i].scale) - m_evolution_avg(M_MIN+(j-0.5)*INV_BPDEX, 1.0/(1.0+starting_z), steps[i].scale);


      dm *= BPDEX;
      fm = (m-M_MIN)*BPDEX-0.5;
      if (fm < 0) { sms[j] = j - M_BINS; ssfrs[j] = 0; mars[j] = 0; continue; }
      if (m > max_m) { sms[j] = j + M_BINS; ssfrs[j] = 0; mars[j] = 0; continue; }
      int64_t mb = fm;
      fm -= mb;
float mar1 = mar_from_mbins(i, mb);
      float mar2 = mar_from_mbins(i, mb+1);
      float mar = mar1 + fm*(mar2-mar1);
      float sm1 = _total_sm(i,mb); //steps[i].sm_avg[mb];
      float sm2 = _total_sm(i,mb+1); //steps[i].sm_avg[mb+1];
      float tsm = sm1 + fm*(sm2-sm1);
      sm1 = steps[i].sm_avg[mb];
      sm2 = steps[i].sm_avg[mb+1];
      float sm = sm1 + fm*(sm2-sm1);
      /*float icl1 = steps[i].sm_icl[mb];
      float icl2 = steps[i].sm_icl[mb+1];
      float icl = icl1 + fm*(icl2-icl1);
      if (icl > 0.2*sm) icl = 0.2*sm;*/
      float sfr1 = steps[i].sfr[mb];
      float sfr2 = steps[i].sfr[mb+1];
      float sfr = sfr1 + fm*(sfr2-sfr1);
      //sfr = calc_ssfr(log10(sm), 1.0/(steps[i].scale)-1.0)*sm;
      float n1 = steps[i].t[mb];
      float n2 = steps[i].t[mb+1];
      float n = (n1+(n2-n1)*fm)*dm;
      //sfr = calc_ssfr(log10(sm), 1.0/steps[i].scale-1.0)*sm;
      //float ssfr = 5.96e-11*pow(10, -0.35*(log10(sm)-11.03))*exp(-pow(10, log10(sm)-11.03));
      //sfr = ssfr*sm;
      //double scatter_corr = exp(pow(0.3*log(10), 2)/2.0);
    
      //mar = ma_rate_avg(mnow, steps[i].scale);//*scatter_corr;
      mar = ma_rate_avg_mnow(m, steps[i].scale);//*scatter_corr;

      float ratio = sfr/tsm/mar*pow(10, m);

      if (!smah_started[j]) { smah_ratio[j] = ratio; }
      double sm_scatter_corr = exp(pow(steps[i].smhm.scatter*log(10), 2)/2.0);

      if (!smah_started[j]) { smah_norm[j] = sm/sm_scatter_corr; }
      if (!smah_started[j]) { smah_mass[j] = m; }
      double pred = log10(smah_norm[j] * pow(10, smah_ratio[j]*(m-smah_mass[j])));
      double pred_sfr = smah_ratio[j]*pow(10, pred)*(mar / pow(10, m))*tsm/sm;
      double pred_ssfr = pred_sfr / pow(10, pred);
      if (pred > max_pred+0.0001) {
	max_pred = pred;
      } else {
	max_pred += 0.0001;
	pred = max_pred;
      }
      ssfrs[j] = pred_ssfr;
      sms[j] = pred;
      mars[j] = mar/pow(10,m);
      csfr += pred_sfr*n;
      smah_started[j] = 1;
    }
    gsl_spline_init(sp, sms, ssfrs, M_BINS);
    gsl_spline_init(sp2, sms, mars, M_BINS);
    fprintf(direct_ssfr_f, "%f", 1.0/steps[i].scale-1.0);
    for (j=16; j<22; j++) {
      float m = (double)j/2.0;
      float pred_ssfr=0, c_ssfr=0, pred_mar=0;
      if (m<max_pred) {
	pred_ssfr = gsl_spline_eval(sp, m, &acc);
	pred_mar = gsl_spline_eval(sp2, m, &acc);
	c_ssfr = calc_ssfr(m, 1.0/steps[i].scale-1.0);
      }
      fprintf(direct_ssfr_f, " %.2e %.2e %.2e", pred_ssfr, c_ssfr, pred_mar);
    }
    fprintf(direct_ssfr_f, "\n");
  }
  fclose(direct_ssfr_f);
  gsl_spline_free(sp);
  return 0;
}
