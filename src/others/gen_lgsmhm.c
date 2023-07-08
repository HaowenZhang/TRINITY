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


float find_m_at_sm(float sm, struct smf c) {
  double m = c.m_1;
  double sm_trial = calc_sm_at_m(m, c);
  while (fabs(sm_trial-sm) > 0.001) {
    double sm_trial2 = calc_sm_at_m(m+0.01, c);
    m = m + (sm-sm_trial)*(0.01/(sm_trial2-sm_trial));
    sm_trial = calc_sm_at_m(m, c);
  }
  return m;
}

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
  int64_t i;

  if (argc<4+NUM_PARAMS) {
    fprintf(stderr, "Usage: %s mass_cache stellar_mass (mcmc output)\n", argv[0]);
    exit(1);
  }

  for (i=0; i<NUM_PARAMS; i++)
    the_smf.params[i] = atof(argv[i+3]);

  double starting_sm = atof(argv[2]);
  gen_exp10cache();
  setup_psf(1);
  load_mf_cache(argv[1]);
  init_timesteps();
  INVALID(the_smf) = 0;
  calc_sfh(&the_smf);

  for (i=0; i<num_outputs; i++) {
    float m = find_m_at_sm(starting_sm, steps[i].smhm);
    float mar = ma_rate_avg_mnow(m, steps[i].scale);
    float sm2 = calc_sm_at_m(m+0.01,steps[i].smhm);
    float sm1 = calc_sm_at_m(m-0.01,steps[i].smhm);
    float mb = (m-M_MIN)*BPDEX-0.5;
    int64_t mb1 = mb;
    float mf = mb-mb1;
    float sfr = steps[i].sfr[mb1] + mf*(steps[i].sfr[mb1+1]-steps[i].sfr[mb1]);
    float tsm1 = _total_sm(i, mb1);
    float tsm2 = _total_sm(i, mb1+1);
    float tsm = tsm1 + mf*(tsm2-tsm1);
    float ssfr = sfr / tsm;
    if (!isfinite(ssfr)) ssfr = 0;
    float sm_av =  steps[i].sm_avg[mb1] + mf*(steps[i].sm_avg[mb1+1]-steps[i].sm_avg[mb1]);
    float icl = steps[i].sm_icl[mb1] + mf*(steps[i].sm_icl[mb1+1]-steps[i].sm_icl[mb1]);
    float smar = mar/pow(10,m);
    float obs_ssfr = calc_ssfr(starting_sm, 1.0/steps[i].scale-1.0);
    printf("%f %f %f %e %f %f %e\n", 1.0/steps[i].scale-1.0, 
	   (sm2-sm1)/0.02, 
	   ssfr/smar, 
	   smar, 
	   tsm/sm_av, 
	   obs_ssfr / smar * sm_av/tsm/pow(10, steps[i].smhm.kappa), 
	   smar*tsm/sm_av);
  }
  return 0;
}
