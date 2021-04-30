#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
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

#define NUM_ZS 9
extern int64_t num_outputs;
extern struct timestep *steps;
double smar_threshold[M_BINS] = {0};
double weights[M_BINS] = {0};
#define MAR_SCATTER 0.3

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
  float m;
  //float *zs = {0.1, 1, 2, 3, 4, 5, 6, 7, 8};
  struct smf_fit the_smf;
  int64_t i, j;
  float fm, ft, t, fb = 0.17;
  float *max_sfr = NULL;
  FILE *mar_thresh;
  float *passive = NULL;
  char buffer[1024];

  if (argc<3+NUM_PARAMS) {
    fprintf(stderr, "Usage: %s mass_cache (mcmc output)\n", argv[0]);
    exit(1);
  }
  for (i=0; i<NUM_PARAMS; i++)
    the_smf.params[i] = atof(argv[i+2]);

  gen_exp10cache();
  setup_psf(1);
  load_mf_cache(argv[1]);
  init_timesteps();
  INVALID(the_smf) = 0;
  calc_sfh(&the_smf);
  passive = check_realloc(passive, sizeof(float)*M_BINS*num_outputs, "");
  memset(passive, 0, sizeof(float)*M_BINS*num_outputs);

  max_sfr = malloc(sizeof(float)*num_outputs);
  for (i=0; i<num_outputs; i++) {
    max_sfr[i] = 0;
    float max_m = m_at_a(15.77, steps[i].scale);
    for (j=0; j<M_BINS; j++) {
      float m = M_MIN+j*INV_BPDEX;
      if (m>max_m) break;
      if (steps[i].sfr[j] > max_sfr[i]) max_sfr[i] = steps[i].sfr[j];
    }
  }

  for (m=8; m<M_MAX; m+=0.05) {
    j = (m-M_MIN-INV_BPDEX*0.5)*BPDEX;
    fm = (m-M_MIN)*BPDEX - j;
    float avg_sfr_ma = 0;
    float total_sf = 0;
    for (t=1; t<(num_outputs-1)*3; t++) {
      ft = (float)(((int)t)%3) / 3.0;
      i = t/3;
      float tscale = steps[i].scale + ft*(steps[i+1].scale - steps[i].scale);
      //float scale = 1.0/tscale;
      float max_m = m_at_a(15.77, tscale);
      /*float m_sfr = biterp(max_sfr[i], max_sfr[i+1],
			   max_sfr[i], max_sfr[i+1],
			   ft, fm);*/
      float ma_rate = biterp(mar_from_mbins(i,j),mar_from_mbins(i+1,j),
			     mar_from_mbins(i,j+1),mar_from_mbins(i+1,j+1),
			     ft, fm);
      float sfr = biterp(steps[i].sfr[j], steps[i+1].sfr[j],
			 steps[i].sfr[j+1], steps[i+1].sfr[j+1],
			 ft, fm);
      if (m>max_m) sfr -= 12.0*(m-max_m);
      float nd = biterp(steps[i].t[j], steps[i+1].t[j],
                        steps[i].t[j+1], steps[i+1].t[j+1],
                        ft, fm);
      nd += log10(BPDEX);
      ma_rate += log10(fb);
      float sfr_ma = sfr - ma_rate;
      if (m<max_m && sfr_ma > -10) {
	total_sf += 1;
	avg_sfr_ma += pow(10, sfr_ma);
      }
      //float ma_nd = ma_rate + nd;
      //float sfr_hm = sfr - m;
      float sm1 = calc_sm_at_m(m, steps[i].smhm);
      float sm2 = calc_sm_at_m(m, steps[i+1].smhm);
      float sm = sm1 + (ft*(sm2-sm1));
      double sm_scatter_corr = exp(pow(steps[i].smhm.scatter*log(10), 2)/2.0);
      float sfr_sm = sfr - sm - log10(sm_scatter_corr);
      //float ma_m = ma_rate - m;
      float pf = 0.5*(1.0+erf((-10.5 - sfr_sm)/sqrt(2.0*MAR_SCATTER*MAR_SCATTER)));
      sfr_ma -= log10(sm_scatter_corr);
      passive[i+j*num_outputs] = pf;
      if (m<max_m && isfinite(sfr_ma)) {
	//ssfr = -11 = SFR - SM = MAR + SFR_MA - SM
	//MAR - M = -11 + SM - SFR_MA - M
	smar_threshold[j] += -11.0 + sm - sfr_ma - m - log10(fb);
	weights[j]++;
      }
    }
  }

  mar_thresh = check_fopen("results/mar_thresh.dat", "w");
  fprintf(mar_thresh, "#M SMAR_THRESH\n");
  for (j=0; j<M_BINS; j++) {
    if (!weights[j]) continue;
    if (!isfinite(smar_threshold[j])) continue;
    fprintf(mar_thresh, "%e %e\n", pow(10, M_MIN+j*INV_BPDEX), pow(10, smar_threshold[j]/weights[j]));
  }
  fclose(mar_thresh);

  for (i=0; i<num_outputs; i += 10) {
    snprintf(buffer, 1024, "results/passive_frac_%f.dat", steps[i].scale);
    FILE *passive_f = check_fopen(buffer, "w");
    fprintf(passive_f, "#SM Passive Passive_Exp\n");
    double pm = steps[i].smhm.passive_mass;
    double sm;
    for (sm=7; sm<12; sm += 0.1) {
      double w = 0;
      double pf = 0;
      for (j=0; j<M_BINS; j++) {
	double c = steps[i].t[j]*exp(-0.5*pow((sm-calc_sm_at_m(M_MIN+(j+0.5)*INV_BPDEX, steps[i].smhm))/steps[i].smhm.combined_scatter, 2));
	if (!isfinite(passive[i + j*num_outputs])) continue;
	pf += passive[i + j*num_outputs]*c;
	w += c;
      }
      if (!w) continue;
      fprintf(passive_f, "%f %f %f\n", sm, pf/w, 1.0/(pow(10, -1.3*(sm-pm))+1.0));
    }
    fclose(passive_f);
  }
  return 0;
}
