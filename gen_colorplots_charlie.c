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
#include "mah.h"
#include "check_syscalls.h"
#include "expcache2.h"

#define NUM_ZS 9
extern int64_t num_outputs;
extern struct timestep *steps;

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
  FILE *sfr_ma_f, *ma_f, *cond_sfr_f, *ma_nd_f, *sfr_f;
  char buffer[1024];

  if (argc<4+NUM_PARAMS) {
    fprintf(stderr, "Usage: %s mass_cache extension (mcmc output)\n", argv[0]);
    exit(1);
  }
  for (i=0; i<NUM_PARAMS; i++)
    the_smf.params[i] = atof(argv[i+3]);

  gen_exp10cache();
  setup_psf(1);
  load_mf_cache(argv[1]);
  init_timesteps();
  INVALID(the_smf) = 0;
  calc_sfh(&the_smf);

  sprintf(buffer, "sfr_ma.dat");
  sfr_ma_f = check_fopen(buffer, "w");
  sprintf(buffer, "sfr_charlie.dat");
  sfr_f = check_fopen(buffer, "w");
  sprintf(buffer, "ma_nd.dat");
  ma_nd_f = check_fopen(buffer, "w");
  sprintf(buffer, "ma.dat");
  ma_f = check_fopen(buffer, "w");
  sprintf(buffer, "cond_sfr.dat");
  cond_sfr_f = check_fopen(buffer, "w");

  FILE *sfr_sm_f = check_fopen("sfr_over_sm.dat", "w");
  FILE *sm_weighted_sfr_ma_f = check_fopen("sm_weighted_sfr_ma.dat", "w");
  FILE *sfr_hm_f = check_fopen("sfr_over_hm.dat", "w");
  FILE *ma_m_f = check_fopen("ma_over_m.dat", "w");

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

  for (m=9; m<15.1; m+=0.05) {
    j = (m-M_MIN-INV_BPDEX*0.5)*BPDEX;
    fm = (m-M_MIN)*BPDEX - j;
    float avg_sfr_ma = 0;
    float total_sf = 0;
    for (t=1; t<(num_outputs-1)*3; t++) {
      //for (t=0; t<max_t; t+=0.007) {
      //float tscale = pow(10, -t);
      //for (i=num_outputs-2; i>=0; i--)
      //	if (steps[i].scale < tscale) break;
      //ft = (t - log10(1.0/steps[i].scale))/log10(steps[i].scale/steps[i+1].scale);
      ft = (float)(((int)t)%3) / 3.0;
      i = t/3;
      float tscale = steps[i].scale + ft*(steps[i+1].scale - steps[i].scale);
      float scale = 1.0/tscale;
      float max_m = m_at_a(15.77, tscale);
      float m_sfr = biterp(max_sfr[i], max_sfr[i+1],
			   max_sfr[i], max_sfr[i+1],
			   ft, fm);
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
      float ma_nd = ma_rate + nd;
      float sfr_hm = sfr - m;
      float sm1 = calc_sm_at_m(m, steps[i].smhm);
      float sm2 = calc_sm_at_m(m, steps[i+1].smhm);
      float sm = sm1 + (ft*(sm2-sm1));
      float sfr_sm = sfr - sm;
      float ma_m = ma_rate - m;
     
      fprintf(sfr_ma_f, "%f %f %f\n", scale, m, sfr_ma);
      fprintf(ma_f, "%f %f %f\n", scale, m, ma_rate);
      fprintf(ma_nd_f, "%f %f %f\n", scale, m, ma_nd);
      fprintf(cond_sfr_f, "%f %f %f\n", scale, m, sfr-m_sfr);
      fprintf(sfr_hm_f, "%f %f %f\n", scale, m, sfr_hm);
      fprintf(sfr_sm_f, "%f %f %f\n", scale, m, sfr_sm);
      fprintf(ma_m_f, "%f %f %f\n", scale, m, ma_m);
      if (!(sfr >= -5)) sfr = -5;
      fprintf(sfr_f, "%f %f %f\n", scale, m, sfr);
    }
    fprintf(sfr_ma_f, "\n");
    fprintf(sfr_f, "\n");
    fprintf(ma_f, "\n");
    fprintf(cond_sfr_f, "\n");
    if (total_sf)
      fprintf(sm_weighted_sfr_ma_f, "%f %f\n", m, log10(avg_sfr_ma/total_sf));
  }

  fclose(sfr_f);
  fclose(sfr_ma_f);
  fclose(ma_f);
  fclose(cond_sfr_f);
  fclose(sfr_sm_f);
  fclose(sfr_hm_f);
  fclose(ma_m_f);
  fclose(sm_weighted_sfr_ma_f);
  return 0;
}
