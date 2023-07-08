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

#define NUM_ZS 100
extern int64_t num_outputs;
extern struct timestep *steps;

static void read_params(char *buffer, double *data, int max_n) {
  int num_entries = 0;
  char *cur_pos = buffer, *end_pos;
  double val = strtod(cur_pos, &end_pos);
  while (cur_pos != end_pos && num_entries < max_n) {
    data[num_entries] = val;
    num_entries++;
    cur_pos=end_pos;
    while (*cur_pos==' ' || *cur_pos=='\t' || *cur_pos=='\n') cur_pos++;
    val = strtod(cur_pos, &end_pos);
  }
}

struct smf_fit parse_smf_fit(char *buffer) {
  struct smf_fit a;
  read_params(buffer, a.params, NUM_PARAMS);
  return a;
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
  float m;
  //float *zs = {0.1, 1, 2, 3, 4, 5, 6, 7, 8};
  struct smf_fit the_smf;
  int64_t i, j;
  float fm, ft, t, fb = 0.17;
  float max_t = log10(9);
  float *max_sfr = NULL;
  FILE *sfe_f, *input;
  char buffer[1024];
  float *sfe[NUM_ZS]={0};

  if (argc<4+NUM_PARAMS) {
    fprintf(stderr, "Usage: %s mass_cache num_runs mcmc_output\n", argv[0]);
    exit(1);
  }
  
  int64_t num_runs = atol(argv[2]);
  for (i=0; i<NUM_ZS; i++)
    sfe[i] = check_realloc(NULL, sizeof(float)*num_runs, "SFEs");
  input = check_fopen(argv[3]);
  gen_exp10cache();
  setup_psf(1);
  load_mf_cache(argv[1]);
  init_timesteps();
  
  while (fgets(buffer, 1024, input)) {
    the_smf = parse_smf_fit(buffer);

  INVALID(the_smf) = 0;
  calc_sfh(&the_smf);

  sprintf(buffer, "sfe_11.dat");
  sfe_f = check_fopen(buffer, "w");

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
      float ma_nd = ma_rate + nd;
     
      fprintf(sfr_ma_f, "%f %f %f\n", scale, m, sfr_ma);
      fprintf(ma_f, "%f %f %f\n", scale, m, ma_rate);
      fprintf(ma_nd_f, "%f %f %f\n", scale, m, ma_nd);
      fprintf(cond_sfr_f, "%f %f %f\n", scale, m, sfr-m_sfr);
      if (!(sfr >= -5)) sfr = -5;
      fprintf(sfr_f, "%f %f %f\n", scale, m, sfr);
    }
    fprintf(sfr_ma_f, "\n");
    fprintf(sfr_f, "\n");
    fprintf(ma_f, "\n");
    fprintf(cond_sfr_f, "\n");
  }

  fclose(sfr_f);
  fclose(sfr_ma_f);
  fclose(ma_f);
  fclose(cond_sfr_f);
  return 0;
}
