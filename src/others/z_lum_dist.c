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
#include "all_luminosities.h"
#include "mt_rand.h"
#include "universe_time.h"
#include <string.h>
#include <gsl/gsl_sf.h>

extern int64_t num_outputs;
extern struct timestep *steps;

#define NUM_L_THRESHES 100

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
  int64_t i, j, lt;
  enum lum_types ltype;
  float fm, ft, t;
  char buffer[1024];
  float *lums_r = NULL;
  float *lums_i = NULL;
  float *lums_m = NULL;

  if (argc<4+NUM_PARAMS) {
    fprintf(stderr, "Usage: %s mass_cache filter thresh (mcmc output)\n", argv[0]);
    exit(1);
  }
  for (i=0; i<NUM_PARAMS; i++)
    the_smf.params[i] = atof(argv[i+4]);

  for (ltype=0; ltype<LUM_TYPES; ltype++) {
    if (!strcasecmp(argv[2], lum_names[ltype])) break;
  }
  if (ltype == LUM_TYPES) {
    fprintf(stderr, "Unknown filter type %s!\n", argv[2]);
    fprintf(stderr, "Allowable types:");
    for (i=0; i<LUM_TYPES; i++) fprintf(stderr, " %s", lum_names[i]);
    fprintf(stderr, "\n");
    exit(1);
  }

  float mthresh = atof(argv[3]);

  r250_init(87L);
  gen_exp10cache();
  setup_psf(1);
  load_mf_cache(argv[1]);
  init_timesteps();
  INVALID(the_smf) = 0;
  calc_sfh(&the_smf);
  load_luminosities("../smf_mcmc2/data/mags_ab_zall_dust.dat");
  gen_all_luminosities(L_sloan_r, &lums_r, -1, 0);
  gen_all_luminosities(L_sloan_i, &lums_i, -1, 0);
  gen_all_luminosities(ltype, &lums_m, -1, 0);

  sprintf(buffer, "z_lums_r.dat");
  FILE *f_r = check_fopen(buffer, "w");
  sprintf(buffer, "z_lums_i.dat");
  FILE *f_i = check_fopen(buffer, "w");

  float *nd_r[NUM_L_THRESHES];
  float *nd_i[NUM_L_THRESHES];
  size_t arr_size = sizeof(float)*(num_outputs*3);
  for (lt = 0; lt<NUM_L_THRESHES; lt++) {
    nd_r[lt] = check_realloc(NULL, arr_size, "arr");
    memset(nd_r[lt], 0, arr_size);
    nd_i[lt] = check_realloc(NULL, arr_size, "arr");
    memset(nd_i[lt], 0, arr_size);
  }
  float *nd_m = check_realloc(NULL, arr_size, "arr");
  memset(nd_m, 0, arr_size);

  for (t=1; t<(num_outputs-1)*3; t++) {
    ft = (float)(((int)t)%3) / 3.0;
    i = t/3;
    float tscale = steps[i].scale + ft*(steps[i+1].scale - steps[i].scale);
    float a1 = tscale - 0.5*(steps[i+1].scale - steps[i].scale)/3.0;
    float a2 = tscale + 0.5*(steps[i+1].scale - steps[i].scale)/3.0;
    float dz = fabs(1.0/a1 - 1.0/a2);
    float vol = comoving_volume(1.0/a1 - 1.0) - comoving_volume(1.0/a2 - 1.0);
    float scale = 1.0/tscale;

    for (m=8; m<15.1; m+=0.05) {
      j = (m-M_MIN-INV_BPDEX*0.5)*BPDEX;
      fm = (m-M_MIN)*BPDEX - j;
      //float sfr = biterp(steps[i].sfr[j], steps[i+1].sfr[j],
      //                   steps[i].sfr[j+1], steps[i+1].sfr[j+1],
      //                   ft, fm);
      //float sm1 = calc_sm_at_m(m, steps[i].smhm)+steps[i].smhm.mu;
      //float sm2 = calc_sm_at_m(m, steps[i+1].smhm)+steps[i+1].smhm.mu;
      //float sm = sm1 + ft*(sm2-sm1);
      //if (!isfinite(sfr)) sfr = -1000;

      float nd = biterp(steps[i].t[j], steps[i+1].t[j],
                        steps[i].t[j+1], steps[i+1].t[j+1],
                        ft, fm);
      nd += log10(BPDEX);
      
      float lum_r = lum_at_hm_z(lums_r, m, scale-1.0);
      float lum_i = lum_at_hm_z(lums_i, m, scale-1.0);
      float scatter = (2.5*steps[i].smhm.scatter);
      float norm = 1.0/(scatter*sqrt(2.0*M_PI));
      norm *= pow(10, nd)*0.05*vol/dz;
      if (!isfinite(norm)) continue;

      for (lt=0; lt<NUM_L_THRESHES; lt++) {
	float thresh = 20.0+lt*0.1;
	float dl = (lum_r-thresh)/scatter;
	if (isfinite(dl))
	  nd_r[lt][(int)t] += exp(-0.5*dl*dl)*norm;
	dl = (lum_i-thresh)/scatter;
	if (isfinite(dl))
	  nd_i[lt][(int)t] += exp(-0.5*dl*dl)*norm;
      }
      float lum_m = lum_at_hm_z(lums_m, m, scale-1.0);
      float dl = (lum_m - mthresh)/scatter;
      if (isfinite(dl)) nd_m[(int)t] += exp(-0.5*dl*dl)*norm;
    }
    nd_m[0] += nd_m[(int)t]*dz;
    for (lt=0; lt<NUM_L_THRESHES; lt++) {
      nd_r[lt][0] += nd_r[lt][(int)t]*dz;
      nd_i[lt][0] += nd_i[lt][(int)t]*dz;
    }
  }

  fprintf(f_r, "#Z P(z|L) for L=20, L=20.1, L=20.2, ...\n");
  fprintf(f_i, "#Z P(z|L) for L=20, L=20.1, L=20.2, ...\n");
  printf("#Z P(z|%s=%f)\n", argv[2], mthresh);
  for (t=1; t<(num_outputs-1)*3; t++) {
    ft = (float)(((int)t)%3) / 3.0;
    i = t/3;
    float tscale = steps[i].scale + ft*(steps[i+1].scale - steps[i].scale);
    float z = 1.0/tscale - 1.0;
    printf("%f %e\n", z, nd_m[(int)t]/nd_m[0]);
    fprintf(f_r, "%f", z);
    for (lt=0; lt<NUM_L_THRESHES; lt++) fprintf(f_r, " %e", nd_r[lt][(int)t]/nd_r[lt][0]);
    fprintf(f_r, "\n");
    fprintf(f_i, "%f", z);
    for (lt=0; lt<NUM_L_THRESHES; lt++) fprintf(f_i, " %e", nd_i[lt][(int)t]/nd_i[lt][0]);
    fprintf(f_i, "\n");
  }
  fclose(f_r);
  fclose(f_i);
  return 0;
}
