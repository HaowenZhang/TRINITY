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

#define NUM_ZS 9
extern int64_t num_outputs;
extern struct timestep *steps;

float m_to_bin(float m) {
  m -= INV_BPDEX/2.0;
  if (m<M_MIN) return 0;
  if (m>(M_MAX-INV_BPDEX)) return (M_BINS-2);
  return ((m-M_MIN)*BPDEX);
}

float sm_at_m(float m, int i) {
  float f = m_to_bin(m);
  int b = f;
  f-=b;
  if (!b) return 0;
  float sm = steps[i].sm_avg[b] + f*(steps[i].sm_avg[b+1]-steps[i].sm_avg[b]);
  return sm;
}

float icl_at_m(float m, int i) {
  float f = m_to_bin(m);
  int b = f;
  f-=b;
  if (!b) return 0;
  float sm = steps[i].sm_icl[b] + f*(steps[i].sm_icl[b+1]-steps[i].sm_icl[b]);
  return sm;
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

#define NUM_MS 5

int main(int argc, char **argv)
{
  float m;
  double fb = 0.0455/(0.0455+0.226);
  struct smf_fit the_smf;
  int i, j;
  FILE *sm_f, *sm_f2, *sfrma, *smhm_hist, *icl_f, *sfrma2, *plain_sfh;
  char buffer[1024] = {0};

  if (argc<4+NUM_PARAMS) {
    fprintf(stderr, "Usage: %s mass_cache suffix (mcmc output)\n", argv[0]);
    exit(1);
  }
  for (i=0; i<NUM_PARAMS; i++)
    the_smf.params[i] = atof(argv[i+3]);

  setup_psf(1);
  load_mf_cache(argv[1]);
  init_timesteps();
  calc_sfh(&the_smf);
  
#define fopen_b(fh, m, b,...) {sprintf(b, __VA_ARGS__); fh = fopen(buffer, m);}
  fopen_b(sm_f, "w", buffer, "plots/sm_hist_%s.dat", argv[2]);
  fopen_b(sm_f2, "w", buffer, "plots/sm_hist_unnorm_%s.dat", argv[2]);
  fopen_b(icl_f, "w", buffer, "plots/icl_hist_%s.dat", argv[2]);
  fopen_b(sfrma, "w", buffer, "plots/sfr_ma_%s.dat", argv[2]);
  fopen_b(sfrma2, "w", buffer, "plots/sfr_ma_atz_%s.dat", argv[2]);
  fopen_b(smhm_hist, "w", buffer, "plots/smhm_hist_%s.dat", argv[2]);
  fopen_b(plain_sfh, "w", buffer, "plots/plain_sfh_%s.dat", argv[2]);

  for (j=11; j<16; j++) {
    m = j;
    float total_sm = 0;
    for (i=1; i<num_outputs-1; i++) {
      float m_then = m_at_a(m, steps[i].scale);
      float f = m_to_bin(m_then);
      int b = f;
      f-=b;
      if (!b) continue;
      float sm = pow(10, calc_sm_at_m(m_then, steps[i].smhm));
      float icl = icl_at_m(m_then, i);
      float sm2 = sm_at_m(m_then, i);
//steps[i].sm[b] + f*(steps[i].sm[b+1]-steps[i].sm[b]);
      float sfr = steps[i].sfr[b] + f*(steps[i].sfr[b+1]-steps[i].sfr[b]);
      float sfr2 = sm_at_m(m_at_a(m, steps[i+1].scale), i+1) -
	sm_at_m(m_at_a(m, steps[i-1].scale), i-1);
      sfr2 /= 0.7*(steps[i].dt + steps[i-1].dt);
      int b_now = (j-M_MIN)*BPDEX;
      total_sm += steps[num_outputs-1].smloss[i]*steps[num_outputs-1].sm_hist[b_now*num_outputs+i];
      float mar = ma_rate(m, steps[i].scale);
      fprintf(sm_f, "%f %f\n", 1.0/steps[i].scale, sm/pow(10,calc_sm_at_m(m, steps[num_outputs-1].smhm)));
      fprintf(sm_f2, "%f %f\n", 1.0/steps[i].scale, total_sm);
      if (!sm2) sm2 = 1;
      fprintf(icl_f, "%f %f\n", 1.0/steps[i].scale, icl/sm2);
      fprintf(sfrma, "%f %f\n", 1.0/steps[i].scale, sfr/mar/fb);
      fprintf(smhm_hist, "%f %f\n", 1.0/steps[i].scale, sm/pow(10, m_then));
      if (steps[i].sfr[(j-M_MIN)*BPDEX]) {
	fprintf(sfrma2, "%f %f\n", 1.0/steps[i].scale, steps[i].sfr[(j-M_MIN)*BPDEX]/mar_from_mbins(i,(j-M_MIN)*BPDEX)/fb);
      }
    }
    fprintf(icl_f, "\n");
    fprintf(sm_f, "\n");
    fprintf(sm_f2, "\n");
    fprintf(sfrma, "\n");
    fprintf(sfrma2, "\n");
    fprintf(smhm_hist, "\n");
  }

  for (j=0; j<M_BINS; j++) {
    for (i=1; i<num_outputs-1; i++) {
      fprintf(plain_sfh, "%f %e #%f\n", steps[i].scale, steps[num_outputs-1].sm_hist[num_outputs*j+i], j*INV_BPDEX+M_MIN);
    }
    fprintf(plain_sfh, "\n");
  }

  fclose(icl_f);
  fclose(sm_f);
  fclose(sm_f2);
  fclose(sfrma);
  fclose(sfrma2);
  fclose(smhm_hist);
  fclose(plain_sfh);
  return 0;
}
