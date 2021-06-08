#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "all_smf.h"
#include "observations.h"
#include "calc_sfh.h"
#include "expcache.h"
#include "sm_limits.h"
#include "mah.h"
#include "check_syscalls.h"
#include "mt_rand.h"
#include "mlist.h"
#include <gsl/gsl_spline.h>

#define LOG_A_BINS 40
#define A_END 1.0
#define A_START (1.0/9.0)
#define ONE_SIGMA 0.682689492137
#define SIGMA_UP ((int)(num_entries*(0.5+ONE_SIGMA/2.0)))
#define SIGMA_DOWN ((int)(num_entries*(0.5-ONE_SIGMA/2.0)))
#define DEFAULT_H0 0.7

void read_params(char *buffer, double *params, int max_n);

#define NUM_HM_BINS 16
#define HM_PER_BIN 0.25
#define HM_START 11
float *sfh[NUM_HM_BINS];
float *mass_at_z[NUM_HM_BINS];
float *dmdt_at_z[NUM_HM_BINS];
gsl_spline *dsm_spline, *dm_spline;
double *merger_table[M_BINS];

double calc_rem_sm(float *the_sfh, int64_t k) {
  int64_t i;
  double sm = 0;
  for (i=0; i<=k; i++)
    sm += the_sfh[i]*steps[i].dt*steps[k].smloss[i];
  return sm;
}

void print_sfh(int64_t i) {
  int64_t j,k;
  int64_t invalid = 1;
  double x[3], y[3], sm_from_icl[NUM_HM_BINS];
  gsl_interp_accel sm_accel = {0}, m_accel = {0};
  x[0] = steps[0].scale;
  x[1] = 0.5;
  x[2] = steps[num_outputs-1].scale;

  printf("#Scale ");
  for (j=0; j<NUM_HM_BINS; j++) printf("Mh,Sm,dMh/dt,SFR(Mh=%.2f@z=0) ", HM_START + (double)j*HM_PER_BIN);
  printf("\n");

  while (invalid) {
    invalid = 0;
    memset(sm_from_icl, 0, sizeof(double)*NUM_HM_BINS);
    memset(y, 0, sizeof(double)*3);
    if (i>0) for (j=0; j<2; j++) y[j] = normal_random(0, 1);
    for (j=0; j<NUM_HM_BINS; j++) {
      float m = j*HM_PER_BIN + HM_START;
      if (m_at_a(m, 0.5)+y[1]*m_scatter_at_a(m, 0.5) > m) {
	invalid = 1;
	continue;
      }
    }
  }

  gsl_spline_init(dm_spline, x, y, 3);
  memset(&m_accel, 0, sizeof(gsl_interp_accel));

  if (i>0)  for (j=0; j<3; j++) y[j] = normal_random(0, 1);
  gsl_spline_init(dsm_spline, x, y, 3);
  memset(&sm_accel, 0, sizeof(gsl_interp_accel));
  
  for (j=0; j<num_outputs; j++) {
    for (k=0; k<NUM_HM_BINS; k++) {
      double m = k*HM_PER_BIN + HM_START;
      double m_at_z = m_at_a(m, steps[j].scale) + m_scatter_at_a(m, steps[j].scale)*gsl_spline_eval(dm_spline, steps[j].scale, &m_accel);
      mass_at_z[k][j] = pow(10, m_at_z);
      double dm = mass_at_z[k][j];
      if (j) dm -= mass_at_z[k][j-1];
      dmdt_at_z[k][j] = dm / steps[k].dt;
      double mbin_at_a = BPDEX*(m_at_z-M_MIN);
      if (mbin_at_a<0) mbin_at_a = 0;
      int64_t b = mbin_at_a;
      double f = mbin_at_a - b;
      sfh[k][j] = (1.0-f)*steps[j].sfr[b] + f*steps[j].sfr[b+1];
      sfh[k][j] /= exp(0.5*pow(steps[j].smhm.scatter*log(10),2));
      sfh[k][j] *= pow(10, steps[j].smhm.scatter*gsl_spline_eval(dsm_spline, steps[j].scale, &sm_accel));
    }
  }
  double merger_corr[NUM_HM_BINS];
  for (k=0; k<NUM_HM_BINS; k++) merger_corr[k] = 1;
  for (j=num_outputs-1; j>=0; j--) {
    for (k=0; k<NUM_HM_BINS; k++) {
      double m = k*HM_PER_BIN + HM_START;
      double m_at_z = m_at_a(m, steps[j].scale) + m_scatter_at_a(m, steps[j].scale)*gsl_spline_eval(dm_spline, steps[j].scale, &m_accel);
      double mbin_at_a = BPDEX*(m_at_z-M_MIN);
      if (mbin_at_a<0) mbin_at_a = 0;
      int64_t b = mbin_at_a;
      double f = mbin_at_a - b;
      merger_corr[k] = (1.0-f)*merger_table[b][j] + f*merger_table[b+1][j];
      sfh[k][j] /= merger_corr[k];
    }
  }
  
  for (j=0; j<num_outputs; j++) {
    printf("%f ", steps[j].scale);
    for (k=0; k<NUM_HM_BINS; k++) {
      printf("%e %e %e %e ", mass_at_z[k][j], calc_rem_sm(sfh[k], j), dmdt_at_z[k][j], sfh[k][j]);
    }
    printf("\n");
  }
}


int main(int argc, char **argv) {
  FILE *input;
  char buffer[1024];
  struct smf_fit smf_fit;
  int64_t i,j,k, num_entries;

  if (argc<3) {
    printf("Usage: %s num_entries mf_cache.dat red2_smf_mcmc.dat\n", argv[0]);
    exit(1);
  }
  
  num_entries = atol(argv[1]);
  r250_init(87L);
  setup_psf(1);
  load_mf_cache(argv[2]);
  init_timesteps();
  input = check_fopen(argv[3], "r");
  gen_exp10cache();
  for (i=0; i<NUM_HM_BINS; i++) {
    sfh[i] = (float *)check_realloc(NULL, sizeof(float)*num_outputs, "Af");
    mass_at_z[i] = (float *)check_realloc(NULL, sizeof(float)*num_outputs, "Af");
    dmdt_at_z[i] = (float *)check_realloc(NULL, sizeof(float)*num_outputs, "Af");
  }
  

  for (i=0; i<M_BINS; i++)
    merger_table[i] = (double *)check_realloc(NULL, sizeof(double)*num_outputs, "Mt");
  for (i=0; i<num_outputs; i++) {
    for (j=0; j<M_BINS; j++) {
      merger_table[j][i] = 0;
      double merged = 0, mmp=0;
      for (k=0; k<M_BINS; k++) merged+=steps[j].merged[j*M_BINS + k];
      for (k=0; k<M_BINS; k++) mmp+=steps[j].mmp[j*M_BINS + k];
      merger_table[j][i] = mmp/(merged+mmp);
    }
  }


  dsm_spline = gsl_spline_alloc(gsl_interp_cspline, 3);
  dm_spline = gsl_spline_alloc(gsl_interp_cspline, 3);
  for (i=0; i<num_entries; i++) {
    if (!fgets(buffer, 1024, input)) break;
    read_params(buffer, smf_fit.params, NUM_PARAMS+2);
    calc_sfh(&smf_fit);
    print_sfh(0); //i);
  }
  fclose(input);
  return 0;
}


void read_params(char *buffer, double *params, int max_n) {
  int num_entries = 0;
  char *cur_pos = buffer, *end_pos;
  float val = strtod(cur_pos, &end_pos);
  while (cur_pos != end_pos && num_entries < max_n) {
    params[num_entries] = val;
    num_entries++;
    cur_pos=end_pos;
    while (*cur_pos==' ' || *cur_pos=='\t' || *cur_pos=='\n') cur_pos++;
    val = strtod(cur_pos, &end_pos);
  }
}