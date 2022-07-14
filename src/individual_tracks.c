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

float *sfh[4];
gsl_spline *dsm_spline, *dm_spline;
double *merger_table[M_BINS];

void print_sfh(int64_t i) {
  int64_t j,k;
  int64_t invalid = 1;
  double x[3], y[3], sm_from_icl[4];
  gsl_interp_accel sm_accel = {0}, m_accel = {0};
  x[0] = steps[0].scale;
  x[1] = 0.5;
  x[2] = steps[num_outputs-1].scale;

  printf("#Scale ");
  for (j=11; j<15; j++) printf("%.1f ", (double)j);
  printf("\n");

  while (invalid) {
    invalid = 0;
    memset(sm_from_icl, 0, sizeof(double)*4);
    memset(y, 0, sizeof(double)*3);
    if (i>0) for (j=0; j<2; j++) y[j] = normal_random(0, 1);
    for (j=0; j<4; j++) 
      if (m_at_a(j+11, 0.5)+y[1]*m_scatter_at_a(j+11, 0.5) > (j+11)) {
	invalid = 1;
	continue;
      }
    gsl_spline_init(dm_spline, x, y, 3);
    memset(&m_accel, 0, sizeof(gsl_interp_accel));

    if (i>0)  for (j=0; j<3; j++) y[j] = normal_random(0, 1);
    gsl_spline_init(dsm_spline, x, y, 3);
    memset(&sm_accel, 0, sizeof(gsl_interp_accel));
 
    for (j=0; j<num_outputs; j++) {
      for (k=0; k<4; k++) {
	double m_at_z = m_at_a(k+11, steps[j].scale) + m_scatter_at_a(k+11, steps[j].scale)*gsl_spline_eval(dm_spline, steps[j].scale, &m_accel);
	double mbin_at_a = BPDEX*(m_at_z-M_MIN);
	if (mbin_at_a<0) mbin_at_a = 0;
	int64_t b = mbin_at_a;
	double f = mbin_at_a - b;
	sfh[k][j] = (1.0-f)*steps[j].sfr[b] + f*steps[j].sfr[b+1];
	sfh[k][j] /= exp(0.5*pow(steps[j].smhm.scatter*log(10),2));
	sfh[k][j] *= pow(10, steps[j].smhm.scatter*gsl_spline_eval(dsm_spline, steps[j].scale, &sm_accel));
      }
    }
    double merger_corr[4];
    for (k=0; k<4; k++) merger_corr[k] = 1;
    for (j=num_outputs-1; j>=0; j--) {
      for (k=0; k<4; k++) {
	double m_at_z = m_at_a(k+11, steps[j].scale) + m_scatter_at_a(k+11, steps[j].scale)*gsl_spline_eval(dm_spline, steps[j].scale, &m_accel);
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
      for (k=0; k<4; k++) {
	printf("%e ", sfh[k][j]);
      }
      printf("\n");
    }
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
  for (i=0; i<4; i++)
  sfh[i] = (float *)check_realloc(NULL, sizeof(float)*num_outputs, "Af");
  

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
    print_sfh(i);
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
