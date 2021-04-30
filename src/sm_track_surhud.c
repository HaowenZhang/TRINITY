#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "all_smf.h"
#include "observations.h"
#include "calc_sfh.h"
#include "expcache.h"
#include "sm_limits.h"
#include "mah.h"

#define LOG_A_BINS 40
#define A_END 1.0
#define A_START (1.0/9.0)
#define ONE_SIGMA 0.682689492137
#define SIGMA_UP ((int)(num_entries*(0.5+ONE_SIGMA/2.0)))
#define SIGMA_DOWN ((int)(num_entries*(0.5-ONE_SIGMA/2.0)))
#define DEFAULT_H0 0.7

void read_params(char *buffer, double *params, int max_n);
int float_compare(const void *a, const void *b);

float find_m_at_sm(float sm, struct smf_fit f, float z) {
  struct smf c = smhm_at_z(z, f);
  double m = c.m_1;
  double sm_trial = calc_sm_at_m(m, c);
  while (fabs(sm_trial-sm) > 0.001) {
    double sm_trial2 = calc_sm_at_m(m+0.01, c);
    m = m + (sm-sm_trial)*(0.01/(sm_trial2-sm_trial));
    sm_trial = calc_sm_at_m(m, c);
  }
  return m;
}


int main(int argc, char **argv) {
  FILE *input;
  char buffer[1024];
  struct smf_fit *fits;
  struct smf_fit smf_fit;
  struct smf *fits_at_z;
  float *masses;
  float *vals, best_val;
  int64_t i, num_entries;
  float m;

  if (argc<3) {
    printf("Usage: %s m z num_entries red_smf_mcmc.dat\n", argv[0]);
    exit(1);
  }
  
  m = atof(argv[1]);
  float z0 = atof(argv[2]);
  num_entries = atol(argv[3]);
  fits = (struct smf_fit *)malloc(sizeof(struct smf_fit)*num_entries);
  fits_at_z = (struct smf *)malloc(sizeof(struct smf)*num_entries);
  masses = (float *)malloc(sizeof(float)*num_entries);
  vals = (float *)malloc(sizeof(float)*num_entries);
  gen_exp10cache();

  if (!(input = fopen(argv[4], "r"))) {
    printf("Couldn't open file %s for reading!\n", argv[3]);
    exit(1);
  }

  for (i=0; i<num_entries; i++) {
    if (!fgets(buffer, 1024, input)) break;
    read_params(buffer, smf_fit.params, NUM_PARAMS);
    fits[i] = smf_fit;
    //masses[i] = find_m_at_sm(sm, smf_fit);
  }
  fclose(input);
  num_entries = i;

  //  EFF_0(fits[0]) += MU(fits[0]);
  //  MU(fits[0]) = 0;

  for (i=0; i<num_entries; i++) {
    fits_at_z[i] = smhm_at_z(z0, fits[i]);
    masses[i] = calc_sm_at_m(m, fits_at_z[i]);
  }

  float m0 = m_now(1.0/(1.0+z0), m);

  sprintf(buffer, "smtracks/sm_hist_rel_%.1f_z%.2f.dat", m, z0);
  FILE *output = fopen(buffer, "w");
  fprintf(output, "#a SM/SM(z=%f) Err+ Err-\n", z0);

  for (i=0; i<LOG_A_BINS+1; i++) {
    double a = A_START*pow(A_END/A_START, (double)i/LOG_A_BINS);
    double z = 1.0/a - 1.0;
    if (z < z0) continue;
    int64_t j;
    for (j=0; j<num_entries; j++) fits_at_z[j] = smhm_at_z(z, fits[j]);
    for (j=0; j<num_entries; j++) vals[j] = calc_sm_at_m(m_at_a(m0,a), fits_at_z[j]) - masses[j];
    best_val = vals[0]; //(int)(0.5*num_entries)];
    qsort(vals, num_entries, sizeof(float), float_compare);
    fprintf(output, "%e %e %e %e\n", a, pow(10,best_val), pow(10, vals[SIGMA_UP]) - pow(10,best_val),
	    pow(10,best_val) - pow(10,vals[SIGMA_DOWN]));
  }

  fclose(output);
  return 0;
}

int float_compare(const void *a, const void *b) {
  float c = *((float *)a);
  float d = *((float *)b);
  if (c<d) return -1;
  if (c>d) return 1;
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
