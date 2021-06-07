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

float find_m_at_sm(float sm, struct smf_fit f) {
  struct smf c = smhm_at_z(0.01, f);
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
  float sm;

  if (argc<3) {
    printf("Usage: %s sm num_entries red_smf_mcmc.dat\n", argv[0]);
    exit(1);
  }
  
  sm = atof(argv[1]);
  num_entries = atol(argv[2]);
  fits = (struct smf_fit *)malloc(sizeof(struct smf_fit)*num_entries);
  fits_at_z = (struct smf *)malloc(sizeof(struct smf)*num_entries);
  masses = (float *)malloc(sizeof(float)*num_entries);
  vals = (float *)malloc(sizeof(float)*num_entries);
  gen_exp10cache();

  if (!(input = fopen(argv[3], "r"))) {
    printf("Couldn't open file %s for reading!\n", argv[3]);
    exit(1);
  }

  for (i=0; i<num_entries; i++) {
    if (!fgets(buffer, 1024, input)) break;
    read_params(buffer, smf_fit.params, NUM_PARAMS);
    fits[i] = smf_fit;
    masses[i] = find_m_at_sm(sm, fits[i]);
  }
  fclose(input);
  num_entries = i;

  EFF_0(fits[0]) += MU(fits[0]);
  MU(fits[0]) = 0;
  //  KAPPA(fits[0]) = 0;
  masses[0] = find_m_at_sm(sm, fits[0]);

  sprintf(buffer, "smtracks/sm_%.1f.dat", sm);
  FILE *output = fopen(buffer, "w");

  for (i=0; i<LOG_A_BINS+1; i++) {
    double a = A_START*pow(A_END/A_START, (double)i/LOG_A_BINS);
    double z = 1.0/a - 1.0;
    int64_t j;
    for (j=0; j<num_entries; j++) fits_at_z[j] = smhm_at_z(z, fits[j]);
    for (j=0; j<num_entries; j++) {
      vals[j] = calc_sm_at_m(m_at_a(masses[j],a), fits_at_z[j]);
      double passive_frac = 1.0/(exp10(-1.3*(vals[j]-fits_at_z[j].passive_mass))+1);
      vals[j] += fits_at_z[j].mu + (1.0-passive_frac)*fits_at_z[j].kappa;
    }
    qsort(vals, num_entries, sizeof(float), float_compare);
    best_val = vals[(int)(0.5*(num_entries))];
    fprintf(output, "%e %e %e %e\n", z, pow(10,best_val), pow(10, vals[SIGMA_UP]) - pow(10,best_val),
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
