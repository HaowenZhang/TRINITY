#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "all_smf.h"
#include "obs_smf.h"
#include "expcache.h"
#include "sm_limits.h"

#define SM_START 8
#define SM_END 12
#define SM_MAX 16
#define SM_MIN 0
#define MASS_MAX 15.2
#define ONE_SIGMA 0.682689492137
#define SIGMA_UP ((int)(num_entries*(0.5+ONE_SIGMA/2.0)))
#define SIGMA_DOWN ((int)(num_entries*(0.5-ONE_SIGMA/2.0)))
#define DEFAULT_H0 0.7

float *smmrs, *chi2;

void read_params(char *buffer, double *params, int max_n);
int float_compare(const void *a, const void *b);
int smmr_compare(const void *a, const void *b);
int chi2_compare(const void *a, const void *b);

float find_sm(float m, struct smf smf) {
#define hm(x) sm_to_log_halo_mass(x, smf, mul)
  float mul = smf.m_1 - 0.5; //smf.gamma*(float)(M_LN2/M_LN10);
  float sm = smf.sm_0;
  float trial_m = hm(sm);
  float step = 0.5;
  while (trial_m < m && sm<SM_MAX) { sm++; trial_m = hm(sm); }
  while (trial_m > m && sm>SM_MIN) { sm--; trial_m = hm(sm); }
  //Binary search time
  while (step > 1e-6) {
    if (trial_m < m) sm += step;
    else sm -= step;
    step *= 0.5;
    trial_m = hm(sm);
  }
  return sm;
#undef hm
}

int main(int argc, char **argv) {
  float z, z2, zmax, m, mass_start, mass_end, mass_step = 0.125;
  float mul;
  FILE *input;
  char buffer[1024];
  struct smf *fits, *fits2;
  struct smf_fit smf_fit;
  float best_smmr, chi2_limit;
  int *indices, best_index;
  float *counts, best_count, count;
  int i, num_entries, num_one_sigma;
  float sm_min, sm_max;
  struct smf (*model_cb)(float, void *) = NULL;

  if (argc<6) {
    printf("Usage: %s model z z2 num_entries red_smf_mcmc.dat\n", argv[0]);
    exit(1);
  }
  
  i = atol(argv[1]);
  if (i == 0) model_cb = single_smf_cb;
  else if (i == 1) model_cb = linear_smf_cb;

  z = atof(argv[2]);
  z2 = atof(argv[3]);
  zmax = (z>z2) ? z : z2;
  sm_min = sm_min_mass(zmax);
  sm_max = sm_max_mass(zmax);
  num_entries = atol(argv[4]);
  fits = (struct smf *)malloc(sizeof(struct smf)*num_entries);
  fits2 = (struct smf *)malloc(sizeof(struct smf)*num_entries);
  smmrs = (float *)malloc(sizeof(float)*num_entries);
  chi2 = (float *)malloc(sizeof(float)*num_entries);
  indices = (int *)malloc(sizeof(int)*num_entries);
  counts = (float *)malloc(sizeof(float)*num_entries);

  if (!(input = fopen(argv[5], "r"))) {
    printf("Couldn't open file %s for reading!\n", argv[5]);
    exit(1);
  }

  for (i=0; i<num_entries; i++) {
    if (!fgets(buffer, 1024, input)) break;
    read_params(buffer, smf_fit.params, NUM_PARAMS+1);
    fits[i] = model_cb(z, (void *)(&smf_fit));
    fits2[i] = model_cb(z2, (void *)(&smf_fit));
    chi2[i] = CHI2(smf_fit);
    indices[i] = i;
  }
  fclose(input);
  num_entries = i;
  num_one_sigma = num_entries*ONE_SIGMA;

  qsort(indices, num_entries, sizeof(int), chi2_compare);
  chi2_limit = chi2[indices[num_one_sigma]];

  gen_exp10cache();
  mul = fits[0].m_1 - 0.5;
  fits[0].sm_0 += fits[0].mu;
  fits[0].mu = 0;
  mass_start = sm_to_log_halo_mass(sm_min, fits[0], mul);
  mass_end = sm_to_log_halo_mass(sm_max, fits[0], mul);
  if (mass_end > MASS_MAX) mass_end = MASS_MAX;
  //if (mass_end > 15.25 - z*0.5) mass_end = 15.25 - z*0.5;

  for (m=mass_start; m<mass_end+mass_step; m+=mass_step) {
    for (i=0; i<num_entries; i++) {
      smmrs[i] = find_sm(m, fits[i]) - find_sm(m, fits2[i]);
    }
    best_smmr = smmrs[0];
    qsort(indices, num_entries, sizeof(int), smmr_compare);
    counts[0] = smmrs[indices[0]]; //(chi2[indices[0]] < chi2_limit) ? 1 : 0;
    for (i=1; i<num_entries; i++)
      counts[i] = counts[i-1] + smmrs[indices[i]];
	//((chi2[indices[i]] < chi2_limit) ? 1 : 0);

    best_index = 0;
    best_count = smmrs[indices[num_one_sigma]]-smmrs[indices[0]];
    for (i=1; i<(num_entries-num_one_sigma); i++) {
      count = smmrs[indices[num_one_sigma+i]]-smmrs[indices[i]];
      if (count < best_count) {
	best_count = count;
	best_index = i;
      }
    }

    printf("%f %f %f %f\n", m, best_smmr, smmrs[indices[best_index+num_one_sigma]] - best_smmr,
	   best_smmr - smmrs[indices[best_index]]);
  }
  return 0;
}

int float_compare(const void *a, const void *b) {
  float c = *((float *)a);
  float d = *((float *)b);
  if (c<d) return -1;
  if (c>d) return 1;
  return 0;
}

int smmr_compare(const void *a, const void *b) {
  float c = smmrs[*((int *)a)];
  float d = smmrs[*((int *)b)];
  if (c<d) return -1;
  if (c>d) return 1;
  return 0;
}

int chi2_compare(const void *a, const void *b) {
  float c = chi2[*((int *)a)];
  float d = chi2[*((int *)b)];
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
