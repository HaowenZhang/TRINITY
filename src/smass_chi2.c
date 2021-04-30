#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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

void read_params(char *buffer, double *params, int max_n);
int float_compare(const void *a, const void *b);


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
  float z, m, mass_start, mass_end, mass_step = 0.125;
  float mul;
  FILE *input;
  char buffer[1024];
  struct smf *fits;
  struct smf_fit smf_fit;
  float chi2_lim;
  float best_smmr, up_smmr, down_smmr, smmr;
  float *chi2_list, *chi2_list_sorted;
  float sm_min, sm_max;
  int i, num_entries;
  struct smf (*model_cb)(float, void *) = NULL;

  if (argc<5) {
    printf("Usage: %s model z num_entries red_smf_mcmc.dat\n", argv[0]);
    exit(1);
  }
  
  i = atol(argv[1]);
  if (i == 0) model_cb = single_smf_cb;
  else if (i == 1) model_cb = linear_smf_cb;

  z = atof(argv[2]);
  sm_min = sm_min_mass(z);
  sm_max = sm_max_mass(z);
  num_entries = atol(argv[3]);
  fits = (struct smf *)malloc(sizeof(struct smf)*num_entries);
  chi2_list = (float *)malloc(sizeof(float)*num_entries);
  chi2_list_sorted = (float *)malloc(sizeof(float)*num_entries);

  if (!(input = fopen(argv[4], "r"))) {
    printf("Couldn't open file %s for reading!\n", argv[3]);
    exit(1);
  }

  for (i=0; i<num_entries; i++) {
    if (!fgets(buffer, 1024, input)) break;
    read_params(buffer, smf_fit.params, NUM_PARAMS+1);
    fits[i] = model_cb(z, (void *)(&smf_fit));
    chi2_list[i] = CHI2(smf_fit);
    if (!i) { fits[i].sm_0 += MU(smf_fit); }
  }
  fclose(input);
  num_entries = i;

  memcpy(chi2_list_sorted, chi2_list, sizeof(float)*num_entries);
  qsort(chi2_list_sorted, num_entries, sizeof(float), float_compare);
  chi2_lim = chi2_list_sorted[(int)(num_entries*ONE_SIGMA)];

  gen_exp10cache();
  mul = fits[0].m_1 - 0.5;
  mass_start = sm_to_log_halo_mass(sm_min, fits[0], mul);
  mass_end = sm_to_log_halo_mass(sm_max, fits[0], mul);
  if (mass_end > MASS_MAX) mass_end = MASS_MAX;
  //if (mass_end > 15.25 - z*0.5) mass_end = 15.25 - z*0.5;

  for (m=mass_start; m<mass_end+mass_step; m+=mass_step) {
    up_smmr = down_smmr = best_smmr = find_sm(m, fits[0]) - m;
    for (i=1; i<num_entries; i++) {
      if (chi2_list[i] > chi2_lim) continue;
      smmr = find_sm(m, fits[i]) - m;
      if (smmr < down_smmr) down_smmr = smmr;
      if (smmr > up_smmr) up_smmr = smmr;
    }
    printf("%f %f %f %f\n", m, best_smmr, up_smmr - best_smmr,
	   best_smmr - down_smmr);
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
