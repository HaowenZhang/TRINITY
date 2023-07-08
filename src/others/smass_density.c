#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "all_smf.h"
#include "observations.h"
#include "calc_sfh.h"
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

int main(int argc, char **argv) {
 //  float z, m, mass_start, mass_end, mass_step = 0.125;
 //  FILE *input;
 //  char buffer[1024];
 //  struct smf *fits;
 //  struct smf_fit smf_fit;
 //  float best_smmr, chi2_limit;
 //  int *indices, best_index;
 //  float *counts, best_count, count;
 //  int i, num_entries, num_one_sigma;
 //  float sm_min, sm_max;

 //  if (argc<4) {
 //    printf("Usage: %s z num_entries red_smf_mcmc.dat\n", argv[0]);
 //    exit(1);
 //  }
  
 //  z = atof(argv[1]);
 //  sm_min = sm_min_mass(z);
 //  sm_max = sm_max_mass(z);
 //  num_entries = atol(argv[2]);
 //  fits = (struct smf *)malloc(sizeof(struct smf)*num_entries);
 //  smmrs = (float *)malloc(sizeof(float)*num_entries);
 //  chi2 = (float *)malloc(sizeof(float)*num_entries);
 //  indices = (int *)malloc(sizeof(int)*num_entries);
 //  counts = (float *)malloc(sizeof(float)*num_entries);
 //  gen_exp10cache();

 //  if (!(input = fopen(argv[3], "r"))) {
 //    printf("Couldn't open file %s for reading!\n", argv[3]);
 //    exit(1);
 //  }

 //  for (i=0; i<num_entries; i++) {
 //    if (!fgets(buffer, 1024, input)) break;
 //    read_params(buffer, smf_fit.params, NUM_PARAMS+1);
 //    fits[i] = smhm_at_z(z, smf_fit);
 //    chi2[i] = CHI2(smf_fit);
 //    if (!i) { fits[i].sm_0 += MU(smf_fit); }
 //    indices[i] = i;
 //  }
 //  fclose(input);
 //  num_entries = i;
 //  num_one_sigma = num_entries*ONE_SIGMA;

 //  qsort(indices, num_entries, sizeof(int), chi2_compare);
 //  chi2_limit = chi2[indices[num_one_sigma]];

 //  mass_start = 10.0;
 //  mass_end = 15.0;

 //  for (m=mass_start; m<mass_end+mass_step; m+=mass_step) {
 //    for (i=0; i<num_entries; i++) {
 //      smmrs[i] = calc_sm_at_m(m, fits[i]) - m;
 //    }
 //    best_smmr = smmrs[0];
 //    qsort(indices, num_entries, sizeof(int), smmr_compare);
 //    counts[0] = smmrs[indices[0]]; //(chi2[indices[0]] < chi2_limit) ? 1 : 0;
 //    for (i=1; i<num_entries; i++)
 //      counts[i] = counts[i-1] + smmrs[indices[i]];
	// //((chi2[indices[i]] < chi2_limit) ? 1 : 0);

 //    best_index = 0;
 //    best_count = smmrs[indices[num_one_sigma]]-smmrs[indices[0]];
 //    for (i=1; i<(num_entries-num_one_sigma); i++) {
 //      count = smmrs[indices[num_one_sigma+i]]-smmrs[indices[i]];
 //      if (count < best_count) {
	// best_count = count;
	// best_index = i;
 //      }
 //    }

 //    if ((m + best_smmr > sm_min) && (m + best_smmr < sm_max)) {
 //      printf("%f %f %f %f\n", m, best_smmr, smmrs[indices[best_index+num_one_sigma]] - best_smmr,
	//      best_smmr - smmrs[indices[best_index]]);
 //    }
 //  }
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
