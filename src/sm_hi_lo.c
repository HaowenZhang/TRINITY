#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "all_smf.h"
#include "observations.h"
#include "calc_sfh.h"
#include "expcache2.h"
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


int main(int argc, char **argv) {
  float z, m, mass_start, mass_end, mass_step = 0.125;
  FILE *input;
  char buffer[1024];
  struct smf *fits;
  struct smf_fit smf_fit;
  float *smmrs, best_smmr;
  int i, num_entries;
  float sm_min, sm_max;

  if (argc<4) {
    printf("Usage: %s z num_entries red_smf_mcmc.dat\n", argv[0]);
    exit(1);
  }
  
  z = atof(argv[1]);
  sm_min = sm_min_mass(z);
  sm_max = sm_max_mass(z);
  num_entries = atol(argv[2]);
  fits = (struct smf *)malloc(sizeof(struct smf)*num_entries);
  smmrs = (float *)malloc(sizeof(float)*num_entries);
  gen_exp10cache();

  if (!(input = fopen(argv[3], "r"))) {
    printf("Couldn't open file %s for reading!\n", argv[3]);
    exit(1);
  }

  for (i=0; i<num_entries; i++) {
    if (!fgets(buffer, 1024, input)) break;
    read_params(buffer, smf_fit.params, NUM_PARAMS);
    fits[i] = smhm_at_z(z, smf_fit);
  }
  fclose(input);
  num_entries = i;

  fits[0].sm_0 += fits[0].mu;
  fits[0].mu = 0;
  mass_start = 10;
  mass_end = 15;
  if (mass_end > MASS_MAX) mass_end = MASS_MAX;
  //if (mass_end > 15.25 - z*0.5) mass_end = 15.25 - z*0.5;

  for (m=mass_start; m<mass_end+mass_step; m+=mass_step) {
    for (i=0; i<num_entries; i++) {
      smmrs[i] = calc_sm_at_m(m, fits[i]);
    }
    best_smmr = smmrs[0];
    qsort(smmrs, num_entries, sizeof(float), float_compare);
    if (best_smmr + m > sm_min && best_smmr + m < sm_max)
      printf("%f %f %f %f\n", m, best_smmr, smmrs[SIGMA_UP] - best_smmr,
	     best_smmr - smmrs[SIGMA_DOWN]);
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
