#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "all_smf.h"
#include "observations.h"
#include "calc_sfh.h"
#include "expcache.h"
#include "sm_limits.h"

#define LOG_A_BINS 40
#define A_END 1.0
#define A_START (1.0/9.0)
#define ONE_SIGMA 0.682689492137
#define SIGMA_UP ((int)(num_entries*(0.5+ONE_SIGMA/2.0)))
#define SIGMA_DOWN ((int)(num_entries*(0.5-ONE_SIGMA/2.0)))
#define DEFAULT_H0 0.7

void read_params(char *buffer, double *params, int max_n);
int float_compare(const void *a, const void *b);


int main(int argc, char **argv) {
  FILE *input;
  char buffer[1024];
  struct smf_fit *fits;
  struct smf_fit smf_fit;
  struct smf *fits_at_z;
  float *vals, best_val, z, m;
  int64_t i, num_entries;

  if (argc<3) {
    printf("Usage: %s num_entries red_smf_mcmc.dat\n", argv[0]);
    exit(1);
  }
  
  num_entries = atol(argv[1]);
  fits = (struct smf_fit *)malloc(sizeof(struct smf_fit)*num_entries);
  fits_at_z = (struct smf *)malloc(sizeof(struct smf)*num_entries);
  vals = (float *)malloc(sizeof(float)*num_entries);
  gen_exp10cache();

  if (!(input = fopen(argv[2], "r"))) {
    printf("Couldn't open file %s for reading!\n", argv[3]);
    exit(1);
  }

  for (i=0; i<num_entries; i++) {
    if (!fgets(buffer, 1024, input)) break;
    read_params(buffer, smf_fit.params, NUM_PARAMS);
    fits[i] = smf_fit;
  }
  fclose(input);
  num_entries = i;

  EFF_0(fits[0]) += MU(fits[0]);
  MU(fits[0]) = 0;

  for (z=0; z<10; z+=0.5) {
    if (z==0) z=0.1;
    for (i=0; i<num_entries; i++) fits_at_z[i] = smhm_at_z(z, fits[i]);
    sprintf(buffer, "systematics/sfr_frac_z%.1f.dat", z);
    FILE *out = fopen(buffer, "w");
    for (m=9; m<16; m+=0.05) {
      for (i=0; i<num_entries; i++) {
	float w = 16.0-fits_at_z[i].icl_m;
	float dm = m - fits_at_z[i].icl_m;
	if (w < 0.1) w = 0.1;
	vals[i] = 1.0/(1.0+exp(4.0*dm/w));
      }
      qsort(vals, num_entries, sizeof(float), float_compare);
      best_val = vals[(int)(0.5*num_entries)];
      fprintf(out, "%f %e %e %e\n", m, best_val, vals[SIGMA_UP] - best_val,
	      best_val - vals[SIGMA_DOWN]);
    }
    fclose(out);
    if (z==0.1) z=0;
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
