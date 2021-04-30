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
  float *vals, best_val;
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
#define oparam(x) FILE *out_ ## x = fopen("systematics/" #x ".dat", "w");
#define cparam(x) fclose(out_ ## x);
#define paramrange(x) { for (j=0; j<num_entries; j++) vals[j] = fits_at_z[j].x;\
    qsort(vals, num_entries, sizeof(float), float_compare); \
    best_val = vals[(int)(0.5*num_entries)];				\
    fprintf(out_ ## x, "%e %e %e %e\n", z, best_val, vals[SIGMA_UP] - best_val,\
	    best_val - vals[SIGMA_DOWN]); }

  oparam(alpha);
  oparam(delta);
  oparam(sm_0);
  oparam(m_1);
  oparam(lm_slope);
  oparam(mpk);
  oparam(combined_scatter);
  oparam(gamma);
  oparam(scatter);
  oparam(obs_scatter);
  oparam(mu);
  oparam(kappa);
  oparam(sm_completeness);
  oparam(csfr_completeness);
  oparam(sfr_sm_corr);

  for (i=0; i<LOG_A_BINS+1; i++) {
    double a = A_START*pow(A_END/A_START, (double)i/LOG_A_BINS);
    double z = 1.0/a - 1.0;
    int64_t j;
    for (j=0; j<num_entries; j++) fits_at_z[j] = smhm_at_z(z, fits[j]);
    paramrange(alpha);
    paramrange(delta);
    paramrange(sm_0);
    paramrange(m_1);
    paramrange(gamma);
    paramrange(scatter);
    paramrange(obs_scatter);
    paramrange(mu);
    paramrange(kappa);
    paramrange(sm_completeness);
    paramrange(csfr_completeness);
    paramrange(sfr_sm_corr);
    paramrange(lm_slope);
    paramrange(mpk);
    paramrange(combined_scatter);
  }

  cparam(alpha);
  cparam(delta);
  cparam(sm_0);
  cparam(m_1);
  cparam(lm_slope);
  cparam(mpk);
  cparam(combined_scatter);
  cparam(gamma);
  cparam(scatter);
  cparam(obs_scatter);
  cparam(mu);
  cparam(kappa);
  cparam(sm_completeness);
  cparam(csfr_completeness);
  cparam(sfr_sm_corr);
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
