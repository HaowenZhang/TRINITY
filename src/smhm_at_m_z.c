#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "all_smf.h"
#include "observations.h"
#include "calc_sfh.h"
#include "sm_limits.h"
#include "expcache2.h"

void read_params(char *buffer, double *params, int max_n);

int main(int argc, char **argv) {
  float z, m; //, mass_start, mass_end;
  FILE *input;
  char buffer[1024];
  struct smf_fit smf_fit;

  if (argc<4) {
    printf("Usage: %s m z red_smf_mcmc.dat\n", argv[0]);
    exit(1);
  }
  
  m = atof(argv[1]);
  z = atof(argv[2]);
  gen_exp10cache();

  if (!(input = fopen(argv[3], "r"))) {
    printf("Couldn't open file %s for reading!\n", argv[3]);
    exit(1);
  }

  fgets(buffer, 1024, input);
  read_params(buffer, smf_fit.params, NUM_PARAMS);

  struct smf c = smhm_at_z(z, smf_fit);
  printf("%f\n", calc_sm_at_m(m, c));
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
