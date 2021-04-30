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
gsl_spline *dsm_spline;
struct smf_fit *smfs = NULL;

int main(int argc, char **argv) {
  FILE *input;
  char buffer[1024];
  struct smf_fit smf_fit;
  int64_t i,j, num_entries;
  int64_t cur_z = 0;
#define NUM_ZS 23
  float zs[NUM_ZS] = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0};
  float last_scale = 0, last_m = 0;
  double x[3], y[3];
  gsl_interp_accel sm_accel = {0};

  if (argc<5) {
    printf("Usage: %s num_entries mf_cache.dat red2_smf_mcmc.dat mah.dat\n", argv[0]);
    exit(1);
  }
  
  num_entries = atol(argv[1]);
  r250_init(87L);
  setup_psf(1);
  load_mf_cache(argv[2]);
  init_timesteps();
  input = check_fopen(argv[3], "r");
  gen_exp10cache();
  smfs = check_realloc(NULL, sizeof(struct smf_fit)*num_entries, "Af");
  x[0] = 0.1;
  x[1] = 0.5;
  x[2] = 1.01;
  
  for (i=0; i<num_entries; i++) {
    if (!fgets(buffer, 1024, input)) break;
    read_params(buffer, smfs[i].params, NUM_PARAMS+2);
  }
  num_entries = i;
  fclose(input);

  dsm_spline = gsl_spline_alloc(gsl_interp_cspline, 3);
  input = check_fopen(argv[4], "r");
  printf("#Mass(z=0) Sub(1)/Cen(0) Scale_of_last_MM Vmax X Y Z");
  for (i=0; i<NUM_ZS; i++)
    printf(" M(z=%f) SM(z=%f)", zs[i], zs[i]);
  printf("\n");
  printf("#Masses in Msun/h; h=0.7\n");

  float cs=1;
  while (fgets(buffer, 1024, input)) {
    float scale, mass, tm;
    if (buffer[0]=='#') {
      buffer[strlen(buffer)-1] = 0;
      printf("%s", buffer+1);
      cur_z = 0;
      cs = 1.0/(1.0+zs[0]);
      last_scale = 1.0;
      last_m = atof(buffer+1);
      smf_fit = smfs[rand()%num_entries];
      for (j=0; j<3; j++) y[j] = normal_random(0, 1);
      gsl_spline_init(dsm_spline, x, y, 3);
      memset(&sm_accel, 0, sizeof(gsl_interp_accel));    
      continue;
    }
    if (sscanf(buffer, "%f %f", &scale, &mass)<2) {
      if (buffer[0]=='\n') printf("\n");
      continue;
    }
    if (scale < cs && last_scale > cs) {
      struct smf c = smhm_at_z(zs[cur_z], smf_fit);
      tm = mass + (cs-scale)*(last_m-mass)/(last_scale-scale);
      printf(" %e", tm);
      tm = log10(tm/0.7);
      float scatter = sqrt(c.obs_scatter*c.obs_scatter + c.scatter*c.scatter);
      float sm = pow(10, calc_sm_at_m(tm, c) + c.mu + scatter*gsl_spline_eval(dsm_spline, cs, &sm_accel));
      printf(" %e", sm*0.7);
      cur_z++;
      if (cur_z < NUM_ZS)
	cs = 1.0/(1.0+zs[cur_z]);
    }
    last_scale = scale;
    last_m = mass;
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
