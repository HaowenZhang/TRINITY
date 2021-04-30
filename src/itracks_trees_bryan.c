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
  double x[4], y[4];
  gsl_interp_accel sm_accel = {0};

  if (argc<5) {
    printf("Usage: %s num_entries mf_cache.dat red2_smf_mcmc.dat mah.dat z_interp\n", argv[0]);
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
  x[0] = 0.0;
  x[1] = 1.0/7.0;
  x[2] = 0.4;
  x[3] = 1.01;
 
  if (argc > 5) x[2] = 1.0/(1.0+atof(argv[5]));
 
  for (i=0; i<num_entries; i++) {
    if (!fgets(buffer, 1024, input)) break;
    read_params(buffer, smfs[i].params, NUM_PARAMS+2);
  }
  num_entries = i;
  fclose(input);

  dsm_spline = gsl_spline_alloc(gsl_interp_cspline, 4);
  input = check_fopen(argv[4], "r");
  printf("#Mass(z_final) Sub(1)/Cen(0) Scale_of_last_MM Vmax X Y Z\n");
  printf("#Scale Mass Mpeak Stellar Mass\n");
  printf("#Masses in Msun/h; h=0.7\n");

  while (fgets(buffer, 1024, input)) {
    float scale, mass, mpeak, tm;
    if (buffer[0]=='#') {
      //buffer[strlen(buffer)-1] = 0;
      //printf("%s", buffer+1);
      printf("%s", buffer);
      smf_fit = smfs[rand()%num_entries];
      for (j=0; j<4; j++) y[j] = normal_random(0, 1);
      gsl_spline_init(dsm_spline, x, y, 4);
      memset(&sm_accel, 0, sizeof(gsl_interp_accel));
      continue;
    }
    if (sscanf(buffer, "%f %f %f", &scale, &mass, &mpeak)<2) {
      if (buffer[0]=='\n') printf("\n");
      continue;
    }

    float z = 1.0/scale - 1.0;
    struct smf c = smhm_at_z(z, smf_fit);
    tm = mass;
    printf("%f %e %e ", scale, tm, mpeak);
    tm = log10(mpeak/0.7);
    float scatter = c.scatter;  //sqrt(c.obs_scatter*c.obs_scatter + c.scatter*c.scatter);
    float sm = pow(10, calc_sm_at_m(tm, c) + c.mu + scatter*gsl_spline_eval(dsm_spline, scale, &sm_accel));
    printf("%e\n", sm*0.7);
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
