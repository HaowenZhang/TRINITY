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
  FILE *input, *output, *output_w;
  char buffer[1024];
  struct smf_fit smf_fit;
  int64_t i, num_entries;

  if (argc<4) {
    printf("Usage: %s mf_cache.dat red2_smf_mcmc.dat icl_raw.dat\n", argv[0]);
    exit(1);
  }
  
  num_entries = 1;
  r250_init(87L);
  setup_psf(1);
  load_mf_cache(argv[1]);
  init_timesteps();
  input = check_fopen(argv[2], "r");
  gen_exp10cache();
  smfs = check_realloc(NULL, sizeof(struct smf_fit)*num_entries, "Af");
  
  for (i=0; i<num_entries; i++) {
    if (!fgets(buffer, 1024, input)) break;
    read_params(buffer, smfs[i].params, NUM_PARAMS+2);
  }
  num_entries = i;
  fclose(input);

  input = check_fopen(argv[3], "r");
  snprintf(buffer, 1024, "%s_sm", argv[3]);
  output = check_fopen(buffer, "w");
  snprintf(buffer, 1024, "%s_sm_weighted", argv[3]);
  output_w = check_fopen(buffer, "w");

  fprintf(output, "#ID X Y Z VX VY VZ HaloID HaloMass SourceMpeak SourceScale Num_ICL ICL_Level SM\n");
  fprintf(output, "#Masses in Msun\n");
  fprintf(output_w, "#ID X Y Z VX VY VZ HaloID HaloMass SourceMpeak SourceScale Num_ICL ICL_Level SM\n");
  fprintf(output_w, "#Masses in Msun\n");

  while (fgets(buffer, 1024, input)) {
    if (buffer[0]=='#') continue;
    int64_t id, hid, num_icl, icl_level;
    float x,y,z,vx,vy,vz,hm, smp, sscale, sm;
    if (sscanf(buffer, "%"SCNd64" %f %f %f %f %f %f %"SCNd64" %f %f %f %"SCNd64" %"SCNd64, &id, &x, &y, &z, &vx, &vy, &vz, &hid, &hm, &smp, &sscale, &num_icl, &icl_level) != 13) continue;
    if (!num_icl) continue;
    smp /= 0.7;
    hm /= 0.7;
    smf_fit = smfs[0];
    struct smf c = smhm_at_z(1.0/sscale-1.0, smf_fit);
    float scatter = sqrt(c.obs_scatter*c.obs_scatter + c.scatter*c.scatter);
    sm = pow(10, calc_sm_at_m(log10(smp), c) + c.mu + normal_random(0, 1)*scatter);
    sm /= (double)num_icl;
    fprintf(output, "%"PRId64" %.6f %.6f %.6f %.3f %.3f %.3f %"PRId64" %.3f %.4e %.5f %"PRId64" %"PRId64" %.4e\n",
	    id, x, y, z, vx, vy, vz, hid, hm, smp, sscale, num_icl,
	    icl_level, sm);
    if (sm>drand48()*1e7)
      fprintf(output_w, "%"PRId64" %.6f %.6f %.6f %.3f %.3f %.3f %"PRId64" %.3f %.4e %.5f %"PRId64" %"PRId64" %.4e\n",
	      id, x, y, z, vx, vy, vz, hid, hm, smp, sscale, num_icl,
	      icl_level, sm);
  }
  fclose(input);
  fclose(output);
  fclose(output_w);
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
