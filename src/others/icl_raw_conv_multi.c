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
#include "../../Rockstar/halo.h"
#include "inthash.h"
#include "luminosities.h"

#define INVALID_THRESH 0.03
#define mpeak min_pos_err
#define mvir mgrav
#define icl_part n_core
#define icl_level flags
#define hscale min_vel_err
#define mmp num_child_particles
#define invalid min_bulkvel_err
#define num_prev_halos num_p_halos
#define sm_pp child_r
#define lum_pp vmax_r
#define real_sm_pp rvmax

struct raw_icl_header {
  float scale, fract;
  int64_t snap, num_icl_halos;
  struct inthash ih;
};

struct halo *icl_halos = NULL;


void read_params(char *buffer, double *params, int max_n);
gsl_spline *dsm_spline;
struct smf_fit *smfs = NULL;

float biterp (float a, float b, float c, float d, float f1, float f2) {
  float al = log10(a);
  float bl = log10(b);
  float cl = log10(c);
  float dl = log10(d);
  float e = al+f1*(bl-al);
  float f = cl+f1*(dl-cl);
  return (e+f2*(f-e));
}


int main(int argc, char **argv) {
  FILE *input, *output;
  char buffer[1024];
  struct smf_fit smf_fit;
  int64_t i, num_entries;
  struct raw_icl_header rh;
  float *lums = NULL;

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
  load_luminosities("../smf_mcmc2/data/mags_vega_z0.dat");
  smfs = check_realloc(NULL, sizeof(struct smf_fit)*num_entries, "Af");
  
  for (i=0; i<num_entries; i++) {
    if (!fgets(buffer, 1024, input)) break;
    read_params(buffer, smfs[i].params, NUM_PARAMS+2);
  }
  num_entries = i;
  fclose(input);
  INVALID(smfs[0]) = 0;
  calc_sfh(smfs);
  gen_all_luminosities(L_cousins_i, &lums);

  input = check_fopen(argv[3], "r");
  fread(&rh, sizeof(struct raw_icl_header), 1, input);
  icl_halos = check_realloc(NULL, sizeof(struct halo)*rh.num_icl_halos, "Ah");
  fread(icl_halos, sizeof(struct halo), rh.num_icl_halos, input);
  fclose(input);
 
  int64_t n;
  double tf;
  calc_step_at_z(1.0/rh.scale-1.0, &n, &tf);
  float obs_scatter2_at_z = steps[n].smhm.obs_scatter*steps[n].smhm.obs_scatter;

  for (i=0; i<rh.num_icl_halos; i++) {
    float sscale = icl_halos[i].hscale;
    float smp = icl_halos[i].mpeak;
    int64_t num_icl = icl_halos[i].icl_part;    
    if (!num_icl) { icl_halos[i].sm_pp = 0; continue; }
    smp /= 0.7;
    smf_fit = smfs[0];
    struct smf c = smhm_at_z(1.0/sscale-1.0, smf_fit);
    float scatter = sqrt(obs_scatter2_at_z + c.scatter*c.scatter);
    float scatter_rand = normal_random(0, 1);
    float scatter_lum = scatter_rand * c.scatter;
    float scatter_sm = scatter_rand * scatter;
    float sm = pow(10, calc_sm_at_m(log10(smp), c) + c.mu + scatter_sm);
    sm /= (double)num_icl;
    icl_halos[i].sm_pp = sm;
    icl_halos[i].real_sm_pp = sm*pow(10, scatter_lum-scatter_sm);

    double tf, mf;
    int64_t n, m_bin;
    calc_step_at_z(1.0/sscale-1.0, &n, &tf);
    if (tf < 0) { n = num_outputs - 2; tf = 0.999; }
    mf = (log10(smp)-M_MIN)*BPDEX-0.5;
    if (mf < 0) mf=0;
    if (mf > M_BINS-2) mf=M_BINS-2;
    m_bin = mf;
    mf -= m_bin;
    float lum = biterp(lums[m_bin*num_outputs + n], lums[m_bin*num_outputs + n + 1],
		       lums[(m_bin+1)*num_outputs + n], lums[(m_bin+1)*num_outputs + n + 1],
		       tf, mf);
    icl_halos[i].lum_pp = pow(10, lum+scatter_lum) / (double)num_icl;
  }

  output = check_fopen(argv[3], "r+");
  fwrite(&rh, sizeof(struct raw_icl_header), 1, output);
  fwrite(icl_halos, sizeof(struct halo), rh.num_icl_halos, output);
  fclose(output);

  int64_t b=0;
  for (i=0; i<rh.num_icl_halos; i++)
    if (icl_halos[i].mpeak > icl_halos[b].mpeak) b=i;
  
  printf("Biggest ICL halo: M=%e; a=%f; Real SM=%e; Obs SM=%e; Num_p=%"PRId64"; Luminosity: %f\n", 
	 icl_halos[b].mpeak, icl_halos[b].hscale,
	 icl_halos[b].real_sm_pp*icl_halos[b].icl_part,
	 icl_halos[b].sm_pp*icl_halos[b].icl_part, icl_halos[b].icl_part,
	 log10(icl_halos[b].lum_pp*icl_halos[b].icl_part)*-2.5);
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
