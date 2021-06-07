#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "observations.h"
#include "smf.h"
#include "all_smf.h"
#include "distance.h"
#include "integrate.h"
#include "mlist.h"
#include "calc_sfh.h"
#include "mah.h"
#include "check_syscalls.h"
#include "expcache2.h"
#include "stringparse.h"
#include "mt_rand.h"

extern int64_t num_outputs;
extern struct timestep *steps;

float biterp (float a, float b, float c, float d, float f1, float f2) {
  float al = log10(a);
  float bl = log10(b);
  float cl = log10(c);
  float dl = log10(d);
  float e = al+f1*(bl-al);
  float f = cl+f1*(dl-cl);
  return (e+f2*(f-e));
}


void add_sm_to_catalog(struct smf_fit f, char *filename) {
  char buffer[1024];
  FILE *input, *output;
  int n, ts;
  float sm, obs_sm;
#define NUM_CAT_INPUTS 22
  int64_t id, pid;
  float scale, m, vmax, mpeak, vpeak, rvir, x, y, z, vx, vy, vz, ra, dec, zcosmo, zlos, nd;
  int64_t rootid, dfid, lpid;
  SHORT_PARSETYPE;
  enum parsetype t[NUM_CAT_INPUTS];
  enum short_parsetype st[NUM_CAT_INPUTS] = 
    {F,D64,D64, F,F,F, F,F,F, F,F,F, F,F,F, F,F,F,F,D64,D64,D64};
  void *data[NUM_CAT_INPUTS] = {&scale, &id, &pid, &m, &vmax, &mpeak, &vpeak,
				&rvir, &x, &y, &z, &vx, &vy, &vz, &ra, &dec,
				&zcosmo, &zlos, &nd, &rootid, &dfid, &lpid};

  strncpy(buffer, filename, 1024);
  n = strlen(buffer);
  strncpy(&(buffer[n-3]), "sm_catalog.list", 1024-n);
  input = check_fopen(filename, "r");
  output = check_fopen(buffer, "w");
  fgets(buffer, 1024, input);  //Header line
  fprintf(output, "#Scale ID PID Mass Vmax Mass_peak Vmax_peak Rvir X Y Z VX VY VZ Ra Dec Z(cosmo) Z(los) Cumulative_ND(Vpeak) Tree_RootID Depthfirst_ID Last_Progenitor_ID Log(SM) Log(SM_obs) Log(SFR_obs)\n");

  ts = num_outputs-2;
  while (1) {
    fgets(buffer, 1024, input);  //Header line
    if (buffer[0] != '#') break;
    fprintf(output, "%s", buffer);
  }
  fprintf(output, "#Log(SM): Log10 of true stellar mass (in Msun)\n");
  fprintf(output, "#Log(SM_obs): Log10 of observed stellar mass (in Msun)\n");
  fprintf(output, "#Log(SFR_obs): Log10 of observed star formation rate (in Msun/yr)\n");

  for (n=0; n<NUM_CAT_INPUTS; n++) t[n] = st[n];

  do {
    n = stringparse(buffer, data, t, NUM_CAT_INPUTS);
    if (n < NUM_CAT_INPUTS) continue;
    while (steps[ts].scale > scale && ts>0) ts--;
    if (!ts) break;

    float tf = (scale - steps[ts].scale) / (steps[ts+1].scale - steps[ts].scale);
    if (!(mpeak>0) || !(m>0)) continue;
    float logm = log10(mpeak/0.7);
    float sm1 = calc_sm_at_m(logm, steps[ts].smhm);
    float sm2 = calc_sm_at_m(logm, steps[ts+1].smhm);
    sm = sm1 + tf*(sm2-sm1);
#define INTERP(y) (steps[ts].smhm.y + tf*(steps[ts+1].smhm.y - steps[ts].smhm.y))
    float ins_scatter = INTERP(scatter);
    float obs_scatter = INTERP(obs_scatter);
    float rho = INTERP(sfr_sm_corr);
    float mu = INTERP(mu);
    float kappa = INTERP(kappa);

    float sm_offset = normal_random(0, 1);
    float sfr_offset = normal_random(0, 1);
    sm = sm + sm_offset*ins_scatter;
    obs_sm = sm + normal_random(0, obs_scatter) + mu;
    
    float mf = (logm - M_MIN)*BPDEX;
    int64_t j = mf;
    mf -= j;
    float sfr = biterp(steps[ts].sfr[j], steps[ts+1].sfr[j],
			 steps[ts].sfr[j+1], steps[ts+1].sfr[j+1],
			 tf, mf);
    float scatter_corr = log10(exp(pow(ins_scatter*log(10), 2)/2.0));
    sfr -= scatter_corr;
    sfr += ins_scatter*(sm_offset*rho + sfr_offset*sqrt(1.0-rho*rho));
    if (pid >= 0) { sfr += log10(m/mpeak); }
    if (sfr - sm > 11) obs_sm += kappa;
    sfr += mu;

    fprintf(output, "%f %"PRId64" %"PRId64" %e %f %e %f %f %f %f %f %f %f %f %f %f %f %f %e %"PRId64" %"PRId64" %"PRId64" %f %f %f\n",
	    scale, id, pid, m, vmax, mpeak, vpeak, rvir, x, y, z, vx, vy, 
	    vz, ra, dec, zcosmo, zlos, nd, rootid, dfid, lpid, sm, obs_sm, sfr);
  } while (fgets(buffer, 1024, input));
  fclose(input);
  fclose(output);
}


int main(int argc, char **argv)
{
  int i;
  enum parsetype t[NUM_PARAMS];
  void *data[NUM_PARAMS];
  char buffer[1024] = {0};
  struct smf_fit base_smf;

  if (argc<3) {
    fprintf(stderr, "Usage: %s mass_cache catalog1 ... < mcmc_output\n", argv[0]);
    exit(1);
  }
  for (i=0; i<NUM_PARAMS; i++) {
    t[i] = PARSE_FLOAT64;
    data[i] = &(base_smf.params[i]);
  }

  fgets(buffer, 1024, stdin);
  if (stringparse(buffer, data, t, NUM_PARAMS) != NUM_PARAMS) {
    fprintf(stderr, "Could not parse input!\n");
    exit(1);
  }

  r250_init(87L);
  gen_exp10cache();
  setup_psf(1);
  load_mf_cache(argv[1]);
  init_timesteps();
  INVALID(base_smf) = 0;
  calc_sfh(&base_smf);

  for (i=2; i<argc; i++)
    add_sm_to_catalog(base_smf, argv[i]);
  return 0;
}
