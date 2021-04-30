#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "observations.h"
#include "smf.h"
#include "all_smf.h"
#include "distance.h"
#include "integrate.h"
#include "mlist.h"
#include "calc_sfh.h"
#include "expcache2.h"
#include "check_syscalls.h"

extern int64_t num_outputs;
extern struct timestep *steps;

#define B_START 6
#define B_END 12
#define B_BPDEX 10
#define B_NB ((B_END-B_START)*B_BPDEX+1)
float conv_coeffs[B_NB][M_BINS];

void gen_conv_coeffs(struct smf_fit f) {
  int64_t i,j;
  struct smf smf = smhm_at_z(1.0/steps[num_outputs-1].scale-1.0, f);
  double scatter2 = smf.scatter*smf.scatter + smf.obs_scatter*smf.obs_scatter;
  for (i=0; i<M_BINS; i++) for (j=0; j<B_NB; j++) conv_coeffs[j][i] = 0;
  for (i=0; i<M_BINS; i++) {
    if (steps[num_outputs-1].sm[i] <= 1) continue;
    if (i*INV_BPDEX + M_MIN > 15) continue;
    double log_hsm = log10(steps[num_outputs-1].sm[i]);
    for (j=0; j<B_NB; j++) {
      double sm = B_START + j/((double)B_BPDEX);
      double delta_m = sm - log_hsm;
      double psf = exp(-0.5*(delta_m*delta_m/scatter2));
      //New Kappa:
      double passive_frac = 1.0/(exp10(-1.3*(sm-smf.passive_mass))+1.0);
      delta_m -= smf.kappa;
      psf = psf*passive_frac + 
	(1.0-passive_frac)*exp(-0.5*(delta_m*delta_m/scatter2));
      conv_coeffs[j][i] = psf;
    }
  }
  for (j=0; j<B_NB; j++) {
    double sum = 0;
    for (i=0; i<M_BINS; i++) sum += conv_coeffs[j][i];
    if (sum>0)
      for (i=0; i<M_BINS; i++) conv_coeffs[j][i] /= sum;
  }
}

float find_fraction_dt(int i, float f, float total_sm) {
  int64_t j=0, j2=0;
  double sum = 0, min_dt;
  for (j2=0; j2<num_outputs; j2++) {
    sum += steps[num_outputs-1].sm_hist[i*num_outputs+j2]*steps[num_outputs-1].smloss[j2];
    if (sum > f*(total_sm)) break;
  }
    min_dt = (j2-j)*steps[num_outputs-1].dt*(f*(total_sm)/sum);
  for (j=0; j<num_outputs; j++) {
    sum -= steps[num_outputs-1].sm_hist[i*num_outputs+j]*steps[num_outputs-1].smloss[j];
    sum -= steps[num_outputs-1].sm_hist[i*num_outputs+j2]*steps[num_outputs-1].smloss[j2];
    for (; j2<num_outputs; j2++) {
      sum += steps[num_outputs-1].sm_hist[i*num_outputs+j2]*steps[num_outputs-1].smloss[j2];
      if (sum > f*(total_sm)) break;
    }
    if (j2==num_outputs) break;
    double dt = (j2-j)*steps[num_outputs-1].dt*(f*(total_sm)/sum);
    if (dt < min_dt) {
      min_dt = dt;
      if (min_dt < 1.2e9) {
	printf("Mass: %f; Scale 1: %f; Scale 2: %f; Sum: %e; Total: %e\n", M_MIN + i*INV_BPDEX, steps[j].scale, steps[j2].scale, sum, total_sm);
      }
    }
  }
  return min_dt;
}

double convolve(float *array, int64_t j) {
  int64_t i;
  double sum=0;
  for (i=0; i<M_BINS; i++) sum += array[i]*conv_coeffs[j][i];
  return sum;
}

int main(int argc, char **argv)
{
  struct smf_fit the_smf;
  int i, j;
  float age;
  float ages[M_BINS], time50[M_BINS], time90[M_BINS];
  
  if (argc<2+NUM_PARAMS) {
    fprintf(stderr, "Usage: %s mass_cache (mcmc output)\n", argv[0]);
    exit(1);
  }
  for (i=0; i<NUM_PARAMS; i++)
    the_smf.params[i] = atof(argv[i+2]);

  setup_psf(1);
  load_mf_cache(argv[1]);
  init_timesteps();

  calc_sfh(&the_smf);
  gen_conv_coeffs(the_smf);
  for (i=0; i<M_BINS; i++) {
    float dt = 0;
    float total_sm = 0;
    if (i*INV_BPDEX + M_MIN > 15) continue;
    age = 0;
    if (steps[num_outputs-1].sm[i] <= 0) continue;
    for (j=num_outputs-1; j>=0; j--) {
      dt += steps[j].dt;
      age += steps[num_outputs-1].sm_hist[i*num_outputs+j]*steps[num_outputs-1].smloss[j]
	*(dt);
      total_sm += steps[num_outputs-1].sm_hist[i*num_outputs+j]*steps[num_outputs-1].smloss[j];
    }
    age /= total_sm;
    ages[i] = age;
    time50[i] = find_fraction_dt(i, 0.5, total_sm);
    time90[i] = find_fraction_dt(i, 0.9, total_sm);
    printf("%e %g %g %g %g %e\n", steps[num_outputs-1].sm[i], age, pow(10, M_MIN+(i+0.5)*INV_BPDEX), time90[i], time50[i], total_sm);
  }
  
  FILE *conv = check_fopen("plots/ages_conv.dat", "w");
  for (i=0; i<B_NB; i++) {
    double sm = pow(10, B_START + i/((double)B_BPDEX));
    fprintf(conv, "%e %g %g %g\n", sm, convolve(ages, i), convolve(time90, i), convolve(time50,i));
  }
  fclose(conv);
  return 0;
}
