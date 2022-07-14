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
#include "mah.h"
#include "check_syscalls.h"
#include "expcache2.h"
#include "sm_limits.h"

#define NUM_ZS 9
extern int64_t num_outputs;
extern struct timestep *steps;

double nd_to_mass(double scale, double nd) {
  double mass = 17;
  double last_tot = 0;
  double tot = 0;
  for (; tot < nd; mass -= 0.01) {
    last_tot = tot;
    tot += pow(10, mf_cache(scale, mass))/100.0;
  }
  mass += 0.015;
  return (mass - 0.01*(nd-last_tot)/(tot-last_tot)) ;
}


float _mar_from_mbins(int64_t n, int64_t j) {
  int64_t i;
  if (!n) return pow(10, M_MIN+(j+0.5)*INV_BPDEX)/steps[n].dt;
  if (n>=num_outputs-1) return _mar_from_mbins(num_outputs-2, j);
  if (j>=M_BINS-1) return 0;
  //Find out average progenitor mass:
  double sum = 0;
  double count = 0;
  for (i=0; i<M_BINS; i++) {
    if (!steps[n].mmp[i*M_BINS + j]) continue;
    sum += pow(10, M_MIN+(i+0.5)*INV_BPDEX)*steps[n].mmp[i*M_BINS+j];
    count += steps[n].mmp[i*M_BINS+j];
  }
  if (!count) return 0;
  sum /= count;
  return ((pow(10, M_MIN+(j+0.5)*INV_BPDEX) - sum)/steps[n].dt);
}

float mar_from_mbins(int64_t n, int64_t j) {
  float mar1 = _mar_from_mbins(n,j);
  float mar2 = _mar_from_mbins(n+1,j);
  return (0.5*(mar1+mar2));
}

float biterp (float a, float b, float c, float d, float f1, float f2) {
  float al = log10(a);
  float bl = log10(b);
  float cl = log10(c);
  float dl = log10(d);
  float e = al+f1*(bl-al);
  float f = cl+f1*(dl-cl);
  return (e+f2*(f-e));
}

float biterp2 (float a, float b, float c, float d, float f1, float f2) {
  float al = a;
  float bl = b;
  float cl = c;
  float dl = d;
  float e = al+f1*(bl-al);
  float f = cl+f1*(dl-cl);
  return log10(e+f2*(f-e));
}

float total_sm(int64_t i, int64_t j) {
  int64_t k;
  double t = 0;
  for (k=0; k<i; k++) {
    t += steps[i].sm_hist[num_outputs*j+k];
  }
  return t;
}

int main(int argc, char **argv)
{
  struct smf_fit the_smf;
  int64_t i, j;
  double ft;
  float z, m;
  char buffer[1024];

  if (argc<3+NUM_PARAMS) {
    fprintf(stderr, "Usage: %s mass_cache (mcmc output)\n", argv[0]);
    exit(1);
  }
  for (i=0; i<NUM_PARAMS; i++)
    the_smf.params[i] = atof(argv[i+2]);

  gen_exp10cache();
  setup_psf(1);
  load_mf_cache(argv[1]);
  init_timesteps();
  INVALID(the_smf) = 0;
  calc_sfh(&the_smf);

  snprintf(buffer, 1024, "results/ssfr_smah_colorplot.dat");
  FILE *output = check_fopen(buffer, "w");
  for (z=4; z<8; z+=0.02) {
    calc_step_at_z(z, &i, &ft);
    for (m=10; m<13.5; m+=0.02) {
      float mnow = m_evolution_avg(m, 1.0/5.0, 1.0/(1.0+z));
      float fm = (mnow - M_MIN)*BPDEX-0.5;
      j = fm;
      fm -= j;
      float ssfr = biterp2(steps[i].sfr[j]/total_sm(i,j), steps[i+1].sfr[j]/total_sm(i+1,j),
			  steps[i].sfr[j+1]/total_sm(i,j+1), steps[i+1].sfr[j+1]/total_sm(i+1, j+1),
                         ft, fm);
      float sm = biterp(total_sm(i,j), total_sm(i+1, j),
			total_sm(i, j+1), total_sm(i+1, j+1),
			ft, fm);
      float mar = ma_rate_avg_mnow(mnow, 1.0/(1.0+z));
      float smar = log10(mar) - mnow;
      fprintf(output, "%f %f %f %f %f %f %f\n", z, m, ssfr, smar, ssfr-smar, sm, mnow);
    }
  }
  fclose(output);
  return 0;
}
