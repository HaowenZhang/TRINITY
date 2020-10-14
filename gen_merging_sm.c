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

#define NUM_ZS 9
extern int64_t num_outputs;
extern struct timestep *steps;


double peri_cdf[1001];
void gen_peri_cdf(double z, double mhost) {
  int64_t i;
  if (z>7) z=7;
  double g1 = pow(1.0+z, -4);
  double mstar = 12.42 - 1.56*z + 0.038*z*z + log10(0.7);
  double ratio = pow(10, mhost-mstar);
  double R1 = 0.450*(1.0-0.395*pow(g1*ratio, 0.109));
  if (R1 < 0.05) R1 = 0.05;
  for (i=0; i<1001; i++) {
    double r = (double)i/(1000.0);
    peri_cdf[i] = exp(-pow(r/R1, 0.85))/2.0;
    if (i>0) {
      peri_cdf[i] += exp(-pow((double)(i-1)/(1000.0*R1), 0.85))/2.0 + peri_cdf[i-1];
    }
  }
  for (i=0; i<1001; i++)
    peri_cdf[i] /= peri_cdf[1000];
}
double get_peri_cdf(double r, double scale, int64_t j) {
  static double last_scale = -1;
  static int64_t last_j = 0;
  if (j!=last_j || scale != last_scale) {
    gen_peri_cdf(1.0/scale - 1.0, M_MIN + (j+0.5)*INV_BPDEX);
  }
  if (r<0) return 0;
  if (r>1) return 1;
  return peri_cdf[(int)(r*1000)];
}


int main(int argc, char **argv)
{
  struct smf_fit the_smf;
  int64_t i, j, k;

  if (argc<4+NUM_PARAMS) {
    fprintf(stderr, "Usage: %s mass_cache extension (mcmc output)\n", argv[0]);
    exit(1);
  }
  for (i=0; i<NUM_PARAMS; i++)
    the_smf.params[i] = atof(argv[i+3]);

  gen_exp10cache();
  setup_psf(1);
  load_mf_cache(argv[1]);
  init_timesteps();
  INVALID(the_smf) = 0;
  calc_sfh(&the_smf);

  for (i=1; i<num_outputs; i++) {
    double msm = 0;
    double msm_major = 0;
    double msm_10 = 0;
    double msm_majhalo = 0;
    double msm_maj_01 = 0;
    double sfr = 0;
    for (j=0; j<M_BINS; j++) {
      msm += steps[i].sm_inc[j]*steps[i].t[j];
      sfr += steps[i].new_sm[j]*steps[i].t[j];
      if (steps[i].t[j] < 1e-9) continue;
      for (k=0; k<M_BINS; k++) {
	double lgsm = log10(steps[i].sm[k]/1e10)/(sqrt(2)*steps[i].smhm.scatter);
	double frac_ab_10 = 0.5+0.5*erf(lgsm/(sqrt(2)*steps[i].smhm.scatter));
	msm_10 += frac_ab_10*steps[i].merged[k*M_BINS + j]*steps[i-1].sm_avg[k];
	double lgsm2 = log10(steps[i].sm_avg[k]/(0.3*steps[i].sm_avg[j]))/(sqrt(2)*steps[i].smhm.scatter);
	double frac_major = 0.5+0.5*erf(lgsm2);
	msm_major += frac_major*steps[i].merged[k*M_BINS + j]*steps[i-1].sm_avg[k];
	if (k>j-3) {
	  msm_majhalo += steps[i].merged[k*M_BINS + j]*steps[i-1].sm_avg[k];
	  msm_maj_01 += steps[i].merged[k*M_BINS + j] * steps[i-1].sm_avg[k] * get_peri_cdf(0.1, steps[i].scale, j);
	}
      }
    }
    msm_10 /= steps[i].dt;
    msm_major /= steps[i].dt;
    msm /= steps[i].dt;
    msm_majhalo /= steps[i].dt;
    msm_maj_01 /= steps[i].dt;
    sfr /= steps[i].dt;
    printf("%f %e %e %e %e %e %e\n", 1.0/steps[i].scale, msm, msm_major, msm_10, msm_majhalo, msm_maj_01, sfr);
  }
  return 0;
}
