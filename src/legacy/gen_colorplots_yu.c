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

int main(int argc, char **argv)
{
  float m;
  //float *zs = {0.1, 1, 2, 3, 4, 5, 6, 7, 8};
  struct smf_fit the_smf;
  int64_t i, j;
  float fm, ft, t, fb = 0.17;
  FILE *sfr_ma_f;
  char buffer[1024];
  gsl_interp_accel ga = {0};

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

  sprintf(buffer, "sfr_ma_as_fn_of_sm.dat");
  sfr_ma_f = check_fopen(buffer, "w");

  float sm;
  for (sm=6; sm<12; sm+=0.05) {
    for (t=1; t<(num_outputs-1)*3; t++) {
      ft = (float)(((int)t)%3) / 3.0;
      i = t/3;
      float tscale = steps[i].scale + ft*(steps[i+1].scale - steps[i].scale);
      float scale = 1.0/tscale;
      float max_m = m_at_a(15.77, tscale);

      if ((sm < steps[i].smhm.sm_min) ||
	  (sm > steps[i].smhm.sm_max) ||
	  (sm < steps[i+1].smhm.sm_min) ||
	  (sm > steps[i+1].smhm.sm_max)) {
	fprintf(sfr_ma_f, "%f %f %f\n", scale, sm, -10.0);
	continue;
      }

      m = gsl_spline_eval(steps[i].spline, sm, &ga);
      float m2 = gsl_spline_eval(steps[i+1].spline, sm, &ga);
      m += (m2-m)*ft;
      if (m>max_m || m<M_MIN) {
	fprintf(sfr_ma_f, "%f %f %f\n", scale, sm, -10.0);
	continue;
      }

      j = (m-M_MIN-INV_BPDEX*0.5)*BPDEX;
      fm = (m-M_MIN)*BPDEX - j;

      float ma_rate = biterp(mar_from_mbins(i,j),mar_from_mbins(i+1,j),
			     mar_from_mbins(i,j+1),mar_from_mbins(i+1,j+1),
			     ft, fm);
      float sfr = biterp(steps[i].sfr[j], steps[i+1].sfr[j],
			 steps[i].sfr[j+1], steps[i+1].sfr[j+1],
			 ft, fm);
      ma_rate += log10(fb);
      float sfr_ma = sfr - ma_rate;
      fprintf(sfr_ma_f, "%f %f %f\n", scale, sm, sfr_ma);
    }
    fprintf(sfr_ma_f, "\n");
  }

  fclose(sfr_ma_f);
  return 0;
}
