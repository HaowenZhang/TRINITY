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
  float fm, ft, t;
  float max_t = 9;
  FILE *sfr_f, *sfh_f;
  char buffer[1024];

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
  calc_sfh(&the_smf);

  argv[2] = "sandy";
  sprintf(buffer, "sfr_%s.dat", argv[2]);
  sfr_f = check_fopen(buffer, "w");
  sprintf(buffer, "sfh_%s.dat", argv[2]);
  sfh_f = check_fopen(buffer, "w");
  sprintf(buffer, "ssfr_%s.dat", argv[2]);
  FILE *ssfr_f = check_fopen(buffer, "w");
  sprintf(buffer, "sfr_sm0_%s.dat", argv[2]);
  FILE *sfr_sm0_f = check_fopen(buffer, "w");
  for (m=9; m<15.1; m+=0.05) {
    j = (m-M_MIN-INV_BPDEX*0.5)*BPDEX;
    fm = (m-M_MIN)*BPDEX - j;
    for (t=0; t<max_t; t+=0.007) {
      float tscale = 1.0/(t+1.0); //pow(10, -t);
      for (i=num_outputs-2; i>=0; i--)
	if (steps[i].scale < tscale) break;
      ft = (tscale - steps[i].scale)/(steps[i+1].scale-steps[i].scale);
      float scale = 1.0/tscale;
      float max_m = m_at_a(15.5, tscale);
      float sfr = biterp(steps[i].sfr[j], steps[i+1].sfr[j],
			 steps[i].sfr[j+1], steps[i+1].sfr[j+1],
			 ft, fm);
      float sm = biterp(steps[i].sm_avg[j], steps[i+1].sm_avg[j],
			 steps[i].sm_avg[j+1], steps[i+1].sm_avg[j+1],
			 ft, fm);
      float ssfr = sfr - sm;
      if (m>max_m) sfr -= 12.0*(m-max_m);
      if (!(ssfr > -12)) ssfr = -12;
      if (!(sfr >= -1)) sfr = -1;
      float sfh = biterp(steps[num_outputs-1].sm_hist[j*num_outputs + i],
			 steps[num_outputs-1].sm_hist[j*num_outputs + i+1],
			 steps[num_outputs-1].sm_hist[(j+1)*num_outputs + i],
			 steps[num_outputs-1].sm_hist[(j+1)*num_outputs + i+1],
			 ft, fm) - log10(steps[i].dt);
      float sm0 = biterp(steps[num_outputs-1].sm_avg[j],
			 steps[num_outputs-1].sm_avg[j],
			 steps[num_outputs-1].sm_avg[j+1],
			 steps[num_outputs-1].sm_avg[j+1],
			 ft, fm);
      float sfr_sm0 = sfh-sm0;
      if (!(sfr_sm0 > -12)) sfr_sm0 = -12;
      if (!(sfh >= -1)) sfh = -1;
      fprintf(sfr_f, "%f %f %f %f %"PRId64" %"PRId64" %"PRId64" %f %f %f %f\n", scale, m, sfr, ft, i, j, num_outputs,
	      steps[i].sfr[j], steps[i+1].sfr[j],
	      steps[i].sfr[j+1], steps[i+1].sfr[j+1]
	      );
      fprintf(sfh_f, "%f %f %f\n", scale, m, sfh);
      fprintf(ssfr_f, "%f %f %f\n", scale, m, ssfr);
      fprintf(sfr_sm0_f, "%f %f %f\n", scale, m, sfr_sm0);
    }
    fprintf(sfr_f, "\n");
    fprintf(sfh_f, "\n");
    fprintf(ssfr_f, "\n");
    fprintf(sfr_sm0_f, "\n");
  }

  fclose(sfr_f);
  fclose(ssfr_f);
  fclose(sfh_f);
  fclose(sfr_sm0_f);
  return 0;
}
