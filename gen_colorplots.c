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
  float max_t = log10(9);
  FILE *sfr_f, *sfh_f, *sfr_fr, *sfh_fr;
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
  INVALID(the_smf) = 0;
  calc_sfh(&the_smf);

  sprintf(buffer, "sfr_%s.dat", argv[2]);
  sfr_f = check_fopen(buffer, "w");
  sprintf(buffer, "sfr_release.dat");
  sfr_fr = check_fopen(buffer, "w");
  sprintf(buffer, "sfh_%s.dat", argv[2]);
  sfh_f = check_fopen(buffer, "w");
  sprintf(buffer, "sfh_release.dat");
  sfh_fr = check_fopen(buffer, "w");
  sprintf(buffer, "sfr_nd_release.dat");
  FILE *sfr_fr_nd = check_fopen(buffer, "w");
  fprintf(sfr_fr, "#(1+z) Log10(Halo Mass) Log10(<SFR>) Log10(Median SM) Log10(Average SM) Log10(Average SM+ICL)\n");
  fprintf(sfr_fr, "#Units are Msun, Msun/yr, and Msun, respectively.\n");
  fprintf(sfr_fr_nd, "#(1+z) Log10(Halo Mass) Log10(<SFR>) Log10(Median SM) Log10(ND)\n");
  fprintf(sfr_fr_nd, "#Units are Msun, Msun/yr, Msun, and 1/(Mpc^3 dex) respectively.\n");
  fprintf(sfh_fr, "#(1+z) Log10(Halo Mass) Log10(<SFH>)\n");
  fprintf(sfh_fr, "#Units are Msun and Msun/yr, respectively.\n");
  for (m=9; m<15.1; m+=0.05) {
    j = (m-M_MIN-INV_BPDEX*0.5)*BPDEX;
    fm = (m-M_MIN-INV_BPDEX*0.5)*BPDEX - j;
    for (t=0; t<max_t; t+=0.007) {
      float tscale = pow(10, -t);
      for (i=num_outputs-2; i>=0; i--)
	if (steps[i].scale < tscale) break;
      ft = (t - log10(1.0/steps[i].scale))/log10(steps[i].scale/steps[i+1].scale);
      float scale = 1.0/tscale;
      float max_m = m_at_a(15.5, tscale);
      float sfr = biterp(steps[i].sfr[j], steps[i+1].sfr[j],
			 steps[i].sfr[j+1], steps[i+1].sfr[j+1],
			 ft, fm);
      //      float sm1 = calc_sm_at_m(m, steps[i].smhm);
      //      float sm2 = calc_sm_at_m(m, steps[i+1].smhm);
      //      float sm = sm1 + ft*(sm2-sm1);
      float sm = biterp(steps[i].sm[j], steps[i+1].sm[j],
			    steps[i].sm[j+1], steps[i+1].sm[j+1],
			    ft, fm);
      float sm_avg = biterp(steps[i].sm_avg[j], steps[i+1].sm_avg[j],
			    steps[i].sm_avg[j+1], steps[i+1].sm_avg[j+1],
			    ft, fm);
      float icl_avg = biterp(steps[i].sm_icl[j], steps[i+1].sm_icl[j],
			     steps[i].sm_icl[j+1], steps[i+1].sm_icl[j+1],
			     ft, fm);
      float nd = biterp(steps[i].t[j], steps[i+1].t[j],
                        steps[i].t[j+1], steps[i+1].t[j+1],
                        ft, fm);
      nd += log10(BPDEX);      
      if (m>max_m) sfr = -1000;
      if (!isfinite(sfr)) sfr = -1000;
      if (!isfinite(sm_avg)) sm_avg = -1000;
      if (!isfinite(icl_avg)) icl_avg = -1000;
      float sm_p_icl = log10(pow(10, sm_avg) + pow(10, icl_avg));

      fprintf(sfr_fr, "%f %f %f %f %f %f\n", scale, m, sfr, sm, sm_avg, sm_p_icl);
      if (isfinite(nd) && sfr > -1000)
	fprintf(sfr_fr_nd, "%f %f %f %f %f\n", scale, m, sfr, sm, nd);
      if (!(sfr >= -1)) sfr = -1;
      float sfh = biterp(steps[num_outputs-1].sm_hist[j*num_outputs + i],
			 steps[num_outputs-1].sm_hist[j*num_outputs + i+1],
			 steps[num_outputs-1].sm_hist[(j+1)*num_outputs + i],
			 steps[num_outputs-1].sm_hist[(j+1)*num_outputs + i+1],
			 ft, fm) - log10(steps[i].dt);
      if (!isfinite(sfh)) sfh = -1000;
      fprintf(sfh_fr, "%f %f %f\n", scale, m, sfh);
      if (!(sfh >= -1)) sfh = -1;
      fprintf(sfr_f, "%f %f %f %f %"PRId64" %"PRId64" %"PRId64" %f %f %f %f\n", scale, m, sfr, ft, i, j, num_outputs,
	      steps[i].sfr[j], steps[i+1].sfr[j],
	      steps[i].sfr[j+1], steps[i+1].sfr[j+1]
	      );
      fprintf(sfh_f, "%f %f %f\n", scale, m, sfh);
    }
    fprintf(sfr_f, "\n");
    fprintf(sfh_f, "\n");
  }

  fclose(sfr_f);
  fclose(sfh_f);
  fclose(sfr_fr);
  fclose(sfh_fr);
  fclose(sfr_fr_nd);
  return 0;
}
