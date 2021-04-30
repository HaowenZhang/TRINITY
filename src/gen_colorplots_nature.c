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
#include "universe_time.h"

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
  FILE *sfr_f, *sfr_nd_f;
  char buffer[1024];

  if (argc<3+NUM_PARAMS) {
    fprintf(stderr, "Usage: %s mass_cache (mcmc output)\n", argv[0]);
    exit(1);
  }
  for (i=0; i<NUM_PARAMS; i++)
    the_smf.params[i] = atof(argv[i+2]);
  INVALID(the_smf) = 0;

  setup_psf(1);
  load_mf_cache(argv[1]);
  init_timesteps();
  calc_sfh(&the_smf);

  sprintf(buffer, "sfr_nature.dat");
  sfr_f = check_fopen(buffer, "w");
  sprintf(buffer, "sfr_nd_nature.dat");
  sfr_nd_f = check_fopen(buffer, "w");
  sprintf(buffer, "halo_nature.dat");
  FILE *halo_f = check_fopen(buffer, "w");

  for (m=8.05; m<15.1; m+=0.05) {
    float total_stars = 0;
    j = (m-M_MIN-INV_BPDEX*0.5)*BPDEX;
    fm = (m-M_MIN)*BPDEX - j;
    for (t=1; t<(num_outputs-1)*3; t++) {
      //float tscale = pow(10, -t);
      ft = (float)(((int)t)%3) / 3.0;
      i = t/3;
      float tscale = steps[i].scale + ft*(steps[i+1].scale - steps[i].scale);
      //float dt = fabs(scale_to_years(pow(10, -(t+0.007/2.0))) - scale_to_years(pow(10, -(t-0.007/2.0))));
      float dt = steps[i].dt/3.0;
      /*for (i=num_outputs-2; i>=0; i--)
	if (steps[i].scale < tscale) break;*/
      //ft = 0; //(t - log10(1.0/steps[i].scale))/log10(steps[i].scale/steps[i+1].scale);
      float scale = 1.0/tscale;
      float sm_max = steps[i].smhm.sm_max;
      if (m==9)
	printf("%f %f\n", tscale, sm_max);
      float max_m = m_at_a(15.5, tscale);
      float sfr = biterp(steps[i].sfr[j], steps[i+1].sfr[j],
			 steps[i].sfr[j+1], steps[i+1].sfr[j+1],
			 ft, fm);
      float nd = biterp(steps[i].t[j], steps[i+1].t[j],
			steps[i].t[j+1], steps[i+1].t[j+1],
			ft, fm);
      nd += log10(BPDEX);
      if (sfr >= -10 && nd >= -12)
	total_stars += dt*pow(10, sfr+nd);
      float sfr_p_nd = sfr+nd;
      if (m>max_m) sfr -= 12.0*(m-max_m);
      if (!(sfr >= 0)) sfr = 0;
      if (!(sfr_p_nd >= -4)) sfr_p_nd = -4;
      fprintf(sfr_f, "%f %f %f\n", scale, m, sfr);
      fprintf(sfr_nd_f, "%f %f %f\n", scale, m, sfr_p_nd);
    }
    fprintf(sfr_f, "\n");
    fprintf(sfr_nd_f, "\n");
    fprintf(halo_f, "%f %f\n", m, total_stars);
  }

  fclose(sfr_f);
  fclose(sfr_nd_f);
  fclose(halo_f);

  sprintf(buffer, "sfr_sm_nature.dat");
  sfr_f = check_fopen(buffer, "w");
  sprintf(buffer, "sfr_sm_nd_nature.dat");
  sfr_nd_f = check_fopen(buffer, "w");
  sprintf(buffer, "sfr_sm_nd_errstrip.dat");
  FILE *galaxy_f = check_fopen("galaxy_nature.dat", "w");
  FILE *madau_f = check_fopen("madau_nature.dat", "w");
  FILE *sfr_nd_errstrip_f = check_fopen(buffer, "w");
  float last_nd = 1;
  for (t=1; t<(num_outputs-1)*3; t++) {
    float total_stars = 0;
    ft = (float)(((int)t)%3) / 3.0;
    i = t/3;
    float tscale = steps[i].scale + ft*(steps[i+1].scale - steps[i].scale);
    float sscale = 1.0/tscale;
    for (m=6.0; m<12.1; m+=0.03) {
      float nd = evaluate_from_step(1.0/tscale - 1.0, m, 1);
      nd = (nd > 0) ? log10(nd) : -10;
      if (nd < -8 && last_nd > -8) {
	fprintf(sfr_nd_errstrip_f, "%f %f %f %f\n", sscale, 12.1, 0.0, 12.1-m);
      }
      last_nd = nd;
    }
    float dt = steps[i].dt/3.0;
    total_stars = calc_cosmic_sfr(1.0/tscale-1.0);
    fprintf(madau_f, "%f %f\n", sscale, total_stars);
  }
  fclose(madau_f);

  for (m=6.0; m<12.1; m+=0.03) {
    float total_stars = 0;
    for (t=1; t<(num_outputs-1)*3; t++) {
      ft = (float)(((int)t)%3) / 3.0;
      i = t/3;
      float tscale = steps[i].scale + ft*(steps[i+1].scale - steps[i].scale);
      float dt = steps[i].dt/3.0;
      float mu = steps[i].smhm.mu + ft*(steps[i+1].smhm.mu-steps[i].smhm.mu);
      double sfr = calc_ssfr(m, 1.0/tscale - 1.0)*pow(10, m-mu);
      float sscale = 1.0/tscale;
      float nd = evaluate_from_step(1.0/tscale - 1.0, m, 1);
      if (sfr >= 1e-10 && nd >= 1e-12)
	total_stars += dt*sfr*nd;
      nd = (nd > 0) ? log10(nd) : -10;
      sfr = (sfr > 0) ? log10(sfr) : 0;
      if (nd < -8.5) sfr = -10;
      float sfr_p_nd = sfr+nd;
      if (!(sfr >= 0)) sfr = 0;
      if (!(sfr_p_nd >= -4)) sfr_p_nd = -4;
      fprintf(sfr_f, "%f %f %f\n", sscale, m, sfr);
      fprintf(sfr_nd_f, "%f %f %f\n", sscale, m, sfr_p_nd);
    }
    fprintf(sfr_f, "\n");
    fprintf(sfr_nd_f, "\n");
    fprintf(galaxy_f, "%f %f\n", m, total_stars);
  }

  fclose(sfr_nd_errstrip_f);
  fclose(sfr_f);
  fclose(sfr_nd_f);
  fclose(galaxy_f);
  return 0;
}
