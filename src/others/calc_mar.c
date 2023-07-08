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

#define NUM_ZS 9
extern int64_t num_outputs;
extern struct timestep *steps;

float m_to_bin(float m) {
  if (m<M_MIN) return 0;
  if (m>(M_MAX-INV_BPDEX)) return (M_BINS-2);
  return ((m-M_MIN)*BPDEX);
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
  float mar1 = mar_from_mbins(n,j);
  float mar2 = mar_from_mbins(n+1,j);
  return (0.5*(mar1+mar2));
}

int main(int argc, char **argv)
{
  float m;
  //float *zs = {0.1, 1, 2, 3, 4, 5, 6, 7, 8};
  struct smf_fit the_smf;
  int i, j;
  FILE *sfr_f, *sfh_f;

  if (argc<3) {
    fprintf(stderr, "Usage: %s mass_cache mass\n", argv[0]);
    exit(1);
  }

  float m = atof(argv[2]);

  setup_psf(1);
  load_mf_cache(argv[1]);
  init_timesteps();
  
  for (j=0; j<M_BINS; j++) {
    m = M_MIN + j*INV_BPDEX;
    if (m < 9 || m > 16) continue;
    for (i=0; i<num_outputs; i++) {
      float scale = 1.0/steps[i].scale;
      float sfr = steps[i].sfr[j];
      if (sfr < 1) sfr = 0;
      else sfr = log10(sfr);
      float sfh = steps[num_outputs-1].sm_hist[j*num_outputs + i]/steps[i].dt;
      if (sfh < 1) sfh = 0;
      else sfh = log10(sfh);
      fprintf(sfr_f, "%f %f %f\n", scale, m, sfr);
      fprintf(sfh_f, "%f %f %f\n", scale, m, sfh);
    }
    fprintf(sfr_f, "\n");
    fprintf(sfh_f, "\n");
  }

  fclose(sfr_f);
  fclose(sfh_f);
  return 0;
}
