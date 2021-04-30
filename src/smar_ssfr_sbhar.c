#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "observations.h"
#include "smf.h"
#include "mah.h"
#include "all_smf.h"
#include "distance.h"
#include "integrate.h"
#include "mlist.h"
#include "calc_sfh.h"
#include "expcache2.h"

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

int main(int argc, char **argv)
{
  int64_t i;
  struct smf_fit smf;
  //double mh, mstar, mbh;
  double mh = 12.1;
  int mb = (mh - M_MIN) * BPDEX - 0.5;
  //mstar = 9;
  //mbh = 4;
  if (argc<2+NUM_PARAMS) {
    fprintf(stderr, "Usage: %s mass_cache (mcmc output)\n", argv[0]);
    exit(1);
  }
  //double z = atof(argv[1]);
  for (i=0; i<NUM_PARAMS; i++)
    smf.params[i] = atof(argv[i+2]);
  gsl_set_error_handler_off();
  setup_psf(1);
  load_mf_cache(argv[1]);
  init_timesteps();
  INVALID(smf) = 0;
  calc_sfh(&smf);
  int64_t step;
  double f;
  //calc_step_at_z(z, &step, &f);
  
  
  printf("#z Mh SMAR SM SSFR Mbh SBHAR\n");


  for (i=0; i<num_outputs; i++) 
  {
    double mbh = steps[i].bh_mass_avg[mb];
    double mstar = steps[i].sm_avg[mb];
    double bhar = steps[i].bh_acc_rate[mb];
    double sfr = steps[i].sfr[mb];
    printf("%f %f %e %e %e %e %e\n", 1 / steps[i].scale - 1.0, mh, ma_rate_avg_mnow(mh, steps[i].scale) / exp10(mh), 
                  mstar, sfr/mstar,
                  mbh, bhar/mbh);
  }
  // fclose(pfile);
  return 0;
}
