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

#define UV_START -25 
#define UV_STOP -15
#define MAG_BPMAG 10
#define MAG_BINS (int)((UV_STOP-UV_START)*(MAG_BPMAG))
#define MAG_STEP (1.0/(double)MAG_BPMAG)

float integrate_uvlf(float z_low, float z_high, double uv, struct smf_fit *fit) {
  float uvlf_val, epsilon;
  double v_high = comoving_volume(z_high);
  double v_low = comoving_volume(z_low);
  double weight = fabs(v_high - v_low);

  if (z_low != z_high) {
      uvlf_val = chi2_err_helper_uv((v_high+v_low)/2.0, &uv);
      //epsilon = chi2_err_helper_uv((v_high+v_low)/2.0, &uv)*weight*1e-5;
      ////if (PHI_HIGH_Z < z_high) epsilon *= 1e1;
      //uvlf_val = adaptiveSimpsons(chi2_err_helper_uv, &uv,
//				 v_low, v_high, epsilon, 10);
      //uvlf_val /= weight;
  }
  else {
    uvlf_val = chi2_err_helper_uv(v_low, &uv);
  }
  return (uvlf_val ? log10f(uvlf_val) : -15);
}

int main(int argc, char **argv)
{
  float z_low, z_high, uv;
  struct smf_fit smfs[4];
  float uvlf_points[4][MAG_BINS];
  int i,j;
  if (argc<4+NUM_PARAMS) {
    fprintf(stderr, "Usage: %s z_low z_high mass_cache (mcmc output)\n", argv[0]);
    exit(1);
  }
  z_low = atof(argv[1]);
  z_high = atof(argv[2]);
  for (i=0; i<NUM_PARAMS; i++)
    smfs[0].params[i] = atof(argv[i+4]);
  smfs[0].params[NUM_PARAMS] = 0;

  smfs[1] = smfs[2] = smfs[3] = smfs[0];

  /*  KAPPA(smfs[0]) = MU(smfs[0]) = SCATTER(smfs[0]) = 0;
  KAPPA(smfs[1]) = MU(smfs[1]) = 0;
  KAPPA(smfs[2]) = MU(smfs[2]) = 0;*/
  setup_psf(1);
  load_mf_cache(argv[3]);
  init_timesteps();
  gsl_set_error_handler_off();
  assert_model(smfs + 3);  
  //all_smf_chi2_err(smfs[3]);
  calc_sfh(smfs + 3);
  for (j=0; j<4; j++) {
    /*for (i=0; i<M_BINS; i++) {
      if (j==3) printf("%f 0 0 0 %f\n", log10fc(steps[num_outputs-1].sm[i]),
		       log10fc(steps[num_outputs-1].sm_nd[i]));
		       }*/
    for (i=0; i<MAG_BINS; i++) {
      uv = UV_START + i*MAG_STEP;
      uvlf_points[j][i] = integrate_uvlf(z_low, z_high, uv, smfs+j);
    }
    if (j==1) setup_psf(1);
  }

  for (i=0; i<MAG_BINS; i++) {
    uv = UV_START + i*MAG_STEP;
    printf("%f %f %f %f %f\n", uv , uvlf_points[0][i], uvlf_points[1][i], uvlf_points[2][i], uvlf_points[3][i]); 
  }

  return 0;
}
