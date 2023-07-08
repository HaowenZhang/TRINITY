#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include "../base/observations.h"
#include "../base/smf.h"
#include "../base/all_smf.h"
#include "../base/distance.h"
#include "../base/integrate.h"
#include "../base/mlist.h"
#include "../base/calc_sfh.h"
#include "../base/expcache2.h"
#include "../base/param_default.h"

extern double param_default[];

#define MASS_START 7
#define MASS_STOP 12.5
#define MASS_BPDEX 10
#define MASS_BINS (int)((MASS_STOP-MASS_START)*(MASS_BPDEX))
#define MASS_STEP (1.0/(double)MASS_BPDEX)

float integrate_qf(float z_low, float z_high, double m, struct smf_fit *fit) {
  float smf_val, epsilon;
  double v_high = comoving_volume(z_high);
  double v_low = comoving_volume(z_low);
  double weight = fabs(v_high - v_low);

  if (z_low != z_high) {
      epsilon = chi2_err_helper_qf((v_high+v_low)/2.0, &m)*weight*1e-5;
      //if (PHI_HIGH_Z < z_high) epsilon *= 1e1;
      smf_val = adaptiveSimpsons(chi2_err_helper_qf, &m,
				 v_low, v_high, epsilon, 10);
      smf_val /= weight;
  }
  else {
    smf_val = chi2_err_helper_qf(v_low, &m);
  }
  return smf_val;
}

int main(int argc, char **argv)
{
  float z_low, z_high, m;
  struct smf_fit smfs[4];
  float qf_points[4][MASS_BINS];
  int i,j;
  // if (argc<4+NUM_PARAMS) {
  //   fprintf(stderr, "Usage: %s z_low z_high mass_cache (mcmc output)\n", argv[0]);
  //   exit(1);
  // }
  if (argc<4) {
    fprintf(stderr, "Usage: %s z_low z_high mass_cache (mcmc output)\n", argv[0]);
    exit(1);
  }
  z_low = atof(argv[1]);
  z_high = atof(argv[2]);



  // for (i=0; i<NUM_PARAMS; i++)
  //   smfs[0].params[i] = atof(argv[i+4]);

  // Read in the model parameter values if provided by the user
  if (argc >= 4+NUM_PARAMS)
    for (i=0; i<NUM_PARAMS; i++)
      smfs[0].params[i] = atof(argv[i+4]);
  // Otherwise use our default values
  else
    for (i=0; i<NUM_PARAMS; i++)
      smfs[0].params[i] = param_default[i];


  smfs[0].params[NUM_PARAMS] = 0;

  smfs[1] = smfs[2] = smfs[3] = smfs[0];

  /*  KAPPA(smfs[0]) = MU(smfs[0]) = SCATTER(smfs[0]) = 0;
  KAPPA(smfs[1]) = MU(smfs[1]) = 0;
  KAPPA(smfs[2]) = MU(smfs[2]) = 0;*/
  gsl_set_error_handler_off();
  setup_psf(1);
  load_mf_cache(argv[3]);
  init_timesteps();
  assert_model(smfs + 3);  
  calc_sfh(smfs + 3);
  for (j=0; j<4; j++) {
    /*for (i=0; i<M_BINS; i++) {
      if (j==3) printf("%f 0 0 0 %f\n", log10fc(steps[num_outputs-1].sm[i]),
		       log10fc(steps[num_outputs-1].sm_nd[i]));
		       }*/
    for (i=0; i<MASS_BINS; i++) {
      m = MASS_START + i*MASS_STEP;
      qf_points[j][i] = integrate_qf(z_low, z_high, m, smfs+j);
    }
    if (j==1) setup_psf(1);
  }

  for (i=0; i<MASS_BINS; i++) {
    m = MASS_START + i*MASS_STEP;
    printf("%f %f %f %f %f\n", m , qf_points[0][i], qf_points[1][i], qf_points[2][i], qf_points[3][i]); 
  }

  return 0;
}
