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




float fq[M_BINS*10] = {0};
float mf[M_BINS*10] = {0};
float smm[M_BINS*10] = {0};
float shortfall[M_BINS*10] = {0};
float shortfall_norm[M_BINS*10] = {0};

float quenched_sm0(float scale) {
  return pow(10, 10.2+0.5*(1.0/scale-1.0)); //Default
  //return pow(10, 10.4+0.5*(1.0/scale-1.0)); //H&W
}

//Modified!!!
float quenched_fraction(float sm, float sm0) {
  if (!(sm>0)) return 0;
  return 1.0/(pow(sm/sm0, -0.5)+1.0);
}


int main(int argc, char **argv)
{
  struct smf_fit smf;
  int i;
  if (argc<3+NUM_PARAMS) {
    fprintf(stderr, "Usage: %s mass_cache niter (mcmc output)\n", argv[0]);
    exit(1);
  }
  for (i=0; i<NUM_PARAMS; i++)
    smf.params[i] = atof(argv[i+3]);
  smf.params[NUM_PARAMS] = 0;
  INVALID(smf) = 0;
  setup_psf(1);
  load_mf_cache(argv[1]);
  init_timesteps();
  calc_sfh(&smf);
  int niter = atol(argv[2]);

  float z;
  for (z=0; z<0.01; z+=0.1) {
    struct smf c = smhm_at_z(z, smf);
    float scale = 1.0/(1.0+z);
    float sm0 = quenched_sm0(scale);
    for (i=0; i<M_BINS*10; i++) {
      float m = M_MIN + i*INV_BPDEX/10.0;
      mf[i] = pow(10, mf_cache(scale, m));
      smm[i] =calc_sm_at_m(m, c);
      fq[i] = quenched_fraction(pow(10, smm[i]), sm0);
    }

    float sm = 0;
    float scatter = sqrt(c.scatter*c.scatter + c.obs_scatter+c.obs_scatter);
    int j=0, max_j = niter;

    for (j=0; j<max_j+1; j++) {    
      for (i=0; i<M_BINS*10; i++) shortfall[i] = shortfall_norm[i] = 0;

      for (sm = 8; sm <12; sm+=0.05) {
	double tot=0, totq=0;
	for (i=0; i<M_BINS*10; i++) {
	  float dm = (smm[i]-sm)/scatter;
	  float nd = mf[i]*exp(-0.5*dm*dm);
	  tot += nd;
	  totq += nd*fq[i];
	}
	float model_qf = totq/tot;
	float actual_qf = quenched_fraction(pow(10, sm), sm0);
	if (j==max_j)
	  printf("%f %f %f\n", sm, totq/tot, quenched_fraction(pow(10, sm), sm0));
	
	for (i=0; i<M_BINS*10; i++) {
	  float dm = (smm[i]-sm)/scatter;
	  float nd = exp(-0.5*dm*dm);
	  shortfall[i] += nd*(actual_qf-model_qf);
	  shortfall_norm[i] += nd;
	}
      }
      
      if (j<max_j) {
	for (i=0; i<M_BINS*10; i++) {
	  fq[i] += shortfall[i]/shortfall_norm[i];
	  if (fq[i] > 1) fq[i] = 1;
	  if (fq[i] < 0) fq[i] = 0;
	}
      }
    }
  }

  /*for (i=0; i<M_BINS*10; i++) {
    float m = M_MIN + i*INV_BPDEX/10.0;
    printf("%f %f\n", m, fq[i]);
    }*/

  return 0;
}
