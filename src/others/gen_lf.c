#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <strings.h>
#include "observations.h"
#include "smf.h"
#include "all_smf.h"
#include "distance.h"
#include "integrate.h"
#include "mlist.h"
#include "calc_sfh.h"
#include "expcache2.h"
#include "all_luminosities.h"


#define MASS_START 7
#define MASS_STOP 16
#define MASS_BPDEX 10
#define MASS_BINS (int)((MASS_STOP-MASS_START)*(MASS_BPDEX))
#define MASS_STEP (1.0/(double)MASS_BPDEX)

int main(int argc, char **argv)
{
  float z;
  struct smf_fit the_smf;
  float *lums = NULL;
  int i,j;
  if (argc<4+NUM_PARAMS) {
    fprintf(stderr, "Usage: %s z lum_name 0(abs)/1(apparent) mass_cache (mcmc output)\n", argv[0]);
    exit(1);
  }
  z = atof(argv[1]);
  for (i=0; i<NUM_PARAMS; i++)
    the_smf.params[i] = atof(argv[i+5]);
  INVALID(the_smf) = 0;
gsl_set_error_handler_off();
  setup_psf(1);
  load_mf_cache(argv[4]);
  init_timesteps();  
  calc_sfh(&the_smf);
  load_luminosities("dust/mags_ab_zall_jwst_nodust.dat");

  enum lum_types ltype;
  for (ltype=0; ltype<LUM_TYPES; ltype++) {
    if (!strcasecmp(argv[2], lum_names[ltype])) break;
  }
  if (ltype == LUM_TYPES) {
    fprintf(stderr, "Unknown filter type %s!\n", argv[2]);
    fprintf(stderr, "Allowable types:");
    for (i=0; i<LUM_TYPES; i++) fprintf(stderr, " %s", lum_names[i]);
    fprintf(stderr, "\n");
    exit(1);
  }

  float z_src = (atof(argv[3]) == 0) ? z : -1;
  float z_obs = (atof(argv[3]) == 0) ? z : 0;
  gen_all_luminosities(ltype, &lums, z_src, z_obs);

  float t[MASS_BINS];
  float l[MASS_BINS];
  int64_t max_n = 0;
  float a = 1.0/(1.0+z);
  for (max_n=0; max_n<MASS_BINS; max_n++) {
    float m = MASS_START + max_n*MASS_STEP;
    t[max_n] = pow(10, mf_cache(a, m));
    if (t[max_n] < 3e-10) break;
    l[max_n] = lum_at_hm_z(lums, m, z);
  }
  struct smf c = smhm_at_z(z, the_smf);
  float scatter = c.scatter*2.5;
  float ds = -1.0/(2.0*scatter*scatter);

  int64_t min_l = (z_src < 0) ? 200 : -250;
  int64_t max_l = (z_src < 0) ? 320 : -120;
  for (i=min_l; i<=max_l; i++) {
    double nd = 0;
    float lum = i/10.0;
    float norm = 1.0/sqrt(2.0*M_PI*scatter*scatter);
    norm /= (double)MASS_BPDEX;
    for (j=0; j<max_n; j++) {
      double dl = lum - l[j];
      nd += t[j]*exp(dl*dl*ds);
    }
    printf("%f %e\n", lum, nd*norm);
  }
  return 0;
}
