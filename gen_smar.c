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

#define MASS_START 7.5
#define MASS_STOP 12.5
#define MASS_BPDEX 10
#define MASS_BINS (int)((MASS_STOP-MASS_START)*(MASS_BPDEX))
#define MASS_STEP (1.0/(double)MASS_BPDEX)

double calc_m_at_sm(double sm, double m, struct smf c) {
  double sm2 = calc_sm_at_m(m, c);
  double dm_max = 5*INV_BPDEX;
  double dm = dm_max/10.0;
  if (m<=0) {
    m = M_MIN;
    dm_max = 10;
  }
  int64_t count = 0;
  while (fabs(sm2-sm)>0.0001 && count<30) {
    double sm3 = calc_sm_at_m(m+dm,c);
    double f = (sm - sm2)/(sm3-sm2);
    double move_m = f*dm;
    if (fabs(move_m) > dm_max)
      return m;
    dm = move_m / 10.0;
    m+=move_m;
    sm2 = calc_sm_at_m(m,c);
    count++;
  }
  if (count>=25) return -1;
  return m;
}

int main(int argc, char **argv)
{
  float z, m;
  struct smf_fit the_smf;
  int i;
  if (argc<3+NUM_PARAMS) {
    fprintf(stderr, "Usage: %s m mass_cache (mcmc output)\n", argv[0]);
    exit(1);
  }
  m = atof(argv[1]);
  for (i=0; i<NUM_PARAMS; i++)
    the_smf.params[i] = atof(argv[i+3]);
  the_smf.params[NUM_PARAMS] = 0;

  setup_psf(1);
  load_mf_cache(argv[2]);
  init_timesteps();
  calc_sfh(&the_smf);
  
  float hm = m + 1.5;
  for (i=0; i<150; i++) {
    z = 0.001 + i*0.1;
    struct smf smf = smhm_at_z(z, the_smf);
    hm = calc_m_at_sm(m, hm, smf);
    float mnow = m_now(1.0/(1.0+z), hm);
    float mar = ma_rate_avg(mnow, 1.0/(1.0+z));
    printf("%f %g\n", z, mar/pow(10, hm));
  }

  return 0;
}
