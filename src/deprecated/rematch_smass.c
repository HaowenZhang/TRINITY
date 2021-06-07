#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "obs_smf.h"
#include "smf.h"
#include "all_smf.h"
#include "distance.h"
#include "integrate.h"
#include "sm_limits.h"

#define MASS_START 8
#define MASS_STOP 12.5
#define MASS_BPDEX 100
#define MASS_BINS (int)((MASS_STOP-MASS_START)*(MASS_BPDEX))
#define MASS_STEP (1.0/(double)MASS_BPDEX)
#define HM_MAX 17
#define HM_MIN 8
#define HM_BPDEX 200
#define HM_BINS (2+(int)(HM_MAX-HM_MIN)*(HM_BPDEX))

double hm_nd[HM_BINS];

float integrate_smf(float z_low, float z_high, float m, struct smf_fit *fit) {
  struct chi2_err_helper_data helper;
  float smf_val;
  helper.cb = linear_smf_cb;
  helper.extra_data = (void *) fit;
  helper.m = m;
  smf_val = chi2_err_helper(comoving_volume(z_low), &helper);
  return (smf_val);
}

void gen_hm_nd_cache(float scale) {
  int i;
  float m;
  for (i=0; i<HM_BINS; i++) {
    m = HM_MAX - ((double)i)/((double)HM_BPDEX);
    hm_nd[i] = pow(10,mf_cache(scale, m))/((double)HM_BPDEX);
    if (i) hm_nd[i]+=hm_nd[i-1];
  }
  for (i=0; i<HM_BINS; i++)
    hm_nd[i] = (hm_nd[i]>0) ? log10(hm_nd[i]) : -1000;
}

float hm_nd_cache(float m) {
  float f = ((float)HM_MAX - m)*(float)HM_BPDEX;
  int i = f;
  if (m>HM_MAX) return -1000;
  if (m<=HM_MIN) return (hm_nd[HM_BINS-1]);
  f -= i;
  return (hm_nd[i] + f*(hm_nd[i+1]-hm_nd[i]));
}

float nd_to_halo_mass(float nd) {
  float m = 13.5;
  float step = 0.5;
  float trial_nd = hm_nd_cache(m);
  while (trial_nd < nd && m > HM_MIN) { m--; trial_nd = hm_nd_cache(m); }
  while (trial_nd > nd && m < HM_MAX) { m++; trial_nd = hm_nd_cache(m); }
  while (step > 1e-6) {
    if (trial_nd < nd) m -= step;
    else m += step;
    step *= 0.5;
    trial_nd = hm_nd_cache(m);
  }
  return m;
}

int main(int argc, char **argv)
{
  float z, m, scale;
  float sm_min, sm_max;
  struct smf_fit base_mp_smf, base_m_smf, base_p_smf, base_smf;
  float bsmf[MASS_BINS], bmsmf[MASS_BINS], bmpsmf[MASS_BINS], bpsmf[MASS_BINS];
  int i;
  if (argc<3+NUM_PARAMS) {
    fprintf(stderr, "Usage: %s z mass_cache (mcmc output)\n", argv[0]);
    exit(1);
  }
  z = atof(argv[1]);
  scale = 1.0 / (z+1.0);
  sm_min = sm_min_mass(z);
  sm_max = sm_max_mass(z);
  for (i=0; i<NUM_PARAMS; i++)
    base_mp_smf.params[i] = atof(argv[i+3]);

  KAPPA(base_mp_smf) = MU(base_mp_smf) = 0;
  base_m_smf = base_mp_smf;
  base_p_smf = base_mp_smf;
  base_smf = base_mp_smf;

  SCATTER(base_m_smf) = 0;
  SCATTER(base_mp_smf) = 0;

  load_psf_caches(1);
  load_mf_cache(argv[2]);
  gen_hm_nd_cache(scale);

  for (i=0; i<MASS_BINS; i++) {
    m = MASS_START + (MASS_BINS-i)*MASS_STEP;
    bsmf[i] = integrate_smf(z, z, m, &base_smf)*MASS_STEP;
    bmsmf[i] = integrate_smf(z, z, m, &base_m_smf)*MASS_STEP;
    if (i) { bmsmf[i]+=bmsmf[i-1]; bsmf[i]+=bsmf[i-1]; }
  }

  load_psf_caches(0);

  printf("#m0p0    m1p0   m0p1   m1p1\n");
  for (i=0; i<MASS_BINS; i++) {
    m = MASS_START + (MASS_BINS-i)*MASS_STEP;
    bmpsmf[i] = integrate_smf(z, z, m, &base_mp_smf) * MASS_STEP;
    bpsmf[i] = integrate_smf(z, z, m, &base_p_smf) * MASS_STEP;
    if (i) { bpsmf[i]+=bpsmf[i-1]; bmpsmf[i]+=bmpsmf[i-1]; }
#define hm(x) nd_to_halo_mass(x[i] ? log10(x[i]) : -15)
    if (m<sm_min) continue;
    //if (m>sm_max) continue;
    printf("%f %f %f %f %f\n", m , hm(bsmf), hm(bmsmf), hm(bpsmf), hm(bmpsmf));
#undef hm
  }

  return 0;
}
