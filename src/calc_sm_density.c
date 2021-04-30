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
#include <string.h>

extern int64_t num_outputs;
extern struct timestep *steps;
struct obs_smf_point *cur_smfs = NULL;
extern int64_t num_points;

struct int_helper {
  float z1, z2;
};

int sort_smf_points(const void *a, const void *b) {
  const struct obs_smf_point *c = a;
  const struct obs_smf_point *d = b;
  if (c->mass < d->mass) return -1;
  if (c->mass > d->mass) return 1;
  return 0;
}

double eval_single_point(struct obs_smf_point *cur_smf) {
  double smf_val, epsilon;
  double v_high = comoving_volume(cur_smf->z_high);
  double v_low = comoving_volume(cur_smf->z_low);
  double weight = v_high - v_low;
  double m;

  m = cur_smf->mass;
  if (cur_smf->type == SMF_TYPE) {
    if (cur_smf->z_low != cur_smf->z_high) {
      if (cur_smf->z_low < 0.2) {
	//epsilon *= 0.01;
	smf_val = adaptiveGauss(chi2_err_helper, &m, v_low, v_high,
				PHI_INTEGRAL_PRECISION*0.01, 1);
      }
      else {
	epsilon = chi2_err_helper((v_high+v_low)/2.0, &m)*weight*PHI_INTEGRAL_PRECISION;
	smf_val = adaptiveSimpsons(chi2_err_helper, &m,
				   v_low, v_high, epsilon, 10);
      }
      smf_val /= weight;
    }
    else {
      smf_val = chi2_err_helper(v_low, &m);
    }
  }
  else if (cur_smf->type == SSFR_TYPE) {
    smf_val = calc_ssfr(m, cur_smf->z_low);
  }
  else { //if (cur_smf->type == CSFR_TYPE) {
    smf_val = calc_cosmic_sfr(cur_smf->mass); //Actually redshift for CSFR
  }
  return smf_val;
}

double sm_density_helper(double sm, void *extra_data) {
  struct int_helper *ih = extra_data;
  struct obs_smf_point tp = {0};
  tp.mass = sm;
  tp.type = SMF_TYPE;
  tp.z_low = ih->z1;
  tp.z_high = ih->z2;
  return (pow(10, sm)*eval_single_point(&tp));
}


int main(int argc, char **argv)
{
  struct smf_fit the_smf;
  struct int_helper the_helper;
  int i;
  char buffer[1024] = {0};

  if (argc<3+NUM_PARAMS) {
    fprintf(stderr, "Usage: %s mass_cache (mcmc output)\n", argv[0]);
    exit(1);
  }
  for (i=0; i<NUM_PARAMS; i++)
    the_smf.params[i] = atof(argv[i+2]);

  setup_psf(1);
  load_mf_cache(argv[1]);
  init_timesteps();
  INVALID(the_smf) = 0;
  calc_sfh(&the_smf);
  while (fgets(buffer, 1024, stdin)) {
    buffer[strlen(buffer)-1] = 0;
    num_points = 0;
    load_real_smf(&cur_smfs, &num_points, buffer);
    qsort(cur_smfs, num_points, sizeof(struct obs_smf_point), sort_smf_points);
    float m_min = cur_smfs[0].mass;
    float m_max = cur_smfs[num_points-1].mass;
    float hbin_size = (m_max-m_min)*((double)num_points)/((double)num_points-1.0)/2.0;
    m_min -= hbin_size;
    m_max += hbin_size;
    the_helper.z1 = cur_smfs[0].z_low;
    the_helper.z2 = cur_smfs[0].z_high;
    float total_sm1 = adaptiveGauss(sm_density_helper, &the_helper, 6, m_min,
				  PHI_INTEGRAL_PRECISION*0.01, 0);
    float total_sm2 = adaptiveGauss(sm_density_helper, &the_helper, m_min, m_max,
				  PHI_INTEGRAL_PRECISION*0.01, 0);
    float total_sm3 = adaptiveGauss(sm_density_helper, &the_helper, m_max, 12.5,
				  PHI_INTEGRAL_PRECISION*0.01, 0);
    int64_t j, j2;
    float total_sm_obs = 0;
    float fake_sm_obs = 0;
    for (j=0; j<num_points; j++) {
      total_sm_obs += pow(10, cur_smfs[j].mass+cur_smfs[j].val);
      fake_sm_obs += sm_density_helper(cur_smfs[j].mass, &the_helper);
    }
    total_sm_obs *= total_sm2 / fake_sm_obs;
    float v_avg = 0.5*(comoving_volume(cur_smfs[0].z_high)+
		       comoving_volume(cur_smfs[0].z_low));
    float z_avg = comoving_volume_to_redshift(v_avg);
    float total_sm = total_sm1 + total_sm2 + total_sm3;
    printf("%f %e %e\n", z_avg, total_sm, (total_sm_obs / total_sm2)*total_sm);
  }
  return 0;
}