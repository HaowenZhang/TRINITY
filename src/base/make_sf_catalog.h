#ifndef _MAKE_SF_CATALOG_H_
#define _MAKE_SF_CATALOG_H_
#include <inttypes.h>
#include "all_smf.h"

#define MAX_THREADS 128
#define H0 0.7

#define OBS_SIGMA_Q 0.3
#define OBS_SIGMA_SF 0.2

#define Q_SSFR -11.8
#define INTR_SCATTER 0.2
#define MAX_OUTPUTS 300

extern struct timestep *steps;
extern int64_t num_outputs;

struct packed_halo {
  int64_t id, descid;
  float pos[6], m, v, mp, vp, r, sat, vmp;
  float acc_q, acc_qne, a_half, a_4p, a_12, c, pc;
  float peak_scale, acc_scale, spin_acc, va, a_v05;
};

struct catalog_halo {
  int64_t id, descid;
  float pos[6], m, v, mp, vp, r, sat, vmp;
  float acc_q, acc_qne, a_half, a_4p, a_12, c, pc, age;
  float peak_scale, acc_scale, spin_acc, va, a_v05;
  float sm, icl, sfr, q, qf, prog_mp, obs_sfr;
  float mf, smf;
};

void read_params(char *buffer, double *params, int max_n);
void load_halos(char *filename, int num_inputs);
int64_t id_to_index(int64_t id);
void calc_losses(float scale, int64_t num_inputs);
void calc_sm_histories(int num_old_inputs);
void quench_galaxies(int64_t n);
void quench_galaxies_sm(int64_t n);
void calc_sf_from_accretion_rates(int64_t n);
void write_catalog(int64_t n);
float quenched_fraction(float sm, float sm0);
void calc_quenched_fraction(struct smf_fit *t, float scale);
int sort_by_age(const void *a, const void *b);


#endif /* _MAKE_SF_CATALOG_H_ */
