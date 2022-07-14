#ifndef SMASS_CV_COSMO_H
#define SMASS_CV_COSMO_H
#include "all_smf.h"
#include "obs_smf.h"

#define SM_START 8
#define SM_END 12
#define SM_STEP 0.125
#define MAX_RESULTS ((int)((SM_END-SM_START)/SM_STEP+2) * 10)


#define Z_CV_LIMIT 0.2
#define SCALE_CV_LIMIT (1.0/(Z_CV_LIMIT+1.0))

#define MASS_CONVERT_MAX_MASS 17
#define MASS_CONVERT_MIN_MASS 3
#define MASS_CONVERT_STEPS_PER_DEX 100
#define MASS_CONVERT_INV_STEPS_PER_DEX (1.0/(double)MASS_CONVERT_STEPS_PER_DEX)
#define MASS_CONVERT_SIZE ((MASS_CONVERT_MAX_MASS-MASS_CONVERT_MIN_MASS) \
			   * MASS_CONVERT_STEPS_PER_DEX + 1)

#define FITTING_ROUNDS 4
#define CHI2_LIMIT 1e-2
#define STEP_LIMIT 8
#define MIN_FIT_LIMIT 10
#define STEP_RATIO (1/2.0) //How quickly step sizes shrink

extern int no_matching_scatter;
extern int no_systematics;
extern int no_obs_scatter;
extern struct smf (*model_cb)(float, void *);
extern int higher_z;


void read_params(char *buffer, double *params, int max_n);
int float_compare(const void *a, const void *b);
int calculate_smhm_relation(struct smf_fit fit, float *results, int max_calcs,
			     float *scales, int num_scales, int convert_mass);

struct cosmology {
  float mstar;
  struct z_mf_cache *mf;
};

struct cosmic_variance {
  float dphi, dalpha, dmstar;
};

float mass_convert(float m);
void find_mstar(struct cosmology *c, float scale);
void gen_mass_convert_table(struct cosmology *c1, struct cosmology *c2,
			    float scale, float dphi, float dmstar, 
			    float dalpha);
void load_cvars(char *filename);
float find_m(float sm, struct smf smf, int convert_mass);
void fitting_round(double *params, double *steps);


#endif /* SMASS_CV_COSMO */
