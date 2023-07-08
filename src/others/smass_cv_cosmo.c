#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "all_smf.h"
#include "obs_smf.h"
#include "expcache.h"
#include "mt_rand.h"
#include "smass_cv_cosmo.h"


int num_cvars = 0;
struct cosmic_variance *cvars = NULL;
float mass_convert_table[MASS_CONVERT_SIZE];
struct cosmology *cosmologies;
int num_cosmos = 0;
float new_smhm[MAX_RESULTS];
float scale_list[4] = {1,0.5,0.33,0.2};
int num_scales = 2;

int main(int argc, char **argv) {
  FILE *input;
  char buffer[1024];
  struct smf_fit *fits;
  int i, j, num_entries;
  double steps[NUM_PARAMS] = {0.3, 1, 1, 0.2, 1, 0.09, 1, 1, 1, 1, 1, 0, 0, 0, 0};

  if (argc<6) {
    printf("Usage: %s model num_entries red_smf_mcmc.dat cv_file mf1 [mf2]...\n", argv[0]);
    exit(1);
  }
  
  r250_init(87L);

  i = atol(argv[1]);
  model_cb = linear_smf_cb;
  num_scales = (i==2) ? 4 : 2;
  steps[2] = (i==2) ? 1 : 0;
  no_matching_scatter = no_systematics = no_obs_scatter = 1;

  num_entries = atol(argv[2]);
  fits = (struct smf_fit *)malloc(sizeof(struct smf_fit)*num_entries);

  load_cvars(argv[4]);
  num_cosmos = argc-5;
  cosmologies = (struct cosmology *)calloc(num_cosmos,sizeof(struct cosmology));
  for (i=5; i<argc; i++) {
    all_mf_caches = NULL;
    load_mf_cache(argv[i]);
    cosmologies[i-5].mf = all_mf_caches;
  }

  if (!(input = fopen(argv[3], "r"))) {
    printf("Couldn't open file %s for reading!\n", argv[3]);
    exit(1);
  }

  for (i=0; i<num_entries; i++) {
    if (!fgets(buffer, 1024, input)) break;
    read_params(buffer, fits[i].params, NUM_PARAMS+1);
  }
  fclose(input);
  num_entries = i;

  /*  for (i=0; i<100; i++) {
    int cvar, cosmo;
    cvar = (((int)(dr250()*num_cvars))%num_cvars);
    cosmo = (((int)(dr250()*(num_cosmos-1)+1))%num_cosmos);
    find_mstar(&(cosmologies[cosmo]), 0.5);
    find_mstar(cosmologies, 0.5);
    gen_mass_convert_table(cosmologies, &(cosmologies[cosmo]),
			   0.5, cvars[cvar].dphi, cvars[cvar].dmstar,
			   cvars[cvar].dalpha);
    for (j=100; j<160; j++) {
      float m = (float)j/10.0;
      printf("%f %f\n", m, mass_convert(m)-m);
    }
    printf("\n");
  }
  exit(0);
  */

  //Build converted relation for each entry.
  for (i=0; i<num_entries; i++) {
    int cosmo = calculate_smhm_relation(fits[i], new_smhm, MAX_RESULTS,
				    scale_list, num_scales, i);
    fitting_round(fits[i].params, steps);
    fitting_round(fits[i].params, steps);
    //Adjust SM_0 for change in h0:
    if (cosmologies[0].mf->h0>0 && cosmologies[cosmo].mf->h0>0)
      EFF_0(fits[i])+=2*log10(cosmologies[0].mf->h0/cosmologies[cosmo].mf->h0);

    for (j=0; j<NUM_PARAMS; j++) {
      printf("%f ", fits[i].params[j]);
    }
    printf("%f\n", fits[i].params[NUM_PARAMS]);
  }
  return 0;
}

int calculate_smhm_relation(struct smf_fit fit, float *results, int max_calcs,
			     float *scales, int num_scales, int convert_mass)
{
  int i, j=0, cvar, cosmo, highz_cvar;
  float sm, scale;
  struct smf the_smf;
  highz_cvar = convert_mass ? (((int)(dr250()*num_cvars))%num_cvars) : 0;
  cosmo = convert_mass ? (((int)(dr250()*(num_cosmos-1)+1))%num_cosmos) : 0;
  for (i=0; i<num_scales; i++) {
    scale = scales[i];
    cvar = (scale > SCALE_CV_LIMIT) ? 0 : highz_cvar;
    if (convert_mass) {
      find_mstar(&(cosmologies[cosmo]), scale);
      find_mstar(cosmologies, scale);
      gen_mass_convert_table(cosmologies, &(cosmologies[cosmo]),
			     scale, cvars[cvar].dphi, cvars[cvar].dmstar,
			     cvars[cvar].dalpha);
    }
    the_smf = model_cb(1.0/scale - 1.0, (void*)(&fit));
    for (sm=SM_START; (sm<SM_END+SM_STEP)&&(j<max_calcs); sm+=SM_STEP,j++)
      results[j] = find_m(sm, the_smf, convert_mass);
  }
  for (; j<max_calcs; j++) results[j] = 0;
  return cosmo;
}


float mass_convert(float m) {
  float f = (m-MASS_CONVERT_MIN_MASS)*MASS_CONVERT_STEPS_PER_DEX;
  int i = f;
  f -= i;
  if (i<0) return (m-MASS_CONVERT_MIN_MASS + mass_convert_table[0]);
  if (i>(MASS_CONVERT_SIZE-2)) 
    return (m-MASS_CONVERT_MAX_MASS + mass_convert_table[MASS_CONVERT_SIZE-1]);
  return (mass_convert_table[i] + 
	  f*(mass_convert_table[i+1]-mass_convert_table[i]));
}


void find_mstar(struct cosmology *c, float scale) {
  float m;
  float s;
  float alpha;
  all_mf_caches = c->mf;
  alpha = 1.0 + (mf_cache(scale, MASS_CONVERT_MIN_MASS-0.0005) - 
		 mf_cache(scale, MASS_CONVERT_MIN_MASS+0.0005))*1000;
  for (m=MASS_CONVERT_MIN_MASS; m<MASS_CONVERT_MAX_MASS; m+=0.01) {
    s = (mf_cache(scale, m-0.0005) - mf_cache(scale, m+0.0005))*1000;
    if (s > alpha) break;
  }
  c->mstar = m;
}

//For conversion from c1 to c2
void gen_mass_convert_table(struct cosmology *c1, struct cosmology *c2,
			    float scale, float dphi, float dmstar, 
			    float dalpha) {
  float nd1[MASS_CONVERT_SIZE];
  float nd2[MASS_CONVERT_SIZE];
  double nd=0;
  float m;
  int i,j,k;

  all_mf_caches = c1->mf;
  for (i=MASS_CONVERT_SIZE-1; i>=0; i--) {
    m = i*MASS_CONVERT_INV_STEPS_PER_DEX + MASS_CONVERT_MIN_MASS;
    nd+=exp10fc(mf_cache(scale, m))*MASS_CONVERT_INV_STEPS_PER_DEX;
    nd1[i] = (nd>0) ? log10fc(nd) : -1000;
  }

  nd=0;
  all_mf_caches = c2->mf;
  for (i=MASS_CONVERT_SIZE-1; i>=0; i--) {
    m = i*MASS_CONVERT_INV_STEPS_PER_DEX + MASS_CONVERT_MIN_MASS - dmstar;
    nd+=exp10fc(mf_cache(scale, m) + dphi - dalpha*(m - c2->mstar))
      * MASS_CONVERT_INV_STEPS_PER_DEX;
    nd2[i] = (nd>0) ? log10fc(nd) : -1000;
  }

  //Find first comparable points:
  for (i=MASS_CONVERT_SIZE-1; i>=0; i--)
    if (nd1[i]>-1000) break;

  for (j=MASS_CONVERT_SIZE-1; j>=0; j--)
    if (nd2[j]>-1000) break;

  m = j*MASS_CONVERT_INV_STEPS_PER_DEX + MASS_CONVERT_MIN_MASS;
  for (k=i; k<MASS_CONVERT_SIZE; k++)
    mass_convert_table[k] = m+(k-i)*MASS_CONVERT_INV_STEPS_PER_DEX;

  if (j==MASS_CONVERT_SIZE-1) j=MASS_CONVERT_SIZE-2;
  //Step through remaining number densities:
  for (i = i-1; i>=0; i--) {
    while (j && nd2[j] < nd1[i]) { j--; }
    if (j>0) {
      mass_convert_table[i] = (j + (nd1[i]-nd2[j])/(nd2[j+1]-nd2[j])) *
	MASS_CONVERT_INV_STEPS_PER_DEX + MASS_CONVERT_MIN_MASS;
      k=i;
    } 
    else { //Beyond end of list
      mass_convert_table[i] = 
	(i-k)*MASS_CONVERT_INV_STEPS_PER_DEX + MASS_CONVERT_MIN_MASS;
    }
  }
}

void add_to_cvars(float dphi, float dmstar, float dalpha)
{
  if (!(num_cvars % 1000)) {
    cvars = (struct cosmic_variance *) 
      realloc(cvars, sizeof(struct cosmic_variance)*(num_cvars+1000));
  }
  cvars[num_cvars].dmstar = dmstar;
  cvars[num_cvars].dalpha = dalpha;
  cvars[num_cvars].dphi = dphi;
  num_cvars++;
}

void load_cvars(char *filename) {
  char buffer[1000];
  FILE *input;
  int n;
  float dmstar, dalpha, dphi;
  if (!(input = fopen(filename, "r"))) {
    printf("Unable to open file %s for reading!\n", filename);
    exit(3);
  }
  add_to_cvars(0, 0, 0);
  while (fgets(buffer, 1000, input)) {
    n = sscanf(buffer, "%f %f %f", &dphi, &dmstar, &dalpha);
    if (n != 3) continue;
    add_to_cvars(dphi, dmstar, dalpha);
  }
  fclose(input);
}

float find_m(float sm, struct smf smf, int convert_mass) {
#define hm(x) sm_to_log_halo_mass(x, smf, mul)
  float mul = smf.m_1 - 0.5; //smf.gamma*(float)(M_LN2/M_LN10);
  float m = hm(sm);
  if (convert_mass) return mass_convert(m);
  return m;
#undef hm
}


void read_params(char *buffer, double *params, int max_n) {
  int num_entries = 0;
  char *cur_pos = buffer, *end_pos;
  float val = strtod(cur_pos, &end_pos);
  while (cur_pos != end_pos && num_entries < max_n) {
    params[num_entries] = val;
    num_entries++;
    cur_pos=end_pos;
    while (*cur_pos==' ' || *cur_pos=='\t' || *cur_pos=='\n') cur_pos++;
    val = strtod(cur_pos, &end_pos);
  }
}


float calc_chi2(double *params) {
  struct smf_fit the_fit;
  float smhm[MAX_RESULTS];
  float chi2 = 0, dx;
  int i;
  for (i=0; i<NUM_PARAMS; i++)
    the_fit.params[i] = params[i];

  calculate_smhm_relation(the_fit, smhm, MAX_RESULTS,
			  scale_list, num_scales, 0);
  for (i=0; i<MAX_RESULTS; i++) {
    dx = smhm[i] - new_smhm[i];
    chi2 += dx*dx;
  }
  return chi2;
}

inline float minimum_delta(float chi2_l, float chi2, float chi2_r)
{
  float d_chi2_dx = (chi2_r - chi2_l) / 2.0;
  float d2_chi2_dx2 = (chi2_r + chi2_l - 2*chi2);
  if (chi2 > chi2_r && chi2 > chi2_l) // We are close to a local maximum
    return ((chi2_l < chi2_r) ? -1 : 1);
  return (d2_chi2_dx2 ? (-d_chi2_dx / d2_chi2_dx2) : 0);
}


float improve_fit(double *params, double *steps, int i, float chi2)
{
  float chi2_l, chi2_r, p, dx, chi2_new;
  p = params[i];
  params[i] = p + steps[i];
  chi2_r = calc_chi2(params);
  params[i] = p - steps[i];
  chi2_l = calc_chi2(params);

  dx = minimum_delta(chi2_l, chi2, chi2_r);
  if (fabs(dx) > STEP_LIMIT) dx = copysign(STEP_LIMIT, dx);
  params[i] = p + dx*steps[i];
  chi2_new = calc_chi2(params);

  if (chi2_l < chi2_new) { chi2_new = chi2_l; params[i] = p-steps[i]; }
  if (chi2_r < chi2_new) { chi2_new = chi2_r; params[i] = p+steps[i]; }
  if (chi2 < chi2_new) { chi2_new = chi2; params[i] = p; }
  steps[i]*=STEP_RATIO;
  //  fprintf(stderr, ".");
  return chi2_new;
}


void fitting_round(double *params, double *steps)
{
  float last_chi2, chi2;
  double cur_steps[NUM_PARAMS];
  int i, j=0;
  
  memcpy(cur_steps, steps, sizeof(double)*NUM_PARAMS);
  chi2 = calc_chi2(params);
  last_chi2 = chi2*(2+CHI2_LIMIT);
  while (chi2*(1+CHI2_LIMIT) < last_chi2 || (j++ < MIN_FIT_LIMIT)) {
    last_chi2 = chi2;
    for (i=0; i<NUM_PARAMS; i++)
      chi2 = improve_fit(params, cur_steps, i, chi2);
  }
}
