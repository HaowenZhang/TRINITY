#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_interp.h>
#include <assert.h>
#include <inttypes.h>
#include "luminosities.h"
#include "mlist.h"
#include "calc_sfh.h"
#include "check_syscalls.h"

gsl_interp *lum_interps[LUM_TYPES][L_NUM_ZS];
gsl_interp_accel lum_accel[LUM_TYPES][L_NUM_ZS];
struct luminosities tlum = {0};

void gen_all_luminosities(enum lum_types t, float **lums) {
  int64_t i,j;
  *lums = check_realloc(*lums, sizeof(float)*num_outputs*M_BINS, "Luminosities");
  for (i=0; i<num_outputs; i++)
    for (j=0; j<M_BINS; j++)
      lums[0][j*num_outputs + i] = luminosity_of_sfh(i, j, t);
}

double luminosity_of_sfh(int64_t n, int64_t m_bin, enum lum_types t) {
  int64_t i;
  double total_sm = 0;
  double total_l = 0;
  double age = 0;
  for (i=num_outputs-1; i>=0; i--) age += steps[i].dt;
  for (i=0; i<=n; i++) {
    double new_sm = steps[n].sm_hist[m_bin*num_outputs + i];
    //Use instantaneous recycling fraction of 30%
    total_sm += new_sm*0.7;
    if (total_sm < 1) total_sm = 1;
    float z = 1.0/steps[i].scale-1.0;
    float Z = metallicity(log10(total_sm), z);
    total_l += lum_at_age_Z(&tlum, t, age-steps[i].dt, age, Z)*new_sm;
    age -= steps[i].dt;
  }
  return total_l;
}



void _luminosities_init(struct luminosities *lum) {
  int64_t i,j;
  for (i=0; i<LUM_TYPES; i++) {
    for (j=0; j<L_NUM_ZS; j++) {
      lum_interps[i][j] = 
	gsl_interp_alloc(gsl_interp_akima, L_NUM_AGES);
      memset(&(lum_accel[i][j]), 0, sizeof(gsl_interp_accel));
      gsl_interp_init(lum_interps[i][j], lum->ages, lum->lum[i][j], L_NUM_AGES);
    }
  }
}

//From Maiolino et al. 2008: http://arxiv.org/abs/0806.2410v2
double metallicity(float sm, float z) {
  float M0 = 11.219+0.466*z;
  float K0 = 9.07-0.0701*z;
  if (sm > M0) sm = M0;
  sm -= M0;
  double Z = -0.0864*(sm*sm) + K0; //12 + log(O/H)
  double Zsol = 8.69; //Asplund+ 2009
  return Z-Zsol;
}

double _lum_at_age_Z(struct luminosities *lum, enum lum_types t, float age1, float age2, int64_t Z_bin) {
  if (age1==age2) return gsl_interp_eval(lum_interps[t][Z_bin], lum->ages, lum->lum[t][Z_bin], age1, &(lum_accel[t][Z_bin]));
  return (gsl_interp_eval_integ(lum_interps[t][Z_bin], lum->ages, lum->lum[t][Z_bin], age1, age2, &(lum_accel[t][Z_bin]))/(age2-age1));
}

double lum_at_age_Z(struct luminosities *lum, enum lum_types t, float age1, float age2, float Z) {
  if (age1 <= lum->ages[0]) age1 = lum->ages[0];
  if (age2 <= age1) age2 = age1;
  if (age2 > lum->ages[L_NUM_AGES-1]) age2 = lum->ages[L_NUM_AGES-1];
  if (age2 <= age1) age1 = age2;
  int64_t Z_bin = Z/0.1+20;
  if (Z_bin < 0) Z_bin=0;
  if (Z_bin >= L_NUM_ZS) Z_bin = L_NUM_ZS - 1;
  while (Z_bin < L_NUM_ZS-1 && lum->Zs[Z_bin] < Z) Z_bin++;
  while (Z_bin > 0 && lum->Zs[Z_bin] > Z) Z_bin--;
  if ((lum->Zs[Z_bin] > Z) || (Z_bin==L_NUM_ZS-1)) 
    return _lum_at_age_Z(lum, t, age1, age2, Z_bin);
  
  double f = (Z-lum->Zs[Z_bin])/(lum->Zs[Z_bin+1]-lum->Zs[Z_bin]);
  double l1 = _lum_at_age_Z(lum, t, age1, age2, Z_bin);
  double l2 = _lum_at_age_Z(lum, t, age1, age2, Z_bin+1);
  return (l1 + f*(l2-l1));
}


void load_luminosities(char *filename) {
  char buffer[1024];
  double age, Z, bands[LUM_TYPES];
  int64_t n;
  int64_t num_ages = 0;
  int64_t num_Zs = 0;
  FILE *input;
  assert(LUM_TYPES==5);

  input = check_fopen(filename, "r");
  while (fgets(buffer, 1024, input)) {
    if (buffer[0] == '#') continue;
    n = sscanf(buffer, "%lf %lf %lf %lf %lf %lf %lf\n", &age, &Z, bands,
	       bands+1, bands+2, bands+3, bands+4);
    if (n != 2+LUM_TYPES) {
      fprintf(stderr, "Syntax error in luminosities input!\n%s", buffer);
      exit(1);
    }
    int64_t Z_bin = 0;
    int64_t age_bin = 0;
    for (age_bin=0; age_bin<num_ages; age_bin++)
      if (tlum.ages[age_bin] == age) break;
    if (age_bin==num_ages) {
      assert(num_ages < L_NUM_AGES);
      num_ages++;
      tlum.ages[age_bin] = age;
      if (age_bin) assert(tlum.ages[age_bin] > tlum.ages[age_bin-1]);
    }

    for (Z_bin=0; Z_bin<num_Zs; Z_bin++)
      if (tlum.Zs[Z_bin] == Z) break;
    if (Z_bin==num_Zs) {
      assert(num_Zs < L_NUM_ZS);
      num_Zs++;
      tlum.Zs[Z_bin] = Z;
      if (Z_bin) assert(tlum.Zs[Z_bin] > tlum.Zs[Z_bin-1]);
    }

    for (n=0; n<LUM_TYPES; n++)
      tlum.lum[n][Z_bin][age_bin] = pow(10, bands[n]/-2.5);
  }
  fclose(input);

  for (n=0; n<L_NUM_AGES; n++) tlum.ages[n] = pow(10, tlum.ages[n]);
  tlum.redshift = 0;

  _luminosities_init(&tlum);
  float test = log10(lum_at_age_Z(&tlum, L_sloan_r, 1e10, 1e10, 0))*-2.5;
  assert(fabs(test-6.436000) < 0.05);
  float test2 = log10(lum_at_age_Z(&tlum, L_sloan_r, 1e9, 2e9, 0))*-2.5;
  assert(fabs(test2-4.550692) < 0.05);
}
