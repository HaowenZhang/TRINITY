#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_interp.h>
#include <assert.h>
#include <inttypes.h>
#include "all_luminosities.h"
#include "mlist.h"
#include "calc_sfh.h"
#include "check_syscalls.h"
#include "distance.h"
#include "universe_time.h"
#include "stringparse.h"

#define MIN_AGE 0

gsl_interp *lum_interps[LUM_TYPES][L_NUM_REDSHIFTS][L_NUM_ZS];
gsl_interp_accel lum_accel[LUM_TYPES][L_NUM_REDSHIFTS][L_NUM_ZS];
struct luminosities tlum = {0};
extern int64_t num_outputs;
char *lum_names[LUM_TYPES] = {"sloan_g",
			      "sloan_r", "sloan_i","cousins_i","f160w","f435w", "f606w", "f775w", "f814w", "f850lp", "2MASS_J", "2MASS_H", "2MASS_Ks", "fors_r",
			      "M1500_UV", "JWST_F200W", "JWST_F277W", "JWST_F356W", "JWST_F444W", "f105w", "f125w"};

double _lum_at_hm_z(float *lums, int64_t i, float z) {
  int64_t j;
  double f;
  calc_step_at_z(z, &j, &f);
  double l1 = lums[i*num_outputs + j];
  /*  if (fabs(z-5.95)<0.01) {
    fprintf(stderr, "L: %f, i: %"PRId64", j: %"PRId64", a: %f\n",
	    log10(l1)*-2.5, i, j, steps[j].scale);
	    }*/
  if (l1<=0) l1 = 1000;
  else l1 = log10(l1)*-2.5;
  if (j==num_outputs-1) return l1;
  double l2 = lums[i*num_outputs + j+1];
  if (l2<=0) return l1;
  else l2 = log10(l2)*-2.5;
  return (l1 + f*(l2-l1));
}

double lum_at_hm_z(float *lums, float hm, float z) {
  double mf, scale = 1.0/(1.0+z);
  mf = (hm-M_MIN)*BPDEX;
  int64_t i = mf;
  mf -= i;
  //  fprintf(stderr, "%f %f %"PRId64" %f\n", hm, z, i, scale);
  if (i<0) return 1000;
  if (i>=M_BINS-1) return 1000;
  if (scale < steps[0].scale) return 1000;
  if (scale > steps[num_outputs-1].scale) z = 1.0/steps[num_outputs-1].scale-1.0;
  double l1 = _lum_at_hm_z(lums, i, z);
  double l2 = _lum_at_hm_z(lums, i+1, z);
  return (l1 + mf*(l2-l1));
}

void gen_all_luminosities(enum lum_types t, float **lums, float z_src, float z_obs) {
  int64_t i,j;
  *lums = check_realloc(*lums, sizeof(float)*num_outputs*M_BINS, "Luminosities");
  for (i=0; i<num_outputs; i++)
    for (j=0; j<M_BINS; j++)
      lums[0][j*num_outputs + i] = luminosity_of_sfh(i, j, t, z_src, z_obs, 0);
}

double distance_factor(double z_src, double z_obs) {
  double com_dist = (comoving_distance(z_src) - comoving_distance(z_obs))/(1.0+z_obs);
  com_dist /= 1e-5; //10 pc in Mpc
  double z_eff = (1.0+z_src)/(1.0+z_obs) - 1.0;
  double d_fac = com_dist*com_dist*(1.0+z_eff);
  //Note the following factors which affect the luminosity:
  //1) Time dilation: 1.0/(1.0+z_eff)
  //2) Frequency redshifting: 1.0/(1.0+z_eff) (from loss of energy)
  //3) Frequency interval compression: (1.0+z_eff)
  //Total: 1.0/(1.0+z_eff)
  if (d_fac == 0) d_fac = 1; //Becomes absolute magnitude
  return 1.0/d_fac;
}

double luminosity_of_sfh(int64_t n, int64_t m_bin, enum lum_types t, float z_src, float z_obs, float dz_rsd) {
  int64_t i;
  double total_sm = 0;
  double total_l = 0;
  double age = 0;
  if (z_src < 0) z_src = 1.0/steps[n].scale-1.0;
  age = scale_to_years(1.0/(1.0+z_src)) - scale_to_years(0);
  if (z_src < z_obs) return 0;
  double d_fac = distance_factor(z_src, z_obs);
  z_src += dz_rsd;  //Redshift-space distortion
  if (z_src < z_obs) return 0;
  z_src = (1.0+z_src)/(1.0+z_obs) - 1.0;
  /*if (n==8 && m_bin==18) {
    fprintf(stderr, "z_src: %f; age: %e; d_fac: %e\n", z_src, age, d_fac);
    }*/
  for (i=0; i<=n; i++) {
    double new_sm = steps[n].sm_hist[m_bin*num_outputs + i];
    //Use instantaneous recycling fraction of 30%
    total_sm += new_sm*0.7;
    if (total_sm < 1) total_sm = 1;
    float z = 1.0/steps[i].scale-1.0;
    float Z = metallicity(log10(total_sm), z);
    float avg_lum = lum_at_age_Z(&tlum, t, age-steps[i].dt, age, Z, z_src);
    total_l += avg_lum*new_sm;
    /*if (n==8 && m_bin==22) {
      fprintf(stderr, "A1: %e; A2: %e; DA: %e; Z: %f; z_src: %f; Avg. L: %f; New_Sm: %e; Total_l: %f; Total_Sm: %e / %e\n", age-steps[i].dt, age, steps[i].dt, Z, z_src, log10(avg_lum*d_fac)*-2.5, new_sm, log10(total_l*d_fac)*-2.5, total_sm, steps[n].sm[m_bin]);
      }*/
    age -= steps[i].dt;
  }
  return total_l*d_fac;
}



void _luminosities_init(struct luminosities *lum) {
  int64_t i,j,k;
  for (i=0; i<LUM_TYPES; i++) {
    for (j=0; j<L_NUM_REDSHIFTS; j++) {
      for (k=0; k<L_NUM_ZS; k++) {
	lum_interps[i][j][k] = 
	  gsl_interp_alloc(gsl_interp_akima, L_NUM_AGES);
	memset(&(lum_accel[i][j][k]), 0, sizeof(gsl_interp_accel));
	gsl_interp_init(lum_interps[i][j][k], lum->ages, lum->lum[i][j][k], L_NUM_AGES);
      }
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
  double sm95 = 9.5 - M0;
  if (sm < sm95) {
    //Z = -0.0864*(sm95*sm95) + K0 + 0.3*(sm-sm95); //Kirby+ 2013: http://arxiv.org/abs/1310.0814
  }
  double Zsol = 8.69; //Asplund+ 2009
  return Z-Zsol;
}

double metallicity_iron(float sm, float z) {
  return (-0.09 + metallicity(sm, z)*1.18);
}

//From Mannucci et al. 2010: http://arxiv.org/pdf/1005.0006v3.pdf
double metallicity_mannucci(float sm, float sfr) {
  sm -= 10.0;
  double sfr_thresh = (-0.14+0.12*sm)/(2.0*0.054);
  if (sfr < sfr_thresh) sfr = sfr_thresh;
  double Z = 8.90 + 0.37*sm - 0.14*sfr - 0.19*sm*sm
    +0.12*sm*sfr - 0.054*sfr*sfr;
  double Zsol = 8.69; //Asplund+ 2009
  return Z-Zsol;
}

double metallicity_moustakas(float sm, float z) {
  double z0 = 9.115;
  double m0 = 9.0+2.043;
  double gamma = 1.41;
  double Z_01 = z0 - log10(1+pow(10, gamma*(sm-m0)));
  double ddz = -0.16 + 0.058*(sm-10.5);
  double Z = Z_01 + ddz*(z-0.1);
  double Zsol = 8.69;
  return Z-Zsol;
}


double __lum_at_age_Z(struct luminosities *lum, enum lum_types t, float age1, float age2, int64_t Z_bin, int64_t redshift_bin) {
  if (age1==age2) return gsl_interp_eval(lum_interps[t][redshift_bin][Z_bin], lum->ages, lum->lum[t][redshift_bin][Z_bin], age1, &(lum_accel[t][redshift_bin][Z_bin]));
  return (gsl_interp_eval_integ(lum_interps[t][redshift_bin][Z_bin], lum->ages, lum->lum[t][redshift_bin][Z_bin], age1, age2, &(lum_accel[t][redshift_bin][Z_bin]))/(age2-age1));
}


double _lum_at_age_Z(struct luminosities *lum, enum lum_types t, float age1, float age2, int64_t Z_bin, float z_src) {
  double f = z_src * 20;
  int64_t rbin = f;
  if (f <= 0) return __lum_at_age_Z(lum, t, age1, age2, Z_bin, 0);
  if (rbin >= L_NUM_REDSHIFTS-1) 
    return __lum_at_age_Z(lum, t, age1, age2, Z_bin, L_NUM_REDSHIFTS-1);
  f -= rbin;
  double l1 = __lum_at_age_Z(lum, t, age1, age2, Z_bin, rbin);
  double l2 = __lum_at_age_Z(lum, t, age1, age2, Z_bin, rbin+1);
  return (l1 + f*(l2-l1));
}

 double lum_at_age_Z(struct luminosities *lum, enum lum_types t, float age1, float age2, float Z, float z_src) {
  if (age1 <= lum->ages[0]) age1 = lum->ages[0];
  if (age1 < MIN_AGE) age1 = MIN_AGE;
  if (age2 <= age1) age2 = age1;
  if (age2 > lum->ages[L_NUM_AGES-1]) age2 = lum->ages[L_NUM_AGES-1];
  if (age2 <= age1) age1 = age2;
  int64_t Z_bin = Z/0.1+20;
  if (Z_bin < 0) Z_bin=0;
  if (Z_bin >= L_NUM_ZS) Z_bin = L_NUM_ZS - 1;
  while (Z_bin < L_NUM_ZS-1 && lum->Zs[Z_bin] < Z) Z_bin++;
  while (Z_bin > 0 && lum->Zs[Z_bin] > Z) Z_bin--;
  if ((lum->Zs[Z_bin] > Z) || (Z_bin==L_NUM_ZS-1)) 
    return _lum_at_age_Z(lum, t, age1, age2, Z_bin, z_src);
  
  double f = (Z-lum->Zs[Z_bin])/(lum->Zs[Z_bin+1]-lum->Zs[Z_bin]);
  double l1 = _lum_at_age_Z(lum, t, age1, age2, Z_bin, z_src);
  double l2 = _lum_at_age_Z(lum, t, age1, age2, Z_bin+1, z_src);
  return (l1 + f*(l2-l1));
}


void load_luminosities(char *filename) {
  char buffer[1024];
  double age, redshift, Z, bands[LUM_TYPES];
  int64_t n;
  int64_t num_ages = 0;
  int64_t num_Zs = 0;
  FILE *input, *output;
  SHORT_PARSETYPE;
  void *data[LUM_TYPES+3] = {&redshift, &age, &Z};
  enum parsetype types[LUM_TYPES+3] = {F64, F64, F64};
  for (n=0; n<LUM_TYPES; n++) {
    data[n+3] = bands+n;
    types[n+3] = F64;
  }

  sprintf(buffer, "%s.bin", filename);
  input = fopen(buffer, "r");
  if (input) {
    if (fread(&tlum, sizeof(tlum), 1, input)<1) {
      fclose(input);
      input = 0;
    } else {
      fclose(input);
    }
  }

  if (!input) {
    output = check_fopen(buffer, "w");
    input = check_fopen(filename, "r");
    for (n=0; n<L_NUM_REDSHIFTS; n++) tlum.redshifts[n] = n*0.05;

    while (fgets(buffer, 1024, input)) {
      if (buffer[0] == '#') continue;
      n = stringparse(buffer, data, types, 3+LUM_TYPES);
      if (n != 3+LUM_TYPES) {
	fprintf(stderr, "Syntax error in luminosities input!\n%s", buffer);
	exit(1);
      }
      int64_t Z_bin = 0;
      int64_t age_bin = 0;
      int64_t redshift_bin = 0;
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
      
      redshift_bin = redshift*20+0.5; //Positive values only
      assert(fabs(redshift - tlum.redshifts[redshift_bin]) < 0.001);

      for (n=0; n<LUM_TYPES; n++)
	tlum.lum[n][redshift_bin][Z_bin][age_bin] = pow(10, bands[n]/-2.5);
    }
    fclose(input);

    for (n=0; n<L_NUM_AGES; n++) tlum.ages[n] = pow(10, tlum.ages[n]);

    check_fwrite(&tlum, sizeof(tlum), 1, output);
    fclose(output);
  }
  _luminosities_init(&tlum);
  /* float test = log10(lum_at_age_Z(&tlum, L_sloan_r, 1e10, 1e10, 0, 0))*-2.5;
  //assert(fabs(test-6.811000) < 0.05); //AB
  float test2 = log10(lum_at_age_Z(&tlum, L_sloan_r, 1e9, 2e9, 0, 0))*-2.5;
  //assert(fabs(test2-4.908) < 0.05);
  float test3 = log10(lum_at_age_Z(&tlum, L_sloan_r, 1612410132.47544, 1612410132.47544, 0, 1))*-2.5;
  float dist = log10(distance_factor(1, 0))*-2.5;
  //fprintf(stderr, "%f\n", test3+dist);
  //assert(fabs(test3+dist - 51.43) < 0.05);
  test3 = log10(lum_at_age_Z(&tlum, L_f160w, 0, 5e7, 0, 15))*-2.5;
  fprintf(stderr, "%f\n", test3);
  test3 = log10(lum_at_age_Z(&tlum, L_jwst_f200w, 0, 5e7, 0, 15))*-2.5;
  fprintf(stderr, "%f\n", test3);
  test3 = log10(lum_at_age_Z(&tlum, L_jwst_f444w, 0, 5e7, 0, 15))*-2.5;
  fprintf(stderr, "%f\n", test3);
  */
}
