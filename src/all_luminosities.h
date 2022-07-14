#ifdef _LUMINOSITIES_H_
#error "luminosities.h and all_luminosities.h cannot be loaded simultaneously!"
#endif

#ifndef _ALL_LUMINOSITIES_H_
#define _ALL_LUMINOSITIES_H_

#define L_NUM_AGES 188
#define L_NUM_ZS 22
#define L_NUM_REDSHIFTS 301
#define LUM_TYPES 21

enum lum_types {
  L_sloan_g=0,
  L_sloan_r,
  L_sloan_i,
  L_cousins_i,
  L_f160w,
  L_f435w,
  L_f606w,
  L_f775w,
  L_f814w,
  L_f850lp,
  L_2mass_j,
  L_2mass_h,
  L_2mass_k,
  L_fors_r,
  L_m1500,
  L_jwst_f200w,
  L_jwst_f277w,
  L_jwst_f356w,
  L_jwst_f444w,
  L_f105w,
  L_f125w,
};

extern char *lum_names[LUM_TYPES];

struct luminosities {
  double unused;
  double ages[L_NUM_AGES];
  double Zs[L_NUM_ZS];
  double redshifts[L_NUM_REDSHIFTS];
  double lum[LUM_TYPES][L_NUM_REDSHIFTS][L_NUM_ZS][L_NUM_AGES];
};

extern struct luminosities tlum;

void gen_all_luminosities(enum lum_types t, float **lums, float z_src, float z_obs);
double lum_at_hm_z(float *lums, float hm, float z);
double luminosity_of_sfh(int64_t n, int64_t m_bin, enum lum_types t, float z_src, float z_obs, float dz_rsd);
double metallicity(float sm, float z);
double metallicity_iron(float sm, float z);
double lum_at_age_Z(struct luminosities *lum, enum lum_types t, float age1, float age2, float Z, float z_src);
void load_luminosities(char *filename);
double metallicity_moustakas(float sm, float z);
double metallicity_mannucci(float sm, float sfr);
double distance_factor(double z_src, double z_obs);


#endif /* _ALL_LUMINOSITIES_H_ */
