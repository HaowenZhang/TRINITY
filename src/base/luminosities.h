#ifndef _LUMINOSITIES_H_
#define _LUMINOSITIES_H_

#define L_NUM_AGES 188
#define L_NUM_ZS 22
#define LUM_TYPES 5

enum lum_types {
  L_sloan_r=0,
  L_sloan_i,
  L_cousins_i,
  L_cousins_i2,
  L_f160w,
};

struct luminosities {
  double redshift;
  double ages[L_NUM_AGES];
  double Zs[L_NUM_ZS];
  double lum[LUM_TYPES][L_NUM_ZS][L_NUM_AGES];
};


void gen_all_luminosities(enum lum_types t, float **lums);
double luminosity_of_sfh(int64_t n, int64_t m_bin, enum lum_types t);
double metallicity(float sm, float z);
double lum_at_age_Z(struct luminosities *lum, enum lum_types t, float age1, float age2, float Z);
void load_luminosities(char *filename);


#endif /* _LUMINOSITIES_H_ */
