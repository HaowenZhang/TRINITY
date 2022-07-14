#ifndef OBS_SMF_H
#define OBS_SMF_H

#include "smf.h"
#include "all_smf.h"
#include "mf_cache.h"
#include <inttypes.h>

#define BOUWENS_SFR_LIMIT (0.7/1.7)

#define SLOPE_SYSTEMATIC_MASS 11.3
//#define PHI_INTEGRAL_PRECISION 1e-4
#define PHI_INTEGRAL_PRECISION 1e-5
#define MAX_SMF_POINTS 200

#define CALCULATION_TOLERANCE 0.01
#define FIT_TOLERANCE 0.02

#define GAUSS_CACHE_SIZE 4097
//#define GAUSS_CACHE_MIN -2
#define GAUSS_CACHE_MIN -8
#define GAUSS_CACHE_MAX 8
#define UV_MAG_OFFSET 6
// #define GAUSS_CACHE_MIN -1
// #define GAUSS_CACHE_MAX 1
#define GAUSS_CACHE_CENTER ((GAUSS_CACHE_SIZE-1)/2)
#define GAUSS_CACHE_INV_STEP ((((double)GAUSS_CACHE_SIZE-1.0) / ((double)(GAUSS_CACHE_MAX-GAUSS_CACHE_MIN))))
#define GAUSS_CACHE_INV_STEPF ((double)GAUSS_CACHE_INV_STEP)
#define GAUSS_CACHE_STEP (1.0/GAUSS_CACHE_INV_STEP)

#define SMF_TYPE 0
#define SSFR_TYPE 1
#define CSFR_TYPE 2
#define QLF_TYPE 3
#define QPDF_TYPE 4
#define ABHMF_TYPE 5
#define QPDF_ETA_TYPE 6
#define QLF_CTN_TYPE 7
#define QF_TYPE 8
#define UVLF_TYPE 9
#define BHMF_TYPEI_TYPE 10

#define SC_ALPHA_MIN -2
#define SC_ALPHA_MAX 2
#define SC_BPUNIT 100
#define SC_BINS ((int64_t)((SC_ALPHA_MAX - SC_ALPHA_MIN)*SC_BPUNIT + 1)


struct real_smf 
{
  double mass, val, err_l, err_h;
};

struct obs_smf 
{
  struct real_smf real_smf[MAX_SMF_POINTS];
  int real_smf_points;
  double z_low, z_high;
  double chi2;
  int linear_errors;
  int type;
};

struct obs_smf_point 
{
  double mass, extra, val, err_l, err_h;
  double z_low, z_high, v_low, v_high;
  double chi2;
  int linear_errors;
  int type;
};

struct step_integral_helper_data 
{
  double f, sm, scatter, obs_scatter, gauss_inv_sigma,corr,s_corr,kappa;
  double passive_mass;
  int64_t step, n, n2;
  double *list, *list2;
  double *sfrac_sm, *sfrac_sm2;
  double *sfr_sm_sf, *sfr_sm_q, *sfr_sm_sf2, *sfr_sm_q2;
  char *ok, *ok2;
};

struct quasar_info 
{
  double l;
  int64_t step;
};


double mf_cache(double scale, double mass);
double evaluate_from_step(double z, double sm, int smf);
double calc_single_chi2_err(struct obs_smf *cur_smf);
// inline double model_chi2(struct obs_smf *m_smf, struct obs_smf *r_smf);
double model_chi2(struct obs_smf *m_smf, struct obs_smf *r_smf);
double chi2_err_helper(double Vc, void *extra_data);
double chi2_err_helper_qf(double Vc, void *extra_data);
double chi2_err_helper_uv(double Vc, void *extra_data);
double distance_weight_helper(double z, void *extra_data);

void load_real_smf(struct obs_smf_point **cur_smf, int64_t *num_points, char *filename);
void setup_psf(int scatter);
void load_mf_cache(char *filename);
extern struct z_mf_cache *all_mf_caches;

double adaptiveGauss(double (*f)(double, void*),   // ptr to function
                     void *extra_data,
                     double a, double b,  // interval [a,b]
                     double rel_err,  // error tolerance
		     int integral_level);
double smf_integral_helper(double sm, void *extra_data);
void setup_psf_cache(double obs_scatter, double scatter);
// inline 
double syst_cache(double dm);
void init_gauss_cache(void);
double find_Lx_at_Lbol(double Lbol);
double find_Lbol_at_Lx(double Lx);
double calc_ssfr(double sm, double z);
double calc_cosmic_sfr(double z);
double calc_single_chi2_err_point(struct obs_smf_point *cur_smf);
double schechter_norm(double alpha, double lower, double upper, struct smf_fit *fit);
double doublePL_norm(double alpha, double delta, double lower, double upper, struct smf_fit *fit);
double schechter_inv_avg(double alpha, double lower, double upper, struct smf_fit *fit);
double schechter_frac_above_thresh(double thresh, double alpha, double norm, struct smf_fit *fit);
double doublePL_frac_above_thresh(double thresh, double alpha, double delta, double norm, struct smf_fit *fit);
double quasar_lf_helper(double m, void *extra_info);
double calc_quasar_lf(double l, double z);
double calc_quasar_lf_new(double Mi, double z);
double calc_quasar_lf_ctn(double l, double z);
double calc_quasar_lf_split(double l, double z, double Mh_min, double Mh_max);
double calc_quasar_lf_split_sm(double l, double z, double SM_min, double SM_max);
double bulge_mass(double sm, double a);
double calc_quasar_lf_mbh(double Mi, double z, double mbh_low, double mbh_high);
double calc_quasar_lf_eta(double Mi, double z, double eta_low, double eta_high);

double calc_bhmf(double m, double z);
double calc_bhmf_mh(double m, double z, double mh_low, double mh_high);
double calc_bhmf_typeI(double m, double z);

double calc_bhar_mbh(double m, double z);
double calc_bhar_mstar(double m, double z);
double calc_bher_mbh(double m, double z);
double calc_bher_mstar(double m, double z);
double calc_bhmr_mbh(double m, double z);
double calc_bhmr_mstar(double m, double z);
double calc_bhar_sfr_mbh(double m, double z);
double calc_bhar_sfr_mstar(double m, double z);
double calc_sbhar_ssfr_mbh(double m, double z);
double calc_sbhar_ssfr_mstar(double m, double z);
double calc_sbhar_mstar(double m, double z);
double calc_sfr_mstar(double m, double z);

double calc_qpdf_at_l_m_z(double l, double m, double z);
double calc_qpdf_at_eta_m_z(double eta, double m, double z);
double calc_qpdf_at_sBHAR_m_z(double sBHAR, double m, double z);
double calc_qpdf_at_sBHAR_m_z_new(double sBHAR, double m, double z);

double calc_active_bhmf(double m, double z);
double cosmic_bh_density(double z, double thresh_ledd, double thresh_lbol, struct smf_fit *fit);
double cosmic_bh_density_split(double z, double Mh_low, double Mh_high, struct smf_fit *fit);
double cosmic_old_bh_density(double z, double thresh_ledd, double thresh_lbol, struct smf_fit *fit);
double cosmic_unmerged_bh_density(double z, double thresh_ledd, double thresh_lbol, struct smf_fit *fit);

double _prob_of_ledd_nonlinear(double ledd, double bher_char, double alpha, double delta, double bh_prob_norm, double bh_eta_crit);
double _prob_of_ledd_kinetic(double ledd, double bher_char, double alpha, double delta, double bh_prob_norm, double bh_eta_crit);
double _prob_of_ledd_linear(double ledd, double bher_char, double alpha, double delta, double bh_prob_norm, double bh_eta_crit);


extern double gauss_inv_sigma;
extern double _gauss_inv_sigma;
extern double gauss_sigma;

#endif /* OBS_SMF_H */
