#ifndef CALC_SFH
#define CALC_SFH

#include "all_smf.h"

extern struct timestep *steps;
extern int64_t num_outputs;

void calc_sfh(struct smf_fit *f);
struct smf smhm_at_z(double z, struct smf_fit f);
// double calc_sm_at_m(double m, struct smf c);
double calc_sm_at_m(double m, struct timestep steps);
void calc_sm(int n, struct smf_fit *fit);
// double calc_smf_at_m(double m, double scale, struct smf c, struct smf_fit *fit);
void calc_smf_and_ssfr(int n, struct smf_fit *fit);
void calc_uvlf(int n, struct smf_fit *fit);
void create_fake_sfr_hist(int n, int i);
void calc_sm_hist(int n, struct smf_fit *fit);
//void calc_old_sm(int n);
// inline void calc_old_sm(int n, int j);
void calc_old_sm(int n, int j);
void calc_new_sm_and_sfr(int n, int i, struct smf_fit *fit);
  //void calc_new_sm_and_sfr(int n, struct smf_fit *fit);
void calc_total_sfr(int n);
void calc_total_bhar(int n);
void calc_observed_bhar(int n);
double calc_smf_at_sm(int64_t n, double sm);
double calc_uvlf_at_uv(int64_t n, double uv);
//double calc_smf_at_sm(int64_t n, int64_t j, double sm);
double recent_sfh_in_massive_halos(void);
double recent_sfh_in_massive_halos_nocorr(void);
double recent_Micl_Mstar_ratio_in_massive_halos(void);
double recent_radiative_power_in_massive_halos(void);
double recent_kinetic_power_in_massive_halos(void);
double recent_kinetic_frac_in_massive_halos(void);
double number_high_z_low_mass_qso(double z_low, double z_high, double lbol, double frac_area);
double ratio_high_z_qso(double z_low, double z_high, double lbol, double frac_area);
void calc_avg_eta_rad(int n);
double rising_sfh_penalty(void);

void calc_bh_integral(struct smf *c);
void calc_bh_eta_avg(int n);
void calc_bh_eta_kin_avg(int n);
void calc_bh_acc_rate_distribution(int n, struct smf_fit *fit);
void calc_bh_acc_rate_distribution_full(int n, struct smf_fit *fit);
void calc_bh_lum_distribution_full(int n, struct smf_fit *fit);
void calc_supEdd_frac_lum(int n, int i, double l);
double calc_bh_at_bm(double bm, struct smf c);
void calc_active_bh_fraction(int n, struct smf_fit *fit);
void calc_active_bh_fraction_lim(int n, struct smf_fit *fit, int ledd_or_lum, double lim);
#endif /* CALC_SFH */
