#ifndef SMF_H
#define SMF_H

#ifndef __APPLE__
#ifndef isfinite
#define isfinite(x) finitef(x)
#endif
#endif

#define MCMC_LENGTH ((1<<22)+100)
#define KAPPA_PRIOR 0.24
#define MU_PRIOR 0.14
#define SCATTER_PRIOR 0.03
#define SCATTER_CENTER 0.20

#define Z_COMPLETENESS_LIMIT 0.8

#define KAPPA_LIMIT 0.6
#define GAMMA_UPPER_LIMIT 50
#define DELTA_UPPER_LIMIT 5

#define BH_MASS_TO_REQUIRE_GROWTH 5
#define BH_MASS_TO_REQUIRE_ND_GROWTH 10.3 //A certain BH mass that we require to hava an increase in its number density from z=1 to z=0.

#define BH_GAMMA_CENTER 1.16
#define BH_GAMMA_CENTER_VD 4.377
#define BH_GAMMA_PRIOR 0.08
#define BH_BETA_CENTER 8.50
#define BH_BETA_CENTER_VD 8.576
//#define BH_BETA_PRIOR 0.05
#define BH_BETA_PRIOR 0.2
#define BH_SCATTER_CENTER 0.4
//#define BH_SCATTER_PRIOR 0.03
#define BH_SCATTER_PRIOR 0.1

#define SIGMA_CENTER 0.07
#define SIGMA_Z_CENTER 0.05
#define SIGMA_Z_PRIOR 0.015
#define SIGMA_Z_START 0.2

#define SM_0_LIMIT 14

#ifndef M_PI
#define M_PI        3.14159265358979323846264338327950288
#endif

struct smf {
  double epsilon, sm_0, v_1, alpha, delta, beta, gamma, lambda, scatter, obs_scatter, 
  				mu, kappa, sm_completeness, csfr_completeness, 
  				ssfr_corr, sfr_sm_corr,sm_max,sm_min, uv_max, uv_min, passive_mass,f_1,lm_slope,mpk, 
  				combined_scatter, 
  				// icl_m, 
  				icl_frac,
  				scatter_corr, valid, qm, qwidth, fqmin,
  				bh_beta, bh_gamma, bh_merge, bh_merge_width, bh_alpha, bh_delta,
  				bh_efficiency, bh_efficiency_rad, bh_duty, bh_scatter, bh_prob_norm, 
  				dc_mbh, dc_mbh_w, eta_mu, rho_bh, bh_eta_crit, abhmf_shift, f_merge_bh, bh_duty_m, bh_duty_alpha,
          f_occ_min, log_bh_scatter_corr;
 int bin_real;
};

#ifndef exp10f
#define exp10f(x) expf((float)M_LN10*(x))
#endif

#ifndef exp10
#define exp10(x) exp(M_LN10*(x))
#endif

#endif /* SMF_H */
