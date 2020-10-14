#ifndef MLIST_H
#define MLIST_H

#include <stdint.h>
#include <gsl/gsl_spline.h>
#include "smf.h"

#define M_MIN 7
#define M_MAX 16

/*#define M_MIN 9
  #define M_MAX 16*/
#define BPDEX 5
#define M_BINS ((int64_t)((M_MAX-M_MIN)*BPDEX + 1))
#define INV_BPDEX (1.0/((double)BPDEX))

#define SM_MIN 5
#define SM_MAX 14
#define SM_BPDEX 100
#define SM_BINS ((int64_t)((SM_MAX-SM_MIN)*SM_BPDEX + 1))
#define SM_INV_BPDEX (1.0/((double)SM_BPDEX))
#define SM_EXTRA 20

#define UV_MIN -25
#define UV_MAX -10
#define UV_BPMAG 50
#define UV_BINS ((int64_t)((UV_MAX-UV_MIN)*UV_BPMAG + 1))
#define UV_INV_BPMAG (1.0/((double)UV_BPMAG))
#define UV_EXTRA 20

// The MBH bin specs to build the bh_unmerged_dist table.
#define MBH_MIN 0
#define MBH_MAX 10
#define MBH_BPDEX 20
#define MBH_BINS 100
#define MBH_INV_BPDEX (1.0/((double)MBH_BPDEX))

#define LBOL_MIN 36
#define LBOL_MAX 50
#define LBOL_BPDEX 10
#define LBOL_INV_BPDEX (1.0 / ((double) LBOL_BPDEX))
#define LBOL_BINS ((int64_t)((LBOL_MAX - LBOL_MIN) * LBOL_BPDEX)) 

#define BHER_MIN -6
#define BHER_MAX 6
#define BHER_EFF_MAX 2
#define BHER_EFF_MIN -4
#define BHER_BPDEX 20
#define BHER_INV_BPDEX (1.0/((double)BHER_BPDEX))
#define BHER_BINS ((int64_t)((BHER_MAX-BHER_MIN)*BHER_BPDEX + 1))
//#define BHER_BINS ((int64_t)((BHER_MAX-BHER_MIN)*20 + 1))


struct timestep {
  float c[M_BINS];
  float merged[M_BINS*M_BINS];
  float mmp[M_BINS*M_BINS];
  float bher_dist_full[M_BINS*BHER_BINS];
  float bher_dist_kin[M_BINS*BHER_BINS];
  float n[M_BINS];
  float t[M_BINS];
  double icl_exp;
  double scale;
  double dt; //In yr
  double *smloss;
  double *sm_hist;
  double *sm_inc_hist; //record the history of the mass that come from the incoming satellite galaxies, 
                        //***NOT*** including ICL from satellites.
  double *icl_stars;
  //double icl_frac[M_BINS];
  double sm_icl[M_BINS];
  double mfrac[M_BINS];
  double new_sm[M_BINS];
  double old_sm[M_BINS];
  double sm[M_BINS];
  double sm_avg[M_BINS];
  double log_bm[M_BINS];
  double re[M_BINS];       //Effective radius
  double vdisp[M_BINS];    //Velocity dispersion
  double vdisp_los[M_BINS]; //Velocity dispersion in line-of-sight direction.
  double bher_dist[M_BINS*BHER_BINS]; //Eddington rate distribution
  double bher_dist_norm[M_BINS];
  double ledd_min[M_BINS];
  double ledd_eff_min[M_BINS];
  double ledd_max[M_BINS];
  double ledd_min_abs; //The minimum/maximum ABSOLUTE (NOT RELATIVE TO THE TYPICAL VALUE) 
  double ledd_max_abs; //RADIATIVE Eddington ratio among all the mass bins at this snapshot. 
  double bh_mass_min; //The minimum/maximum (log) bh mass corresponding to the Eddington 
  double bh_mass_max; //ratio limits defined above. bh_mass_min/max == l_min/max - ledd_max/min_abs
  
  double lum_dist_full[LBOL_BINS * MBH_BINS]; //The mass function/probability 
  double lum_func_full[LBOL_BINS * MBH_BINS]; //distribution of luminosities given Mbh.
  double bhmf[MBH_BINS];
  double ledd_bpdex[M_BINS];
  double lv[M_BINS];
  double sfrac[M_BINS];
  double sfr[M_BINS];
  double mr[M_BINS];
  double sm_nd[M_BINS];
  double sm_inc[M_BINS];
  double merged_frac[M_BINS];
  double frac_supEdd_l[M_BINS];
 
  double smf[SM_BINS+SM_EXTRA];
  double sfrac_sm[SM_BINS+SM_EXTRA];
  double m_at_sm[SM_BINS+SM_EXTRA];
  char smf_ok[SM_BINS+SM_EXTRA];
  double sfr_sm[SM_BINS+SM_EXTRA];
  double sfr_sm_sf[SM_BINS+SM_EXTRA];
  double sfr_sm_q[SM_BINS+SM_EXTRA];
  double cosmic_sfr, observed_cosmic_sfr;
  double cosmic_bhar, observed_cosmic_bhar;
  double bh_growth[M_BINS];
  double log_sm[M_BINS];
  double log_sm_obs[M_BINS];
  double med_hm_at_a[M_BINS];

  double k_uv[M_BINS];
  double b_uv[M_BINS];
  double k_std_uv[M_BINS];
  double b_std_uv[M_BINS];
  double obs_uv[M_BINS];
  double std_uv[M_BINS];
  double std_uvlf[UV_BINS+UV_EXTRA];
  double uvlf[UV_BINS+UV_EXTRA];
  char uvlf_ok[UV_BINS+UV_EXTRA];

  double log_bh_mass[M_BINS];
  double bh_mass_avg[M_BINS];
  double old_bh_mass[M_BINS];
  double old_bh_mass_self[M_BINS];
  double old_bh_mass_next[M_BINS];
  double bh_eta_avg[M_BINS];
  double bh_eta_kin_avg[M_BINS];
  double bh_acc_rate[M_BINS];
  double bh_acc_rate_obs[M_BINS]; //BH accretion rates that induce a luminosity above the observable limit.
  double bh_merge_rate[M_BINS];
  double bh_eta[M_BINS];
  double dc_obs[M_BINS];
  double bh_unmerged[M_BINS];
  // double bh_merged[M_BINS];
  double new_bh_mass[M_BINS];
  double bh_unmerged_avail[M_BINS]; //The total available unmerged BH mass (i.e. bh_unmerged after subtracting merged masses)
  double bh_unmerged_dist[M_BINS*MBH_BINS]; //The total available unmerged BH mass (i.e. bh_unmerged after subtracting merged masses)

  double n_merge10[M_BINS]; //The number of mergers with BH mass ratio > 1:10
  double n_merge100[M_BINS]; //The number of mergers with BH mass ratio > 1:100
  double f_CTK[M_BINS]; //The Compton-thick fractions.
  //double bh_eta_corr; //The correction from the average eta to characteristic eta: average + bh_eta_corr = char.
  double bh_Lbol[M_BINS];
  double bh_sBHAR[M_BINS];
  double f_active[M_BINS]; //Active fraction
  gsl_spline *spline;
  gsl_spline *spline2; //Since the relation may not be monotonic, we have to split the SM-Mh relation into two pieces to do the interpolation.
  gsl_spline *spline_sfr;
  gsl_spline *spline_sfr2;
  gsl_spline *spline_uv;
  gsl_spline *spline_uv2; //Since the relation is not monotonic, we have to split the UV-Mh relation into two pieces to do the interpolation.
  // gsl_spline *spline_std_uv;

  int alloc2smf;
  int alloc2uvlf;

  int flag_alloc;
  double sm_from_icl[M_BINS];
  struct smf smhm;
};

void init_timesteps(void);
void calc_step_at_z(double z, int64_t *step, double *f);

#endif /* MLIST_H */
