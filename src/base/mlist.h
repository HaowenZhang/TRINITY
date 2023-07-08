#ifndef MLIST_H
#define MLIST_H

#include <stdint.h>
#include <gsl/gsl_spline.h>
#include "smf.h"

// Minimum/maximum halo peak masses
#define M_MIN 7
#define M_MAX 16

/*#define M_MIN 9
  #define M_MAX 16*/
#define BPDEX 5 // # of bins per dex in halo peak mass.
#define M_BINS ((int64_t)((M_MAX-M_MIN)*BPDEX + 1)) // total # of halo mass bins.
#define INV_BPDEX (1.0/((double)BPDEX))

#define SM_MIN 5 //minimum/maximum stellar mass for intrinsic stellar mass functions.
#define SM_MAX 14
#define SM_BPDEX 100 // # of bins per dex in stellar mass
#define SM_BINS ((int64_t)((SM_MAX-SM_MIN)*SM_BPDEX + 1))
#define SM_INV_BPDEX (1.0/((double)SM_BPDEX))
#define SM_EXTRA 20

#define UV_MIN -25 // minimum/maximum UV magnitudes for observed UV luminosity functions.
#define UV_MAX -10
#define UV_BPMAG 50 // # of bins per dex in UV magnitude
#define UV_BINS ((int64_t)((UV_MAX-UV_MIN)*UV_BPMAG + 1))
#define UV_INV_BPMAG (1.0/((double)UV_BPMAG))
#define UV_EXTRA 20

#define MBH_MIN 2 // minimum/maximum BH masses for 2-dimensional M_BH--Lbol distributions.
#define MBH_MAX 13
#define MBH_BINS 100 // # of bins in BH mass
#define MBH_BPDEX ((double)(MBH_BINS) / (MBH_MAX - MBH_MIN))
#define MBH_INV_BPDEX (1.0/((double)MBH_BPDEX))

#define LBOL_MIN 36 // minimum/maximum AGN bolometric luminosities for 2-dimensional M_BH--Lbol distributions.
#define LBOL_MAX 52
#define LBOL_BPDEX 10 // # of bins in bolometric luminosity
#define LBOL_INV_BPDEX (1.0 / ((double) LBOL_BPDEX))
#define LBOL_BINS ((int64_t)((LBOL_MAX - LBOL_MIN) * LBOL_BPDEX)) 

#define BHER_MIN -6 // minimum/maximum ***fractional*** Eddington ratios.
#define BHER_MAX 6
#define BHER_EFF_MAX 2 // minimum/maximum ***fractional*** Eddington ratios at which
                       // the Eddington ratio distributions are calculated.
#define BHER_EFF_MIN -4 
#define BHER_BPDEX 20 // # of bins per dex in Eddington ratios
#define BHER_INV_BPDEX (1.0/((double)BHER_BPDEX))
#define BHER_BINS ((int64_t)((BHER_MAX-BHER_MIN)*BHER_BPDEX + 1))
//#define BHER_BINS ((int64_t)((BHER_MAX-BHER_MIN)*20 + 1))


struct timestep {
  double med_hm_at_a[M_BINS]; // The median halo mass for each halo mass bin.
  float c[M_BINS]; // the number density of halos that newly emerged 
                   // (instead of inherited from the previous snapshot.)
  float merged[M_BINS*M_BINS]; // The number density of the halos in a mass bin
                               // that are merged into another bin.
  float mmp[M_BINS*M_BINS]; // The number density of halos in a mass bin that are 
                            // merged into another mass bin.

  float bher_dist_full[M_BINS*BHER_BINS]; // The full radiative Eddington ratio 
                                          // distributions. Now deprecated.
  float bher_dist_kin[M_BINS*BHER_BINS]; // The kinetic Eddington ratio distributions.
                                         // In the MCMC process, this is not calculated
                                         // because no observational data are
                                          // used to provide constraints.
  float n[M_BINS]; // The number density of halos in a certain mass bin that is 
                   // inherited from last snapshot.
  float t[M_BINS]; // The total number density of halos in a certain mass bin.
  double icl_exp; // Deprecated.
  double scale; // The scale factor of the snapshot.
  double dt; // The time interval for each snapshot. In yrs.
  double *smloss; // The ***remaining*** fraction of stellar mass (not lost) due to stellar
                  // evolution.
  double *sm_hist; // The stellar mass history of all halo mass bins.
  double *sm_inc_hist; //record the history of the mass that come from the incoming satellite galaxies, 
                        //***NOT*** including ICL from satellites.
  double *icl_stars; // The intracluster light (ICL) mass history of all halo mass bins.

  double sm_icl[M_BINS]; // The stellar masses in ICL.
  double mfrac[M_BINS]; // The fractional merger contribution to total BH growth.
  double new_sm[M_BINS]; // New stellar mass compared to the old stellar mass.
  double old_sm[M_BINS]; // Old stellar mass.
  double sm[M_BINS]; // ***Median*** stellar mass.
  double sm_avg[M_BINS]; // ***Average*** stellar mass.
  double log_sm[M_BINS]; // the log10 of intrinsic median stellar mass.
  double log_sm_obs[M_BINS]; // the log10 of observed median stellar mass.
  double log_bm[M_BINS]; // the log10 of ***observed*** bulge mass.
  double re[M_BINS];       //Effective radius. Deprecated.
  double vdisp[M_BINS];    //Velocity dispersion. Deprecated.
  double vdisp_los[M_BINS]; //Velocity dispersion in line-of-sight direction. Deprecated.
  double bher_dist[M_BINS*BHER_BINS]; // Eddington ratio distribution.
  double bher_dist_norm[M_BINS]; // The normalization of Eddingron ratio distributions.
  double ledd_min[M_BINS]; // The minimum ***fractional*** Eddington ratios (relative
                           // to the typical Eddington ratio in the double power-law)
                           // This might be different for different halo mass bins, 
                           // due to the non-linear scaling relation between the total 
                           // and radiative Eddington ratio distributions.
  double ledd_eff_min[M_BINS]; // The effective minimum ***fractional*** Eddington ratios 
                               // (relative to the typical Eddington ratio in the double 
                               // power-law) This might be different for different halo 
                               // mass bins, due to the non-linear scaling relation 
                               // between the total and radiative Eddington ratio distributions.
  double ledd_max[M_BINS]; // The maximum ***fractional*** Eddington ratios (relative
                           // to the typical Eddington ratio in the double power-law)
                           // This might be different for different halo mass bins, 
                           // due to the non-linear scaling relation between the total 
                           // and radiative Eddington ratio distributions.
  double ledd_min_abs; //The minimum/maximum ABSOLUTE (NOT RELATIVE TO THE TYPICAL VALUE) 
  double ledd_max_abs; //RADIATIVE Eddington ratio among all the mass bins at this snapshot. 
  double bh_mass_min; //The minimum/maximum (log) bh mass.
  double bh_mass_max; 
  
  double lum_func_norm[MBH_BINS];
  double lum_dist_full[LBOL_BINS * MBH_BINS]; //The luminosity function/probability 
  double lum_func_full[LBOL_BINS * MBH_BINS]; //distribution of luminosities given Mbh.
                                              //The difference between the two is that
                                              //lum_dist_full only covers active BHs,
                                              //whereas lum_func_full is normalized to
                                              //match AGN duty cycles.
  double bhmf[MBH_BINS]; // Total BH mass function, including active+dormant BHs.
  double ledd_bpdex[M_BINS]; // The # of bins per dex in Eddington ratios.
  double lv[M_BINS]; // log10 of maximum circular velocities.
  double sfrac[M_BINS]; // The fraction of star-forming galaxies in each halo mass bin.
  double sfr[M_BINS]; // The average star formation rates in each halo mass bin.
  double mr[M_BINS]; // The average galaxy merger rates in each halo mass bin.
  // double sm_nd[M_BINS]; // Deprecated.
  double sm_inc[M_BINS]; // The incoming stellar mass from satellite galaxies for 
                         // each halo mass bin.
  double merged_frac[M_BINS]; // The fraction of incoming satellite stellar mass
                              // that got merged.
  double frac_supEdd_l[M_BINS]; // The fraction of super-Eddington objects among
                                // quasars above a certain luminosity limit.
 
  double smf[SM_BINS+SM_EXTRA]; // The intrinsic stellar mass function.
  double sfrac_sm[SM_BINS+SM_EXTRA]; // The star-forming fraction of galaxies,
                                     // as a function of stellar mass.
  double m_at_sm[SM_BINS+SM_EXTRA]; // Deprecated.
  char smf_ok[SM_BINS+SM_EXTRA]; // Flag that indicates if the intrinsic SMF
                                 // varies smoothly enough. See calc_smf_ssfr()
                                 // in calc_sfh.c.
  double sfr_sm[SM_BINS+SM_EXTRA]; // Average SFR for ***star-forming+quenched galaxies*** 
                                   // as a function of stellar mass.
  double sfr_sm_sf[SM_BINS+SM_EXTRA]; // Average SFR for ***star-forming galaxies*** 
                                   // as a function of stellar mass.
  double sfr_sm_q[SM_BINS+SM_EXTRA]; // Average SFR for ***quenched galaxies*** 
                                   // as a function of stellar mass.
  double cosmic_sfr, observed_cosmic_sfr;
  double cosmic_bhar, observed_cosmic_bhar;
  // double bh_growth[M_BINS]; // Deprecated.
  
  double k_uv[M_BINS]; // The slopes of the median UV mag--SFR relation for each halo mass bin.
  double b_uv[M_BINS]; // The intercepts of the median UV mag--SFR relation for each halo mass bin.
  double k_std_uv[M_BINS]; // The slopes of the UV mag scatter--halo mass relation for each halo mass bin.
  double b_std_uv[M_BINS]; // The intercepts of the UV mag scatter--halo mass relation for each halo mass bin.
  double obs_uv[M_BINS]; // The median observed UV magnitudes for each halo mass bin.
  double std_uv[M_BINS]; // The scatters of observed UV magnitudes for each halo mass bin.
  double std_uvlf[UV_BINS+UV_EXTRA]; // The scatters of observed UV magnitudes as a function of UV magnitudes.
  double uvlf[UV_BINS+UV_EXTRA]; // The observed UV luminosity function.
  char uvlf_ok[UV_BINS+UV_EXTRA]; // Flags indicating if the UV luminosity function varies smoothly.
                                  // See calc_uvlf() in calc_sfh.c.

  double log_bh_mass[M_BINS]; // The log10 of median BH mass for each halo mass bin.
  double bh_mass_avg[M_BINS]; // The average BH mass for each halo mass bin.
  double old_bh_mass[M_BINS]; // The average old BH mass for each halo mass bin.
  double old_bh_mass_self[M_BINS]; // The average old BH mass from the same halo mass bin int last snapshot.
  double old_bh_mass_next[M_BINS]; // The average old BH mass from the immediate smaller halo mass bin int last snapshot.
  double bh_eta_avg[M_BINS]; // The average total Eddington ratio for each halo mass bin.
  double bh_eta_rad_avg[MBH_BINS]; // The average radiative Eddington ratio for each halo mass bin.
  double bh_eta_kin_avg[M_BINS]; // The average kinetic Eddington ratio for each halo mass bin.
  double bh_acc_rate[M_BINS]; // Average BH accretion rate for each halo mass bin
  double bh_acc_rate_obs[M_BINS]; //BH accretion rates that induce a luminosity above the observable limit.
  double bh_merge_rate[M_BINS]; // Average BH merger rate for each halo mass bin
  double bh_eta[M_BINS]; // Typical BH Eddington ratio (the breaking point in the 
                         // double power-law) for each halo mass bin.
  double bh_f_occ[M_BINS]; //The SMBH occupation fraction for each halo mass bin.
  double bh_duty[M_BINS]; //The SMBH occupation fraction for each halo mass bin.
  double bh_f_occ_max[M_BINS];
  double dc_obs[M_BINS]; // The observed AGN duty cycle for each halo mass bin.
                          // See calc_observed_bhar() in calc_sfh.c.
  double bh_unmerged[M_BINS]; // The average unmerged (or wandering) BH mass for each halo mass bin.
  // double bh_merged[M_BINS];
  double new_bh_mass[M_BINS]; // The average newly grown BH mass for each halo mass bin.
  double bh_unmerged_avail[M_BINS]; //The total available unmerged BH mass (i.e. bh_unmerged after subtracting merged masses)
  double bh_unmerged_dist[M_BINS*MBH_BINS]; //The total available unmerged BH mass (i.e. bh_unmerged after subtracting merged masses)

  double n_merge10[M_BINS]; //The number of mergers with BH mass ratio > 1:10
  double n_merge100[M_BINS]; //The number of mergers with BH mass ratio > 1:100
  double f_CTK[M_BINS]; //The Compton-thick fractions.

  double f_active[MBH_BINS]; //Active BH fraction for each halo mass bin. 
                           // See calc_active_bh_fraction() in calc_sfh.c.
  double abhmf[MBH_BINS];
  // For the use of these gsl_spline objects, see calc_smf_and_ssfr() and calc_uvlf()
  // in calc_sfh.c.
  gsl_spline *spline; // Spline interpolation object for the stellar mass--halo mass 
                      // (SMHM) relation
  gsl_spline *spline2; //Since the SMHM relation may not be monotonic, we have to 
                      // split the SM-Mh relation into two pieces to do the interpolation.
  gsl_spline *spline_sfr; // Spline interpolation for the SFR--stellar mass relation.
  gsl_spline *spline_sfr2; // Similar to the SMHM relation, we might need a second segment.
  gsl_spline *spline_uv; // Spline interpolation for the UV--halo mass (UV--Mh) relation.
  gsl_spline *spline_uv2; //Since the relation is not monotonic, we have to split the UV--Mh relation into two pieces to do the interpolation.

  int alloc2smf; // The flags indicating if 2-segment interpolation have been allocated.
  int alloc2uvlf;

  int flag_alloc;
  double sm_from_icl[M_BINS];
  struct smf smhm;
};

void init_timesteps(void);
void calc_step_at_z(double z, int64_t *step, double *f);

#endif /* MLIST_H */
