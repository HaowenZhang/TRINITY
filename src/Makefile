CC=gcc
GSLFLAGS=-L$(HOME)/Software/gsl/2.6/lib -I$(HOME)/Software/gsl/2.6/include
CFLAGS=-lm -m64 -lgsl -lgslcblas -Wall -fopenmp -L$(GSL_LIB) -I$(GSL_INCLUDE) -L/opt/local/lib -I/opt/local/include -fno-math-errno
#CC=gcc
#CFLAGS=-lm -m64 -lgsl -lgslcblas -Wall -fopenmp
OPT_FLAGS=-O3 -msse3
DEBUG_FLAGS=-g -Wall -O0
PROF_FLAGS=-O3 -pg -ftest-coverage -fprofile-arcs
BASE_FILES=universe_time.c observations.c mt_rand.c jacobi.c integrate.c distance.c expcache2.c all_smf.c gen_mlists.c calc_sfh.c smoothing.c smloss.c check_syscalls.c mah.c merger_rates.c
SSFR_FILES = $(BASE_FILES) gen_ssfr.c
SFR_FILES = $(BASE_FILES) gen_sfr_check.c
SSFR2_FILES = $(BASE_FILES) gen_ssfr2.c
CPLOT_FILES = $(BASE_FILES) gen_colorplots.c
CSFR_FILES = $(BASE_FILES) gen_csfr.c
CATALOG_FILES= $(BASE_FILES) sm_catalog.c stringparse.c inthash.c 
AGN_UM_CATALOG_FILES = $(BASE_FILES) sm_agn_catalog_from_um.c stringparse.c 
ASMF_C_FILES = $(BASE_FILES) 
SMF_FILES = $(BASE_FILES) gen_smf.c 
UVLF_FILES = $(BASE_FILES) gen_uvlf.c
UVLF_TEST_FILES = $(BASE_FILES) gen_uvlf_test.c
QF_FILES = $(BASE_FILES) gen_qf.c
SMF_CHECK_FILES = $(BASE_FILES) gen_smf_check.c
SMASS_FILES = $(BASE_FILES) smass_density.c sm_limits.c
SMASS_COMP_FILES = $(BASE_FILES) smass_comparison_density.c sm_limits.c
SMASS_S_FILES = $(BASE_FILES) smass_hi_lo.c sm_limits.c
SMASS_P_FILES = $(BASE_FILES) smass_chi2.c sm_limits.c
SMASS_CV_FILES = $(BASE_FILES) smass_cv_cosmo.c 
SM2_FILES = $(BASE_FILES) smass_sm_avg.c sm_limits.c
FIT_FILES = $(BASE_FILES) fitter.c 
CHI2_FILES = $(BASE_FILES) calc_chi2.c
REMATCH_FILES = $(BASE_FILES)  rematch_smass.c sm_limits.c
COSMO_FILES = cosmo_volume.c distance.c

CATALOG_NAME=sm_catalog
AGN_UM_CATALOG_NAME=sm_agn_catalog_from_um
REMARG_NAME=marginalize_systematics
ASMF_OBJ_NAME=all_smf_mcmc
SMASS_NAME=smass
SSFR_NAME=gen_ssfr
SFR_NAME=gen_sfr_check
SSFR2_NAME=gen_ssfr2
CSFR_NAME=gen_csfr
SMASS_COMP_NAME=smass_comp
SMASS_S_NAME=smass_syst
SMASS_P_NAME=smass_chi2
SMASS_CV_NAME=smass_cv_cosmo
SM2_NAME=smass_sm_avg
SMF_NAME=gen_smf
UVLF_NAME=gen_uvlf
UVLF_TEST_NAME=gen_uvlf_test
SMF_CHECK_NAME=gen_smf_check
FIT_NAME=fitter
CHI2_NAME=calc_chi2
FIT2_NAME=fitter-mcmc
DEC_NAME=fiducial_deconvolute
COSMO_NAME = cosmo_volume
R_NAME = rematch_smass
QF_NAME = gen_qf

all:
	@make reg EXTRA_FLAGS="$(OPT_FLAGS)"

debug:
	@make reg EXTRA_FLAGS="$(DEBUG_FLAGS)"

prof:
	@make reg EXTRA_FLAGS="$(PROF_FLAGS)"

debug_csfr:
	@make csfr EXTRA_FLAGS="$(DEBUG_FLAGS)"

reg:
	$(CC) $(UVLF_TEST_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(UVLF_TEST_NAME)
	$(CC) $(BASE_FILES) gen_bh_f_occ.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_bh_f_occ
	#$(CC) $(BASE_FILES) gen_icl_map.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_icl_map
	#$(CC) mcmc_client.c $(BASE_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o mcmc_client
	#$(CC) $(CPLOT_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_colorplots
	$(CC) $(SSFR_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(SSFR_NAME)
	#$(CC) $(SFR_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(SFR_NAME)
	#$(CC) $(SSFR2_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(SSFR2_NAME)
	$(CC) $(CSFR_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(CSFR_NAME)
	$(CC) $(FIT_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(FIT_NAME)
	#$(CC) $(CHI2_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(CHI2_NAME)
	$(CC) $(ASMF_C_FILES) $(CFLAGS) $(EXTRA_FLAGS) -o $(ASMF_OBJ_NAME)
	$(CC) $(SMF_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(SMF_NAME)
	$(CC) $(UVLF_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(UVLF_NAME)
	$(CC) $(QF_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(QF_NAME)
	$(CC) $(SMF_CHECK_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(SMF_CHECK_NAME)
	#$(CC) $(BASE_FILES) gen_merger_tree.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_merger_tree
	$(CC) $(BASE_FILES) gen_galaxy_growth_rate.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_galaxy_growth_rate
	$(CC) $(BASE_FILES) gen_edd_r_color.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_edd_r_color
	$(CC) $(BASE_FILES) gen_binding_energy.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_binding_energy	
	$(CC) $(BASE_FILES) bh_host_halo_integral_JT.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bh_host_halo_integral_JT
	$(CC) $(BASE_FILES) bh_host_halo_JT.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bh_host_halo_JT
	$(CC) $(BASE_FILES) bh_host_galaxy_JT.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bh_host_galaxy_JT
	$(CC) $(BASE_FILES) bh_host_galaxy_halo_JT.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bh_host_galaxy_halo_JT
	$(CC) $(BASE_FILES) gen_max_mfrac.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_max_mfrac
	$(CC) $(BASE_FILES) gen_avg_eff_rad.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_avg_eff_rad
	$(CC) $(BASE_FILES) gen_merger_rate.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_merger_rate
	$(CC) $(BASE_FILES) bh_unmerged_dist.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bh_unmerged_dist
	$(CC) $(BASE_FILES) bher_dist.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bher_dist
	$(CC) $(BASE_FILES) test_bh_unmerged.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o test_bh_unmerged
	
	# $(CC) $(BASE_FILES) gen_sbhar_ssfr.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_sbhar_ssfr
	# $(CC) $(BASE_FILES) gen_edd_r_color_low_res.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_edd_r_color_low_res
	# $(CC) $(BASE_FILES) gen_sbhar_ssfr_low_res.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_sbhar_ssfr_low_res
	# $(CC) $(BASE_FILES) gen_qlf_split.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_qlf_split
	$(CC) $(BASE_FILES) gen_qlf.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_qlf
	$(CC) $(BASE_FILES) gen_qlf_mbh.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_qlf_mbh
	$(CC) $(BASE_FILES) gen_qlf_eta.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_qlf_eta
	$(CC) $(BASE_FILES) gen_qpdf.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_qpdf
	$(CC) $(BASE_FILES) bh_density.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bh_density
	$(CC) $(BASE_FILES) bhmf.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bhmf
	$(CC) $(BASE_FILES) gen_duty_cycle.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_duty_cycle
	$(CC) $(BASE_FILES) gen_qpdf_eta.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_qpdf_eta
	$(CC) $(BASE_FILES) gen_merger_check.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_merger_check
	$(CC) $(BASE_FILES) gen_bhar_sfr.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_bhar_sfr

	# $(CC) corr_hlist.c stringparse.c check_syscalls.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o corr_hlist
	# $(CC) corr_catalog.c check_syscalls.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o corr_catalog
	# $(CC) $(BASE_FILES) sm_limits.c smass_sm_med.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o smass_sm_med
	# $(CC) $(BASE_FILES) gen_smhm_smah.c sm_limits.c all_luminosities.c stringparse.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_smhm_smah
	# $(CC) sm_track.c sm_limits.c $(BASE_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o sm_track
	# $(CC) $(BASE_FILES) mf_integrate.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o mf_integrate
	# $(CC) $(BASE_FILES) gen_merging_sm.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_merging_sm
	# $(CC) $(BASE_FILES) gen_colorplots_planets.c  all_luminosities.c stringparse.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_colorplots_planets
	# $(CC) $(BASE_FILES) corr_testcam.c stringparse.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o corr_testcam
	# $(CC) $(BASE_FILES) gen_icl_fraction.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_icl_fraction
	# $(CC) $(BASE_FILES) deconv_fq.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o deconv_fq
	# $(CC) $(BASE_FILES) make_sf_catalog.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o make_sf_catalog
	# $(CC) $(BASE_FILES) make_sf_catalog_nc.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o make_sf_catalog_nc
	# $(CC) $(BASE_FILES) gen_colorplots_molly.c  all_luminosities.c stringparse.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_colorplots_molly
	# $(CC) print_doug_catalog.c check_syscalls.c -o print_doug_catalog
	# $(CC) convert_doug_catalog.c check_syscalls.c -o convert_doug_catalog
	# $(CC) print_sm_catalog.c check_syscalls.c -o print_sm_catalog
	# $(CC) shear_catalog.c check_syscalls.c stringparse.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o shear_catalog
	# $(CC) gen_sfr_nd.c $(BASE_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_sfr_nd
	# $(CC) shear.c check_syscalls.c stringparse.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o shear
	# $(CC) print_acc_list.c check_syscalls.c -o print_acc_list
	# $(CC) $(BASE_FILES) gen_lf.c all_luminosities.c stringparse.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_lf
	$(CC) $(BASE_FILES) gen_obs_ssfr_smah.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_obs_ssfr_smah
	# $(CC) $(BASE_FILES) gen_obs_sbhar_smah.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_obs_sbhar_smah
	# $(CC) $(BASE_FILES) gen_ssfr_hm.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_ssfr_hm
	# $(CC) $(BASE_FILES) cosmic_sfr_integral.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o cosmic_sfr_integral
	$(CC) $(BASE_FILES) gen_ssfr_smah_colorplot.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_ssfr_smah_colorplot
	# $(CC) corr_catalog_fast.c check_syscalls.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o corr_catalog_fast
	# $(CC) sm_track_surhud.c sm_limits.c $(BASE_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o sm_track_surhud
	# $(CC) itracks_trees_bryan.c $(BASE_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o itracks_trees_bryan
	# $(CC) individual_tracks_molly.c $(BASE_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o individual_tracks_molly
	# $(CC) make_graphs.c $(BASE_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -DCALC_ICL -o make_graphs
	# $(CC) $(BASE_FILES) z_lum_delta_dist.c  all_luminosities.c stringparse.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o z_lum_delta_dist
	# $(CC) $(BASE_FILES) gen_colorplots_grb.c all_luminosities.c stringparse.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_colorplots_grb
	# $(CC) $(BASE_FILES) andrey.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o andrey
	# $(CC) $(BASE_FILES) smhm_at_m_z.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o smhm_at_m_z
	# $(CC) $(BASE_FILES) make_fake_halo_catalog.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o make_fake_halo_catalog
	# $(CC) $(COSMO_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(COSMO_NAME)
	# $(CC) $(BASE_FILES) gen_lgsmhm.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_lgsmhm
	# $(CC) $(BASE_FILES) mar_threshold.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o mar_threshold
	# $(CC) $(BASE_FILES) gen_smar.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_smar
	$(CC) $(BASE_FILES) gen_ssfr_smah.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_ssfr_smah
	$(CC) $(BASE_FILES) gen_sbhar_smah.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_sbhar_smah
	# $(CC) $(BASE_FILES) gen_colorplots_grb_mb.c all_luminosities.c stringparse.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_colorplots_grb_mb
	# $(CC) $(BASE_FILES) gen_lum_plot.c stringparse.c all_luminosities.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_lum_plot
	# $(CC) $(BASE_FILES) stringparse.c inthash.c lightcone_catalog.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o lightcone_catalog
	# $(CC) individual_tracks_mark.c $(BASE_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o individual_tracks_mark
	# $(CC) $(SM2_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(SM2_NAME)
	# $(CC) $(BASE_FILES) number_density_conv.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o number_density_conv	
	# $(CC) $(BASE_FILES) icl_raw_conv_multi.c luminosities.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o icl_raw_conv_multi	
	# $(CC) $(CATALOG_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(CATALOG_NAME)
	# $(CC) $(BASE_FILES) icl_raw_conv.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o icl_raw_conv	
	# $(CC) $(BASE_FILES) gen_colorplots_nature.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_colorplots_nature
	# $(CC) gen_smf_intrinsic.c $(BASE_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_smf_intrinsic
	# $(CC) gen_smf_range.c $(BASE_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_smf_range
	# $(CC) smass_hi_lo.c sm_limits.c $(BASE_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o smass_syst
	# $(CC) gen_colorplots_yu.c $(BASE_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_colorplots_yu
	# $(CC) gen_colorplots_charlie.c $(BASE_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_colorplots_charlie
	# $(CC) sfr_frac.c $(BASE_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o sfr_frac
	# $(CC) calc_icl_fraction.c $(BASE_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o calc_icl_fraction
	# $(CC) calc_sm_density.c $(BASE_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o calc_sm_density
	# $(CC) gen_ma.c $(BASE_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_ma
	# $(CC) $(BASE_FILES) gen_colorplots_sandy.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_colorplots_sandy
	# $(CC) test_mah.c $(BASE_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o test_mah
	# $(CC) individual_tracks.c $(BASE_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o individual_tracks
	# $(CC) fixup_sfr_errors.c $(BASE_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o fixup_sfr_errors
	# $(CC) fitter_mcmc.c $(BASE_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o fitter_mcmc
	# $(CC) range_check.c $(BASE_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o range_check
	# $(CC) range_check_alpha.c $(BASE_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o range_check_alpha
	# $(CC) range_check_delta.c $(BASE_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o range_check_delta
	# $(CC) quickie_mcmc.c $(BASE_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o quickie_mcmc
	# $(CC) calc_smage.c $(BASE_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o calc_smage
	# $(CC) gen_ma_mass.c $(BASE_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_ma_mass
	# $(CC) smass_null_evolution.c sm_limits.c $(BASE_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o smass_null_evolution
	# $(CC) parameter_range.c $(BASE_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o parameter_range
	# $(CC) smass_hi_lo_then.c sm_limits.c $(BASE_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o smass_syst_then
	# $(CC) $(SMASS_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(SMASS_NAME)
	# $(CC) gen_sfr_z.c $(BASE_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_sfr_z

#	$(CC) $(REMARG_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(REMARG_NAME)
#	$(CC) $(CACHE_FILES) $(CFLAGS) $(EXTRA_FLAGS) -o $(CACHE_OBJ_NAME)
#	$(CC) $(REMATCH_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(R_NAME)
#	$(CC) $(SMASS_COMP_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(SMASS_COMP_NAME)
#	$(CC) $(SMASS_P_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(SMASS_P_NAME)
#	$(CC) $(SMASS_CV_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(SMASS_CV_NAME)
#	$(CC) $(FIT2_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(FIT2_NAME)
gen_sm_agn_catalog_from_um:
	$(CC) $(AGN_UM_CATALOG_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(AGN_UM_CATALOG_NAME)
bher_dist_check:
	$(CC) $(BASE_FILES) bher_dist_check.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_bher_dist_check
lum_dist_check:
	$(CC) $(BASE_FILES) lum_dist_check.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_lum_dist_check
baryon_eff:
	$(CC) $(BASE_FILES) gen_baryon_eff.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_baryon_eff
bhar_mbh:
	$(CC) $(BASE_FILES) bhar_mbh.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bhar_mbh
bhar_sfr_mbh:
	$(CC) $(BASE_FILES) bhar_sfr_mbh.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bhar_sfr_mbh
sbhar_ssfr_mbh:
	$(CC) $(BASE_FILES) sbhar_ssfr_mbh.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o sbhar_ssfr_mbh
bher_mbh:
	$(CC) $(BASE_FILES) bher_mbh.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bher_mbh
bhmr_mbh:
	$(CC) $(BASE_FILES) bhmr_mbh.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bhmr_mbh
bhar_mstar:
	$(CC) $(BASE_FILES) bhar_mstar.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bhar_mstar
bhar_sfr_mstar:
	$(CC) $(BASE_FILES) bhar_sfr_mstar.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bhar_sfr_mstar
smar_ssfr_sbhar:
	$(CC) $(BASE_FILES) smar_ssfr_sbhar.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o smar_ssfr_sbhar
sbhar_mstar:
	$(CC) $(BASE_FILES) sbhar_mstar.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o sbhar_mstar
sbhar_ssfr_mstar:
	$(CC) $(BASE_FILES) sbhar_ssfr_mstar.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o sbhar_ssfr_mstar
bh_f_occ:
	$(CC) $(BASE_FILES) gen_bh_f_occ.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_bh_f_occ
bher_mstar:
	$(CC) $(BASE_FILES) bher_mstar.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bher_mstar
bhmr_mstar:
	$(CC) $(BASE_FILES) bhmr_mstar.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bhmr_mstar
binding_energy:
	$(CC) $(BASE_FILES) gen_binding_energy.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_binding_energy	
hist2d:
	$(CC) $(BASE_FILES) gen_edd_r_color.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_edd_r_color
gen_sm_sfr_check:
	$(CC) $(BASE_FILES) sm_sfr_check.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o sm_sfr_check
chi2_type:
	$(CC) $(BASE_FILES) gen_chi2_type.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_chi2_type
asymmetry_check:
	$(CC) $(BASE_FILES) gen_asymmetry_check.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_asymmetry_check
gen_energy_density:
	$(CC) $(BASE_FILES) energy_density.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o energy_density
qf:
	$(CC) $(QF_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(QF_NAME)

sfr_map:
	$(CC) $(BASE_FILES) gen_sfr_map.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_sfr_map
uvlf:
	$(CC) $(UVLF_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(UVLF_NAME)
uvlf_test:
	$(CC) $(UVLF_TEST_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(UVLF_TEST_NAME)
uv_map:
	$(CC) $(BASE_FILES) gen_uv_map.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_uv_map

fitter:
	$(CC) $(FIT_FILES) $(CFLAGS) $(OPT_FLAGS) -DGEN_SMF -o $(FIT_NAME)
asmf:
	$(CC) $(ASMF_C_FILES) $(CFLAGS) $(OPT_FLAGS) -o $(ASMF_OBJ_NAME)
csfr:
	$(CC) $(CSFR_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(CSFR_NAME)
ssfr:
	$(CC) $(SSFR_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(SSFR_NAME)
duty_cycle:
	$(CC) $(BASE_FILES) gen_duty_cycle.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_duty_cycle
frac_supEdd_l:
	$(CC) $(BASE_FILES) gen_frac_supEdd_l.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_frac_supEdd_l
f_active:
	$(CC) $(BASE_FILES) gen_active_fraction.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_active_fraction
sbhar_smah:
	$(CC) $(BASE_FILES) gen_sbhar_smah.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_sbhar_smah
ssfr_smah:
	$(CC) $(BASE_FILES) gen_ssfr_smah.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_ssfr_smah
gen_bhmf_all_z:
	$(CC) $(BASE_FILES) bhmf_all_z.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bhmf_all_z
gen_bhmf_mh:
	$(CC) $(BASE_FILES) bhmf_mh.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bhmf_mh
asymmetry_check:
	$(CC) $(BASE_FILES) gen_asymmetry_check.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_asymmetry_check
erdf:
	$(CC) $(BASE_FILES) gen_erdf.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_erdf
gen_bh_density_split:
	$(CC) $(BASE_FILES) bh_density_split.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bh_density_split
gen_bh_density_split_check:
	$(CC) $(BASE_FILES) bh_density_split_check.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bh_density_split_check
gen_quasar_host_halo_Feige:
	$(CC) $(BASE_FILES) quasar_host_halo_Feige.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o quasar_host_halo_Feige
gen_bhmf_quasar_lum:
	$(CC) $(BASE_FILES) bhmf_quasar_lum.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bhmf_quasar_lum
gen_host_smf_quasar_lum:
	$(CC) $(BASE_FILES) host_smf_quasar_lum.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o host_smf_quasar_lum
gen_host_hmf_quasar_lum:
	$(CC) $(BASE_FILES) host_hmf_quasar_lum.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o host_hmf_quasar_lum
gen_lx_sfr_mstar_james:
	$(CC) $(BASE_FILES) lx_sfr_mstar_james.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o lx_sfr_mstar_james
gen_lx_mstar_james:
	$(CC) $(BASE_FILES) lx_mstar_james.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o lx_mstar_james
gen_bhmf_typeI:
	$(CC) $(BASE_FILES) bhmf_typeI.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bhmf_typeI
qpdf_eta:
	$(CC) $(BASE_FILES) gen_qpdf_eta.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_qpdf_eta
qlf_mbh:
	$(CC) $(BASE_FILES) gen_qlf_mbh.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_qlf_mbh
qlf_eta:
	$(CC) $(BASE_FILES) gen_qlf_eta.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_qlf_eta
qlf:
	$(CC) $(BASE_FILES) gen_qlf.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_qlf
lum_dist_mbh_obs:
	$(CC) $(BASE_FILES) gen_lum_dist_mbh_obs.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_lum_dist_mbh_obs
lum_dist_mbh_int:
	$(CC) $(BASE_FILES) gen_lum_dist_mbh_int.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_lum_dist_mbh_int
gen_bh_density:
	$(CC) $(BASE_FILES) bh_density.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bh_density
bhmf_lum_dist:
	$(CC) $(BASE_FILES) output_bhmf_lum_dist.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o output_bhmf_lum_dist
clean:
	rm -f $(CACHE_OBJ_NAME) $(ASMF_OBJ_NAME) $(SMF_NAME) $(R_NAME) $(SMASS_NAME) $(SM2_NAME) $(FIT_NAME) $(DEC_NAME) $(COSMO_NAME) $(SMASS_CV_NAME) $(FIT2_NAME) $(REMARG_NAME)
	rm -f *~
	rm -f gmon.out *.gcda *.gcno 
	rm:wq
	-rf *.dSYM
