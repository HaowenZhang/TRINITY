CC=gcc
GSLFLAGS=-L$(HOME)/Software/gsl/2.6/lib -I$(HOME)/Software/gsl/2.6/include
CFLAGS=-lm -m64 -lgsl -lgslcblas -Wall -fopenmp -L$(GSL_LIB) -I$(GSL_INCLUDE) -L/opt/local/lib -I/opt/local/include -fno-math-errno
#CC=gcc
#CFLAGS=-lm -m64 -lgsl -lgslcblas -Wall -fopenmp
OPT_FLAGS=-O3 -msse3
DEBUG_FLAGS=-g -Wall -O0
PROF_FLAGS=-O3 -pg -ftest-coverage -fprofile-arcs
DIR_BIN = ../bin
DIR_BASE=./base
DIR_MCMC=./mcmc
DIR_OTHER=./others
DIR_PAPER1=./paper1
DIR_PAPER2=./paper2
DIR_PAPER3=./paper3
DIR_PAPER4=./paper4
BASE_FILES=$(DIR_BASE)/universe_time.c $(DIR_BASE)/observations.c $(DIR_BASE)/mt_rand.c $(DIR_BASE)/jacobi.c $(DIR_BASE)/integrate.c $(DIR_BASE)/distance.c $(DIR_BASE)/expcache2.c $(DIR_BASE)/all_smf.c $(DIR_BASE)/gen_mlists.c $(DIR_BASE)/calc_sfh.c $(DIR_BASE)/smoothing.c $(DIR_BASE)/smloss.c $(DIR_BASE)/check_syscalls.c $(DIR_BASE)/mah.c $(DIR_BASE)/merger_rates.c $(DIR_BASE)/sm_limits.c
# BASE_FILES=universe_time.c observations.c mt_rand.c jacobi.c integrate.c distance.c expcache2.c all_smf.c gen_mlists.c calc_sfh.c smoothing.c smloss.c check_syscalls.c mah.c merger_rates.c
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
SMASS_FILES = $(BASE_FILES) smass_density.c
SMASS_COMP_FILES = $(BASE_FILES) smass_comparison_density.c
SMASS_S_FILES = $(BASE_FILES) smass_hi_lo.c
SMASS_P_FILES = $(BASE_FILES) smass_chi2.c
SMASS_CV_FILES = $(BASE_FILES) smass_cv_cosmo.c 
SM2_FILES = $(BASE_FILES) smass_sm_avg.c
FIT_FILES = $(BASE_FILES) fitter.c 
CHI2_FILES = $(BASE_FILES) calc_chi2.c
REMATCH_FILES = $(BASE_FILES)  rematch_smass.c
COSMO_FILES = cosmo_volume.c distance.c

CATALOG_NAME=sm_catalog
AGN_UM_CATALOG_NAME=sm_agn_catalog_from_um
# REMARG_NAME=marginalize_systematics
ASMF_OBJ_NAME=$(DIR_BIN)/all_smf_mcmc
# SMASS_NAME=smass
SSFR_NAME=gen_ssfr
SFR_NAME=gen_sfr_check
SSFR2_NAME=gen_ssfr2
CSFR_NAME=gen_csfr
SMASS_COMP_NAME=smass_comp
SMASS_S_NAME=smass_syst
SMASS_P_NAME=smass_chi2
# SMASS_CV_NAME=smass_cv_cosmo
# SM2_NAME=smass_sm_avg
SMF_NAME=gen_smf
UVLF_NAME=gen_uvlf
UVLF_TEST_NAME=gen_uvlf_test
SMF_CHECK_NAME=gen_smf_check
FIT_NAME=$(DIR_BIN)/fitter
CHI2_NAME=$(DIR_BIN)calc_chi2
# FIT2_NAME=fitter-mcmc
# DEC_NAME=fiducial_deconvolute
# COSMO_NAME = cosmo_volume
# R_NAME = rematch_smass
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
	$(CC) $(BASE_FILES) gen_bh_f_occ.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_bh_f_occ
	$(CC) $(FIT_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(FIT_NAME)
	$(CC) $(ASMF_C_FILES) $(CFLAGS) $(EXTRA_FLAGS) -o $(ASMF_OBJ_NAME)
	$(CC) $(SMF_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(SMF_NAME)
	$(CC) $(UVLF_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(UVLF_NAME)
	$(CC) $(QF_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(QF_NAME)
	$(CC) $(SMF_CHECK_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(SMF_CHECK_NAME)
	$(CC) $(BASE_FILES) gen_galaxy_growth_rate.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_galaxy_growth_rate
	$(CC) $(BASE_FILES) gen_edd_r_color.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_edd_r_color
	$(CC) $(BASE_FILES) gen_binding_energy.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_binding_energy	
	$(CC) $(BASE_FILES) bh_host_halo_integral_JT.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bh_host_halo_integral_JT
	$(CC) $(BASE_FILES) bh_host_halo_JT.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bh_host_halo_JT
	$(CC) $(BASE_FILES) bh_host_galaxy_JT.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bh_host_galaxy_JT
	$(CC) $(BASE_FILES) bh_host_galaxy_halo_JT.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bh_host_galaxy_halo_JT
	$(CC) $(BASE_FILES) gen_max_mfrac.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_max_mfrac
	$(CC) $(BASE_FILES) gen_avg_eff_rad.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_avg_eff_rad
	$(CC) $(BASE_FILES) gen_merger_rate.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_merger_rate
	$(CC) $(BASE_FILES) bh_unmerged_dist.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bh_unmerged_dist
	$(CC) $(BASE_FILES) bher_dist.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bher_dist
	$(CC) $(BASE_FILES) test_bh_unmerged.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o test_bh_unmerged

	$(CC) $(BASE_FILES) $(DIR_PAPER1)/gen_qlf.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(DIR_BIN)/gen_qlf
	$(CC) $(BASE_FILES) gen_qlf_mbh.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_qlf_mbh
	$(CC) $(BASE_FILES) gen_qlf_eta.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_qlf_eta
	$(CC) $(BASE_FILES) gen_qpdf.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_qpdf
	$(CC) $(BASE_FILES) bh_density.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bh_density
	$(CC) $(BASE_FILES) bhmf.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bhmf
	$(CC) $(BASE_FILES) gen_duty_cycle.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_duty_cycle
	$(CC) $(BASE_FILES) gen_qpdf_eta.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_qpdf_eta
	$(CC) $(BASE_FILES) gen_merger_check.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_merger_check
	$(CC) $(BASE_FILES) gen_bhar_sfr.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_bhar_sfr

	$(CC) $(BASE_FILES) gen_obs_ssfr_smah.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_obs_ssfr_smah
	$(CC) $(BASE_FILES) gen_ssfr_smah_colorplot.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_ssfr_smah_colorplot

	$(CC) $(BASE_FILES) gen_ssfr_smah.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_ssfr_smah
	$(CC) $(BASE_FILES) gen_sbhar_smah.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_sbhar_smah



paper2_plots:
	$(CC) $(BASE_FILES) $(DIR_PAPER2)/mbh_perc_mstar_lbol_individual.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(DIR_BIN)/mbh_perc_mstar_lbol_individual
	$(CC) $(BASE_FILES) $(DIR_PAPER2)/mbh_perc_mstar_lbol.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(DIR_BIN)/mbh_perc_mstar_lbol

gen_bhar_kin_rad_ratio:
	$(CC) $(BASE_FILES) bhar_kin_rad_ratio.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bhar_kin_rad_ratio
gen_bhar_hist_eta:
	$(CC) $(BASE_FILES) bhar_hist_eta.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bhar_hist_eta
gen_kin_power_mbh:
	$(CC) $(BASE_FILES) kin_power_mbh.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o kin_power_mbh
gen_mbh_perc_mstar_lbol_individual:
	$(CC) $(BASE_FILES) $(DIR_PAPER2)/mbh_perc_mstar_lbol_individual.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(DIR_BIN)/mbh_perc_mstar_lbol_individual
gen_mbh_perc_mstar_lbol:
	$(CC) $(BASE_FILES) $(DIR_PAPER2)/mbh_perc_mstar_lbol.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(DIR_BIN)/mbh_perc_mstar_lbol
gen_bhmf:
	$(CC) $(BASE_FILES) bhmf.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bhmf
gen_sm_agn_catalog_from_um:
	$(CC) $(AGN_UM_CATALOG_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(AGN_UM_CATALOG_NAME)
bher_dist_check:
	$(CC) $(BASE_FILES) bher_dist_check.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_bher_dist_check
lum_dist_check:
	$(CC) $(BASE_FILES) lum_dist_check.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_lum_dist_check
baryon_eff:
	$(CC) $(BASE_FILES) gen_baryon_eff.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_baryon_eff
qlf_mstar_bias:
	$(CC) $(BASE_FILES) gen_qlf_mstar_bias.c $(CFLAGS) $(OPT_FLAGS) -DGEN_SMF -o gen_qlf_mstar_bias
bhar_mbh:
	$(CC) $(BASE_FILES) bhar_mbh.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bhar_mbh
bhar_sfr_mbh:
	$(CC) $(BASE_FILES) bhar_sfr_mbh.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bhar_sfr_mbh
sbhar_ssfr_mbh:
	$(CC) $(BASE_FILES) sbhar_ssfr_mbh.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o sbhar_ssfr_mbh
bher_mbh:
	$(CC) $(BASE_FILES) bher_mbh.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bher_mbh
gen_bher_mbh_step:
	$(CC) $(BASE_FILES) bher_mbh_step.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bher_mbh_step
gen_bher_med_active_mh:
	$(CC) $(BASE_FILES) bher_med_active_mh.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bher_med_active_mh
gen_bher_active_mstar:
	$(CC) $(BASE_FILES) bher_active_mstar.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bher_active_mstar
gen_average_rad_bher_active_mstar:
	$(CC) $(BASE_FILES) average_rad_bher_active_mstar.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o average_rad_bher_active_mstar
gen_average_rad_bher_active_mh:
	$(CC) $(BASE_FILES) average_rad_bher_active_mh.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o average_rad_bher_active_mh
gen_bher_med_active_mstar:
	$(CC) $(BASE_FILES) bher_med_active_mstar.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bher_med_active_mstar
bhmr_mbh:
	$(CC) $(BASE_FILES) bhmr_mbh.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bhmr_mbh
bhar_mstar:
	$(CC) $(BASE_FILES) bhar_mstar.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bhar_mstar
bhar_sfr_mstar:
	$(CC) $(BASE_FILES) bhar_sfr_mstar.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bhar_sfr_mstar
smar_ssfr_sbhar:
	$(CC) $(BASE_FILES) smar_ssfr_sbhar.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o smar_ssfr_sbhar
sbhar_mstar:
	$(CC) $(BASE_FILES) sbhar_mstar.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o sbhar_mstar
sbhar_ssfr_mstar:
	$(CC) $(BASE_FILES) sbhar_ssfr_mstar.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o sbhar_ssfr_mstar
bh_f_occ:
	$(CC) $(BASE_FILES) gen_bh_f_occ.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_bh_f_occ
bher_mstar:
	$(CC) $(BASE_FILES) bher_mstar.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bher_mstar
bher_mh:
	$(CC) $(BASE_FILES) bher_mh.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bher_mh
bhmr_mstar:
	$(CC) $(BASE_FILES) bhmr_mstar.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bhmr_mstar
binding_energy:
	$(CC) $(BASE_FILES) gen_binding_energy.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_binding_energy	
hist2d:
	$(CC) $(BASE_FILES) gen_edd_r_color.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_edd_r_color
gen_sm_sfr_check:
	$(CC) $(BASE_FILES) sm_sfr_check.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o sm_sfr_check
chi2_type:
	$(CC) $(BASE_FILES) gen_chi2_type.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_chi2_type
gen_hmf_mbh_limit:
	$(CC) $(BASE_FILES) hmf_mbh_limit.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o hmf_mbh_limit
asymmetry_check:
	$(CC) $(BASE_FILES) gen_asymmetry_check.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_asymmetry_check
asymmetry_check_mag:
	$(CC) $(BASE_FILES) gen_asymmetry_check_mag.c $(CFLAGS) $(DEBUG_FLAGS) -DGEN_SMF -o gen_asymmetry_check_mag
asymmetry_check_lum:
	$(CC) $(BASE_FILES) gen_asymmetry_check_lum.c $(CFLAGS) $(DEBUG_FLAGS) -DGEN_SMF -o gen_asymmetry_check_lum
gen_energy_density_mstar:
	$(CC) $(BASE_FILES) energy_density_mstar.c $(CFLAGS) $(OPT_FLAGS) -DGEN_SMF -o energy_density_mstar
gen_energy_density:
	$(CC) $(BASE_FILES) energy_density.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o energy_density
qf:
	$(CC) $(QF_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(QF_NAME)

sfr_map:
	$(CC) $(BASE_FILES) gen_sfr_map.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_sfr_map
uvlf:
	$(CC) $(UVLF_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(UVLF_NAME)
uvlf_test:
	$(CC) $(UVLF_TEST_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(UVLF_TEST_NAME)
uv_map:
	$(CC) $(BASE_FILES) gen_uv_map.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_uv_map

fitter:
	$(CC) $(FIT_FILES) $(CFLAGS) $(OPT_FLAGS) -DGEN_SMF -o $(FIT_NAME)
asmf:
	$(CC) $(ASMF_C_FILES) $(CFLAGS) $(OPT_FLAGS) -o $(ASMF_OBJ_NAME)
csfr:
	$(CC) $(CSFR_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(CSFR_NAME)
ssfr:
	$(CC) $(SSFR_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(SSFR_NAME)

lkin_avg_at_lx:
	$(CC) $(BASE_FILES) gen_lkin_avg_at_lx.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_lkin_avg_at_lx
frac_supEdd_l:
	$(CC) $(BASE_FILES) gen_frac_supEdd_l.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_frac_supEdd_l
f_active:
	$(CC) $(BASE_FILES) gen_active_fraction.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_active_fraction
sbhar_smah:
	$(CC) $(BASE_FILES) gen_sbhar_smah.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_sbhar_smah
ssfr_smah:
	$(CC) $(BASE_FILES) gen_ssfr_smah.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_ssfr_smah
gen_bhmf_all_z:
	$(CC) $(BASE_FILES) bhmf_all_z.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bhmf_all_z
gen_bhmf_mh:
	$(CC) $(BASE_FILES) bhmf_mh.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bhmf_mh
erdf:
	$(CC) $(BASE_FILES) gen_erdf.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_erdf
gen_bh_density_split:
	$(CC) $(BASE_FILES) bh_density_split.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bh_density_split
gen_bh_density_split_check:
	$(CC) $(BASE_FILES) bh_density_split_check.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bh_density_split_check
gen_quasar_host_halo_Feige:
	$(CC) $(BASE_FILES) quasar_host_halo_Feige.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o quasar_host_halo_Feige
gen_bhmf_quasar_lum:
	$(CC) $(BASE_FILES) bhmf_quasar_lum.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bhmf_quasar_lum
gen_host_smf_quasar_lum:
	$(CC) $(BASE_FILES) host_smf_quasar_lum.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o host_smf_quasar_lum
gen_host_hmf_quasar_lum:
	$(CC) $(BASE_FILES) host_hmf_quasar_lum.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o host_hmf_quasar_lum
gen_lx_sfr_mstar_james:
	$(CC) $(BASE_FILES) lx_sfr_mstar_james.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o lx_sfr_mstar_james
gen_lx_mstar_james:
	$(CC) $(BASE_FILES) lx_mstar_james.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o lx_mstar_james
gen_bhmf_typeI:
	$(CC) $(BASE_FILES) bhmf_typeI.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bhmf_typeI
qpdf_eta:
	$(CC) $(BASE_FILES) gen_qpdf_eta.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_qpdf_eta
qlf_mstar:
	$(CC) $(BASE_FILES) gen_qlf_mstar.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_qlf_mstar
frac_agn_lum_mstar:
	$(CC) $(BASE_FILES) gen_frac_agn_lum_mstar.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_frac_agn_lum_mstar
mstar_perc_quasars:
	$(CC) $(BASE_FILES) gen_mstar_perc_quasars.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_mstar_perc_quasars
mstar_avg_quasars:
	$(CC) $(BASE_FILES) gen_mstar_avg_quasars.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_mstar_avg_quasars
qlf_mh:
	$(CC) $(BASE_FILES) gen_qlf_mh.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_qlf_mh
qlf_mbh:
	$(CC) $(BASE_FILES) gen_qlf_mbh.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_qlf_mbh
qlf_eta:
	$(CC) $(BASE_FILES) gen_qlf_eta.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_qlf_eta
qlf:
	$(CC) $(BASE_FILES) $(DIR_PAPER1)/gen_qlf.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(DIR_BIN)/gen_qlf
lum_dist_mbh_obs:
	$(CC) $(BASE_FILES) gen_lum_dist_mbh_obs.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_lum_dist_mbh_obs
lum_dist_mbh_int:
	$(CC) $(BASE_FILES) gen_lum_dist_mbh_int.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_lum_dist_mbh_int
gen_bh_density:
	$(CC) $(BASE_FILES) bh_density.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bh_density
bhmf_lum_dist:
	$(CC) $(BASE_FILES) output_bhmf_lum_dist.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o output_bhmf_lum_dist
clean:
# 	rm -f $(CACHE_OBJ_NAME) $(ASMF_OBJ_NAME) $(SMF_NAME) $(R_NAME) $(SMASS_NAME) $(SM2_NAME) $(FIT_NAME) $(DEC_NAME) $(COSMO_NAME) $(SMASS_CV_NAME) $(FIT2_NAME) $(REMARG_NAME)
	rm -f $(CACHE_OBJ_NAME) $(ASMF_OBJ_NAME) $(SMF_NAME) $(FIT_NAME)

	rm -f *~
	rm -f gmon.out *.gcda *.gcno 
	rm:wq
	-rf *.dSYM
