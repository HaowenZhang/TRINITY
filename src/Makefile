CC=gcc
GSL_LIB=$(HOME)/Software/gsl/2.6/lib 
GSL_INCLUDE=$(HOME)/Software/gsl/2.6/include
CFLAGS=-lm -m64 -lgsl -lgslcblas -Wall -fopenmp -L$(GSL_LIB) -I$(GSL_INCLUDE) -L/opt/local/lib -I/opt/local/include -fno-math-errno
OPT_FLAGS=-O3 -msse3
DEBUG_FLAGS=-g -Wall -O0
PROF_FLAGS=-O3 -pg -ftest-coverage -fprofile-arcs
BASE_FILES=universe_time.c observations.c mt_rand.c jacobi.c integrate.c distance.c expcache2.c all_smf.c gen_mlists.c calc_sfh.c smoothing.c smloss.c check_syscalls.c mah.c merger_rates.c
SSFR_FILES = $(BASE_FILES) gen_ssfr.c
SFR_FILES = $(BASE_FILES) gen_sfr_check.c
SSFR2_FILES = $(BASE_FILES) gen_ssfr2.c
CSFR_FILES = $(BASE_FILES) gen_csfr.c
CATALOG_FILES= $(BASE_FILES) sm_catalog.c stringparse.c inthash.c 
AGN_UM_CATALOG_FILES = $(BASE_FILES) sm_agn_catalog_from_um.c stringparse.c 
ASMF_C_FILES = $(BASE_FILES) 
SMF_FILES = $(BASE_FILES) gen_smf.c 
UVLF_FILES = $(BASE_FILES) gen_uvlf.c
QF_FILES = $(BASE_FILES) gen_qf.c
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
SSFR2_NAME=gen_ssfr2
CSFR_NAME=gen_csfr
SMASS_COMP_NAME=smass_comp
SMASS_S_NAME=smass_syst
SMASS_P_NAME=smass_chi2
SMASS_CV_NAME=smass_cv_cosmo
SM2_NAME=smass_sm_avg
SMF_NAME=gen_smf
UVLF_NAME=gen_uvlf
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
	$(CC) $(SSFR_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(SSFR_NAME)
	$(CC) $(CSFR_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(CSFR_NAME)
	$(CC) $(FIT_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(FIT_NAME)
	$(CC) $(CHI2_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(CHI2_NAME)
	$(CC) $(ASMF_C_FILES) $(CFLAGS) $(EXTRA_FLAGS) -o $(ASMF_OBJ_NAME)
	$(CC) $(SMF_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(SMF_NAME)
	$(CC) $(UVLF_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(UVLF_NAME)
	$(CC) $(QF_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(QF_NAME)
	$(CC) $(BASE_FILES) gen_merger_tree.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_merger_tree
	$(CC) $(BASE_FILES) gen_galaxy_growth_rate.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_galaxy_growth_rate
	$(CC) $(BASE_FILES) gen_edd_r_color.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_edd_r_color
	$(CC) $(BASE_FILES) bh_host_halo_integral_JT.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bh_host_halo_integral_JT
	$(CC) $(BASE_FILES) bh_host_halo_JT.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bh_host_halo_JT
	$(CC) $(BASE_FILES) bh_host_galaxy_JT.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bh_host_galaxy_JT
	$(CC) $(BASE_FILES) gen_max_mfrac.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_max_mfrac
	$(CC) $(BASE_FILES) gen_avg_eff_rad.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_avg_eff_rad
	$(CC) $(BASE_FILES) gen_merger_rate.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_merger_rate
	$(CC) $(BASE_FILES) bh_unmerged_dist.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bh_unmerged_dist
	
	$(CC) $(BASE_FILES) gen_qlf_mstar.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_qlf_mstar
	$(CC) $(BASE_FILES) gen_qlf.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_qlf
	$(CC) $(BASE_FILES) gen_qlf_mbh.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_qlf_mbh
	$(CC) $(BASE_FILES) gen_qlf_eta.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_qlf_eta
	$(CC) $(BASE_FILES) gen_qpdf.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_qpdf
	$(CC) $(BASE_FILES) bh_density.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bh_density
	$(CC) $(BASE_FILES) bhmf.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bhmf
	$(CC) $(BASE_FILES) gen_duty_cycle.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_duty_cycle
	$(CC) $(BASE_FILES) gen_qpdf_eta.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_qpdf_eta
	$(CC) $(BASE_FILES) gen_merger_check.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_merger_check

	$(CC) $(BASE_FILES) gen_ssfr_hm.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_ssfr_hm
	$(CC) $(BASE_FILES) cosmic_sfr_integral.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o cosmic_sfr_integral
	$(CC) $(COSMO_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(COSMO_NAME)
	$(CC) $(BASE_FILES) gen_ssfr_smah.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_ssfr_smah
	$(CC) $(BASE_FILES) gen_sbhar_smah.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_sbhar_smah
	$(CC) gen_smf_intrinsic.c $(BASE_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_smf_intrinsic
	$(CC) gen_smf_range.c $(BASE_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_smf_range
	$(CC) calc_sm_density.c $(BASE_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o calc_sm_density


merger_tree:
	$(CC) $(BASE_FILES) gen_merger_tree.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_merger_tree
max_mfrac:
	$(CC) $(BASE_FILES) gen_max_mfrac.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_max_mfrac
merger_rate:
	$(CC) $(BASE_FILES) gen_merger_rate.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_merger_rate
gen_sm_agn_catalog_from_um:
	$(CC) $(AGN_UM_CATALOG_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(AGN_UM_CATALOG_NAME)
bher_dist_check:
	$(CC) $(BASE_FILES) bher_dist_check.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_bher_dist_check
lum_dist_check:
	$(CC) $(BASE_FILES) lum_dist_check.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_lum_dist_check
baryon_eff:
	$(CC) $(BASE_FILES) gen_baryon_eff.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_baryon_eff
gen_bhmf_mh:
	$(CC) $(BASE_FILES) bhmf_mh.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bhmf_mh
gen_bhar_mbh:
	$(CC) $(BASE_FILES) bhar_mbh.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bhar_mbh
gen_bhar_sfr_mbh:
	$(CC) $(BASE_FILES) bhar_sfr_mbh.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bhar_sfr_mbh
gen_sbhar_ssfr_mbh:
	$(CC) $(BASE_FILES) sbhar_ssfr_mbh.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o sbhar_ssfr_mbh
gen_bher_mbh:
	$(CC) $(BASE_FILES) bher_mbh.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bher_mbh
gen_bhmr_mbh:
	$(CC) $(BASE_FILES) bhmr_mbh.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bhmr_mbh
gen_bhar_mstar:
	$(CC) $(BASE_FILES) bhar_mstar.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bhar_mstar
gen_bhar_sfr_mstar:
	$(CC) $(BASE_FILES) bhar_sfr_mstar.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bhar_sfr_mstar
gen_smar_ssfr_sbhar:
	$(CC) $(BASE_FILES) smar_ssfr_sbhar.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o smar_ssfr_sbhar
gen_sbhar_mstar:
	$(CC) $(BASE_FILES) sbhar_mstar.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o sbhar_mstar
gen_sbhar_ssfr_mstar:
	$(CC) $(BASE_FILES) sbhar_ssfr_mstar.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o sbhar_ssfr_mstar
gen_bher_mstar:
	$(CC) $(BASE_FILES) bher_mstar.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bher_mstar
gen_bhmr_mstar:
	$(CC) $(BASE_FILES) bhmr_mstar.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bhmr_mstar
galaxy_growth_rate:
	$(CC) $(BASE_FILES) gen_galaxy_growth_rate.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_galaxy_growth_rate
hist2d:
	$(CC) $(BASE_FILES) gen_edd_r_color.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_edd_r_color
gen_sm_sfr_check:
	$(CC) $(BASE_FILES) sm_sfr_check.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o sm_sfr_check
chi2_type:
	$(CC) $(BASE_FILES) gen_chi2_type.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_chi2_type
mbh_info:
	$(CC) $(BASE_FILES) gen_mbh_info.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_mbh_info
gen_energy_density:
	$(CC) $(BASE_FILES) energy_density.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o energy_density
avg_eff_rad:
	$(CC) $(BASE_FILES) gen_avg_eff_rad.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_avg_eff_rad
qf:
	$(CC) $(QF_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(QF_NAME)
gen_bh_host_galaxy_JT:
	$(CC) $(BASE_FILES) bh_host_galaxy_JT.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bh_host_galaxy_JT
gen_bh_host_halo_JT:
	$(CC) $(BASE_FILES) bh_host_halo_JT.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bh_host_halo_JT
gen_bh_host_halo_integral_JT:
	$(CC) $(BASE_FILES) bh_host_halo_integral_JT.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bh_host_halo_integral_JT
sfr_map:
	$(CC) $(BASE_FILES) gen_sfr_map.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_sfr_map
uvlf:
	$(CC) $(UVLF_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(UVLF_NAME)
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
asymmetry_check:
	$(CC) $(BASE_FILES) gen_asymmetry_check.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_asymmetry_check
erdf:
	$(CC) $(BASE_FILES) gen_erdf.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o gen_erdf
gen_bh_density_split:
	$(CC) $(BASE_FILES) bh_density_split.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bh_density_split
gen_bh_density:
	$(CC) $(BASE_FILES) bh_density.c sm_limits.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o bh_density
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
clean:
	rm -f $(CACHE_OBJ_NAME) $(ASMF_OBJ_NAME) $(SMF_NAME) $(R_NAME) $(SMASS_NAME) $(SM2_NAME) $(FIT_NAME) $(DEC_NAME) $(COSMO_NAME) $(SMASS_CV_NAME) $(FIT2_NAME) $(REMARG_NAME)
	rm -f *~
	rm -f gmon.out *.gcda *.gcno 
	rm:wq
	-rf *.dSYM
