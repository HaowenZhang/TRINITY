CC=gcc
GSLFLAGS=-L$(HOME)/Software/gsl/2.6/lib -I$(HOME)/Software/gsl/2.6/include
CFLAGS=-lm -m64 -lgsl -lgslcblas -Wall -fopenmp $(GSLFLAGS) -L/opt/local/lib -I/opt/local/include -fno-math-errno
OPT_FLAGS=-O3 -msse3
DEBUG_FLAGS=-g -Wall -O0
PROF_FLAGS=-O3 -pg -ftest-coverage -fprofile-arcs
DIR_BIN = ../bin
DIR_BASE=./base
DIR_MCMC=./fitting_mcmc
DIR_OTHER=./others
DIR_PAPER1=./paper1
DIR_PAPER2=./paper2
DIR_PAPER3=./paper3
DIR_PAPER4=./paper4
DIR_TBD=./paperTBD
DIR_DEBUGGING=./debugging
DIR_MOCK_CATALOG=./mock_catalog
BASE_FILES=$(DIR_BASE)/universe_time.c $(DIR_BASE)/observations.c $(DIR_BASE)/mt_rand.c $(DIR_BASE)/jacobi.c $(DIR_BASE)/integrate.c $(DIR_BASE)/distance.c $(DIR_BASE)/expcache2.c $(DIR_BASE)/all_smf.c $(DIR_BASE)/gen_mlists.c $(DIR_BASE)/calc_sfh.c $(DIR_BASE)/smoothing.c $(DIR_BASE)/smloss.c $(DIR_BASE)/check_syscalls.c $(DIR_BASE)/mah.c $(DIR_BASE)/merger_rates.c $(DIR_BASE)/sm_limits.c
SSFR_FILES = $(BASE_FILES) $(DIR_PAPER1)/gen_ssfr.c
SFR_FILES = $(BASE_FILES) gen_sfr_check.c
CPLOT_FILES = $(BASE_FILES) gen_colorplots.c
CSFR_FILES = $(BASE_FILES) $(DIR_PAPER1)/gen_csfr.c
CATALOG_FILES= $(BASE_FILES) sm_catalog.c stringparse.c inthash.c 
AGN_UM_CATALOG_FILES = $(BASE_FILES) $(DIR_MOCK_CATALOG)/sm_agn_catalog_from_um.c $(DIR_BASE)/stringparse.c 
ASMF_C_FILES = $(BASE_FILES) 
SMF_FILES = $(BASE_FILES) $(DIR_PAPER1)/gen_smf.c 
UVLF_FILES = $(BASE_FILES) $(DIR_PAPER1)/gen_uvlf.c
UVLF_TEST_FILES = $(BASE_FILES) $(DIR_DEBUGGING)/gen_uvlf_test.c
QF_FILES = $(BASE_FILES) $(DIR_PAPER1)/gen_qf.c
SMF_CHECK_FILES = $(BASE_FILES) gen_smf_check.c
SMASS_FILES = $(BASE_FILES) smass_density.c
SMASS_COMP_FILES = $(BASE_FILES) smass_comparison_density.c
SMASS_S_FILES = $(BASE_FILES) smass_hi_lo.c
SMASS_P_FILES = $(BASE_FILES) smass_chi2.c
SMASS_CV_FILES = $(BASE_FILES) smass_cv_cosmo.c 
SM2_FILES = $(BASE_FILES) smass_sm_avg.c
FIT_FILES = $(BASE_FILES) $(DIR_MCMC)/fitter.c 
CHI2_FILES = $(BASE_FILES) calc_chi2.c
REMATCH_FILES = $(BASE_FILES)  rematch_smass.c
COSMO_FILES = cosmo_volume.c distance.c

CATALOG_NAME=$(DIR_BIN)/sm_catalog
AGN_UM_CATALOG_NAME=$(DIR_BIN)/sm_agn_catalog_from_um
ASMF_OBJ_NAME=$(DIR_BIN)/all_smf_mcmc
SSFR_NAME=$(DIR_BIN)/gen_ssfr
SFR_NAME=$(DIR_BIN)/gen_sfr_check
CSFR_NAME=$(DIR_BIN)/gen_csfr
SMF_NAME=$(DIR_BIN)/gen_smf
UVLF_NAME=$(DIR_BIN)/gen_uvlf
UVLF_TEST_NAME=$(DIR_BIN)/gen_uvlf_test
SMF_CHECK_NAME=$(DIR_BIN)/gen_smf_check
FIT_NAME=$(DIR_BIN)/fitter
CHI2_NAME=$(DIR_BIN)/calc_chi2
QF_NAME = $(DIR_BIN)/gen_qf

all:
	@make reg EXTRA_FLAGS="$(OPT_FLAGS)"

debug:
	@make reg EXTRA_FLAGS="$(DEBUG_FLAGS)"

prof:
	@make reg EXTRA_FLAGS="$(PROF_FLAGS)"

debug_csfr:
	@make csfr EXTRA_FLAGS="$(DEBUG_FLAGS)"

reg:
	$(CC) $(FIT_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(FIT_NAME)
	$(CC) $(ASMF_C_FILES) $(CFLAGS) $(EXTRA_FLAGS) -o $(ASMF_OBJ_NAME)
	
	$(CC) $(BASE_FILES) $(DIR_PAPER1)/gen_edd_r_color.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(DIR_BIN)/gen_edd_r_color
	$(CC) $(BASE_FILES) $(DIR_PAPER1)/gen_qlf.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(DIR_BIN)/gen_qlf
	$(CC) $(BASE_FILES) $(DIR_PAPER1)/gen_qpdf_eta.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(DIR_BIN)/gen_qpdf_eta
	$(CC) $(BASE_FILES) $(DIR_PAPER1)/bh_density.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(DIR_BIN)/bh_density
	$(CC) $(BASE_FILES) $(DIR_PAPER1)/bhmf.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(DIR_BIN)/bhmf
	$(CC) $(BASE_FILES) $(DIR_PAPER1)/bhmf_typeI.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(DIR_BIN)/bhmf_typeI
	$(CC) $(BASE_FILES) $(DIR_PAPER1)/gen_ssfr_smah.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(DIR_BIN)/gen_ssfr_smah
	$(CC) $(BASE_FILES) $(DIR_PAPER1)/gen_sbhar_smah.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(DIR_BIN)/gen_sbhar_smah
	$(CC) $(SMF_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(SMF_NAME)
	$(CC) $(UVLF_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(UVLF_NAME)
	$(CC) $(QF_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(QF_NAME)
	$(CC) $(CSFR_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(CSFR_NAME)
	$(CC) $(SSFR_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(SSFR_NAME)
	
	$(CC) $(BASE_FILES) $(DIR_PAPER2)/mbh_perc_mstar_lbol_individual.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(DIR_BIN)/mbh_perc_mstar_lbol_individual
	$(CC) $(BASE_FILES) $(DIR_PAPER2)/mbh_perc_mstar_lbol.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(DIR_BIN)/mbh_perc_mstar_lbol


mcmc:
	$(CC) $(FIT_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(FIT_NAME)
	$(CC) $(ASMF_C_FILES) $(CFLAGS) $(EXTRA_FLAGS) -o $(ASMF_OBJ_NAME)

paper1_plots:
	$(CC) $(BASE_FILES) $(DIR_PAPER1)/gen_edd_r_color.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(DIR_BIN)/gen_edd_r_color
	$(CC) $(BASE_FILES) $(DIR_PAPER1)/gen_qlf.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(DIR_BIN)/gen_qlf
	$(CC) $(BASE_FILES) $(DIR_PAPER1)/gen_qpdf_eta.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(DIR_BIN)/gen_qpdf_eta
	$(CC) $(BASE_FILES) $(DIR_PAPER1)/bh_density.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(DIR_BIN)/bh_density
	$(CC) $(BASE_FILES) $(DIR_PAPER1)/bhmf.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(DIR_BIN)/bhmf
	$(CC) $(BASE_FILES) $(DIR_PAPER1)/bhmf_typeI.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(DIR_BIN)/bhmf_typeI
	$(CC) $(BASE_FILES) $(DIR_PAPER1)/gen_ssfr_smah.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(DIR_BIN)/gen_ssfr_smah
	$(CC) $(BASE_FILES) $(DIR_PAPER1)/gen_sbhar_smah.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(DIR_BIN)/gen_sbhar_smah
	$(CC) $(SMF_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(SMF_NAME)
	$(CC) $(UVLF_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(UVLF_NAME)
	$(CC) $(QF_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(QF_NAME)
	$(CC) $(CSFR_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(CSFR_NAME)
	$(CC) $(SSFR_FILES) $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(SSFR_NAME)

paper2_plots:
	$(CC) $(BASE_FILES) $(DIR_PAPER2)/mbh_perc_mstar_lbol_individual.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(DIR_BIN)/mbh_perc_mstar_lbol_individual
	$(CC) $(BASE_FILES) $(DIR_PAPER2)/mbh_perc_mstar_lbol.c $(CFLAGS) $(EXTRA_FLAGS) -DGEN_SMF -o $(DIR_BIN)/mbh_perc_mstar_lbol

clean:
# 	rm -f $(CACHE_OBJ_NAME) $(ASMF_OBJ_NAME) $(SMF_NAME) $(R_NAME) $(SMASS_NAME) $(SM2_NAME) $(FIT_NAME) $(DEC_NAME) $(COSMO_NAME) $(SMASS_CV_NAME) $(FIT2_NAME) $(REMARG_NAME)
	rm -f $(CACHE_OBJ_NAME) $(ASMF_OBJ_NAME) $(SMF_NAME) $(FIT_NAME)

	rm -f *~
	rm -f gmon.out *.gcda *.gcno 
	rm:wq
	-rf *.dSYM
