#ifndef ALL_SMF_H
#define ALL_SMF_H
#include "smf.h"

// The SFR threshold value for the prior against high recent SFRs among massive
// halos. See all_smf_chi2_err() in all_smf.c
#define SFR_CONSTRAINT 1.0
// The width in SFR for the same prior.
#define SFR_CONSTRAINT_WIDTH 0.3
// The minimum scale factor from which the SFRs are determined as "recent".
#define SFR_A_CONSTRAINT 0.8
// The minimum halo mass from which the SFRs are included in the prior.
#define SFR_M_CONSTRAINT 13.5

// The minimum halo mass and scale factor from which
// the SFH will be examined in rising_sfh_penalty(). See calc_sfh.c 
#define RISING_SFH_M_CONSTRAINT 14
#define RISING_SFH_A_CONSTRAINT 0.5

// The parameters for the prior against too high M_ICL/M_*
// ratios. See recent_Micl_Mstar_ratio_in_massive_halos() in calc_sfh.c
// and all_smf_chi2_err() in all_smf.c.
#define ICL_RATIO_CONSTRAINT 0.0
#define ICL_RATIO_CONSTRAINT_WIDTH 0.3
#define ICL_RATIO_A_CONSTRAINT 0.8
#define ICL_RATIO_M_CONSTRAINT 15

// Constraint on the fraction of massive halos that are kinetically powerful.
// Due to the lack of good observational constraints, this is now deprecated.
#define KIN_POWER_CONSTRAINT_LOW 44
#define KIN_FRAC_CONSTRAINT_LOW 0.25
// #define KIN_POWER_CONSTRAINT_HIGH 45
#define KIN_POWER_CONSTRAINT_WIDTH 0.02
#define KIN_POWER_A_CONSTRAINT_HIGH 0.8
#define KIN_POWER_A_CONSTRAINT_LOW 0.5
#define KIN_POWER_M_CONSTRAINT 14

// The parameters for the prior against too high radiative AGN powers among
// massive halos at low redshift. See recent_radiative_power_in_massive_halos()
// in calc_sfh.c and all_smf_chi2_err() in all_smf.c.
#define RAD_POWER_CONSTRAINT_LOW 41.00
#define RAD_POWER_CONSTRAINT_HIGH 43
#define RAD_POWER_CONSTRAINT_WIDTH 0.3
#define RAD_POWER_A_CONSTRAINT_HIGH 0.8
#define RAD_POWER_A_CONSTRAINT_LOW 0.5
#define RAD_POWER_M_CONSTRAINT 14

// #define SMF_CPUS 8

// The number of steps in the burn-in phase.
#define BURN_IN_LENGTH (1<<16)

// The number of parameters in Trinity. Note that this number may be different
// than is shown in the papers, because lots of the parameters are fixed or
// even deprecated in the code. Please refer to Zhang et al. (2021) for the
// accurate number.
#define NUM_PARAMS 70

// The initial eigenvalues in each dimension.  See all_smf.c.
#define EFF_0_STEP    0.001
#define EFF_0_A_STEP  0.015
#define EFF_0_A2_STEP  0.015
#define EFF_0_A3_STEP  0.005

#define EXPSCALE_STEP 0.015

#define M_1_STEP     0.001
#define M_1_A_STEP   0.015
#define M_1_A2_STEP   0.015
#define M_1_A3_STEP   0.01

#define DELTA_STEP   0.001
#define DELTA_A_STEP 0.01
#define DELTA_A2_STEP 0.005

#define BETA_STEP    0.01
#define BETA_A_STEP  0.01
#define BETA_A2_STEP 0.001

#define ALPHA_STEP    0.001
#define ALPHA_A_STEP  0.01
#define ALPHA_A2_STEP 0.005
#define ALPHA_A3_STEP 0.01

#define GAMMA_STEP   0.005
#define GAMMA_A_STEP 0.005
#define GAMMA_A2_STEP 0.005

#define LAMBDA_STEP   0.005
#define LAMBDA_A_STEP 0.005
#define LAMBDA_A2_STEP 0.005

#define MU_STEP      0.005
#define KAPPA_STEP   0.005
#define SCATTER_STEP 0.005
#define SIGMA_Z_STEP 0.005

#define ICL_STEP        0.01
#define MU_A_STEP       0.01
#define KAPPA_A_STEP    0.01
#define SCATTER_A_STEP  0.01
#define ICL_FRAC_E_STEP 0.05

#define BURST_DUST_STEP 0.01
#define BURST_DUST_AMP_STEP 0.01
#define BURST_DUST_Z_STEP 0.01
#define RHO_05_STEP 0.003

#define BH_BETA_0_STEP 0.01
#define BH_BETA_1_STEP 0.01
#define BH_BETA_2_STEP 0.01
#define BH_GAMMA_0_STEP 0.01
#define BH_GAMMA_1_STEP 0.01
#define BH_GAMMA_2_STEP 0.01
#define BH_SIGMA_0_STEP 0.01
#define BH_SIGMA_1_STEP 0.01

#define BH_ALPHA_0_STEP 0.01
#define BH_ALPHA_1_STEP 0.01
#define BH_ALPHA_2_STEP 0.01
#define BH_DELTA_0_STEP 0.01
#define BH_DELTA_1_STEP 0.01
#define BH_DELTA_2_STEP 0.01

#define BH_ETA_0_STEP 0.01
#define BH_ETA_1_STEP 0.01
#define BH_DUTY_0_STEP 0.01
#define BH_DUTY_1_STEP 0.01

#define BH_MERGE_F_0_STEP 0.05
#define BH_MERGE_F_1_STEP 0.05
#define BH_MERGE_W_0_STEP 0.05
#define BH_MERGE_W_1_STEP 0.05

#define BH_EFFICIENCY_0_STEP 0.03
#define BH_EFFICIENCY_1_STEP 0.01
#define BH_SCATTER_0_STEP 0.001
#define BH_SCATTER_1_STEP 0.001
#define RHO_BH_STEP 0.01
// #define BH_NORM_STEP 0.01

#define CTK_ALPHA_0_STEP 0.001
#define CTK_ALPHA_1_STEP 0.001
#define CTK_ALPHA_2_STEP 0.001

#define DC_MBH_0_STEP 0.05
#define DC_MBH_1_STEP 0.05
#define DC_MBH_W_0_STEP 0.05
#define DC_MBH_W_1_STEP 0.05

#define ETA_MU_0_STEP 0.001
#define ETA_MU_1_STEP 0.001

// Helper function to get model parameters in a human-friendly way.
#define EFF_0(x)     ((x).params[0])
#define EFF_0_A(x)   ((x).params[1])
#define EFF_0_A2(x)  ((x).params[2])
#define EFF_0_A3(x)  ((x).params[3])
//#define EXPSCALE(x)  ((x).params[4])
#define DC_MBH_1(x) ((x).params[4])
//#define Z_CEIL(x) ((x).params[4])
#define M_1(x)       ((x).params[5])
#define M_1_A(x)     ((x).params[6])
#define M_1_A2(x)    ((x).params[7])
// M_1_A3: an additional parameter compared to the SMHM parametrization.
#define M_1_A3(x)    ((x).params[8])
#define ALPHA(x)      ((x).params[9])
#define ALPHA_A(x)   ((x).params[10])
#define ALPHA_A2(x) ((x).params[11])
// ALPHA_A3: another additional parameter compared to the SMHM parametrization.
#define ALPHA_A3(x) ((x).params[12])
#define BETA(x)      ((x).params[13])
#define BETA_A(x)   ((x).params[14])
#define BETA_A2(x) ((x).params[15])
#define SIGMA_A(x) ((x).params[24])
#define DELTA(x)     ((x).params[16])
#define DELTA_A(x)   ((x).params[17])
#define DELTA_A2(x)  ((x).params[18])
#define GAMMA(x)    ((x).params[19])
#define GAMMA_A(x)  ((x).params[20])
#define GAMMA_A2(x)  ((x).params[21])
// #define LAMBDA(x)    ((x).params[22])
// #define LAMBDA_A(x)  ((x).params[23])

#define BH_DUTY_M_0(x)   ((x).params[22])
#define BH_DUTY_M_1(x)  ((x).params[23])

//#define LAMBDA_A2(x)  ((x).params[24])
#define DC_MBH_W_1(x) ((x).params[24])
#define MU(x)       ((x).params[25])
#define KAPPA(x)    ((x).params[26])
#define SCATTER(x)  ((x).params[27])
#define SIGMA_Z(x)  ((x).params[28])
#define ICL_FRAC(x) ((x).params[29])
#define ICL_FRAC_Z(x) ((x).params[17])
#define MU_A(x)     ((x).params[30])
#define KAPPA_A(x)  ((x).params[31])
#define SCATTER_A(x) ((x).params[32])
#define ICL_FRAC_E(x) ((x).params[33])
// #define BURST_DUST(x) ((x).params[34])
// #define BURST_DUST_AMP(x) ((x).params[35])
// #define BURST_DUST_Z(x) ((x).params[36])
#define RHO_BH_0(x) ((x).params[34])
#define RHO_BH_1(x) ((x).params[35])
#define RHO_BH_2(x) ((x).params[36])
#define RHO_05(x) ((x).params[37])



#define BH_BETA_0(x)       ((x).params[38])
#define BH_BETA_1(x)       ((x).params[39])
#define BH_BETA_2(x)       ((x).params[40])
#define BH_GAMMA_0(x)      ((x).params[41])
#define BH_GAMMA_1(x)      ((x).params[42])
#define BH_GAMMA_2(x)      ((x).params[43])
#define BH_ALPHA_0(x)      ((x).params[44])
#define BH_ALPHA_1(x)      ((x).params[45])
#define BH_DELTA_0(x)      ((x).params[46])
#define BH_DUTY_ALPHA_0(x)       ((x).params[47])
#define BH_DUTY_ALPHA_1(x)       ((x).params[48])
#define BH_MERGE_F_0(x)    ((x).params[49])
#define BH_MERGE_F_1(x)    ((x).params[50])
// #define BH_MERGE_W_0(x)    ((x).params[51])
#define F_OCC_MIN_0(x)    ((x).params[51])
// #define BH_MERGE_W_1(x)    ((x).params[52])
#define F_OCC_MIN_1(x)    ((x).params[52])
// #define BH_EFFICIENCY_0(x) ((x).params[53])
// #define BH_EFFICIENCY_1(x) ((x).params[54])
#define BH_ETA_CRIT_0(x) ((x).params[53])
#define BH_ETA_CRIT_1(x) ((x).params[54])
#define BH_SCATTER_0(x)    ((x).params[55])
#define BH_SCATTER_1(x)    ((x).params[56])
// #define BH_NORM(x)         ((x).params[57])
// #define RHO_BH(x)         ((x).params[57])
#define ABHMF_SHIFT(x)         ((x).params[57])
#define DC_MBH_0(x)           ((x).params[58])
#define DC_MBH_W_0(x)           ((x).params[59])
#define ETA_MU_0(x)     ((x).params[60])
// #define ETA_MU_1(x)     ((x).params[61])
#define BH_EFFICIENCY_0(x)     ((x).params[61])
#define BH_DELTA_1(x)      ((x).params[62])


#define QM_0(x) ((x).params[63])
#define QM_1(x) ((x).params[64])
#define QM_2(x) ((x).params[65])
#define QWIDTH_0(x) ((x).params[66])
#define QWIDTH_1(x) ((x).params[67])
#define QWIDTH_2(x) ((x).params[68])


#define MODEL_FLAGS(x)     ((x).params[69]) 

// Flag indicating if the model is invalid (!= 0) or not (== 0).
#define INVALID(x)  ((x).params[NUM_PARAMS])
#define CHI2(x)     ((x).params[NUM_PARAMS+1])
#ifdef __APPLE__
#define INVALIDATE(x,y) { INVALID(*x) = 1;  fprintf(stderr, "%s\n", y); }
#else
#define INVALIDATE(x,y) { INVALID(*x) = 1; }
#endif

// The structure containing all the model parameters, 
// including the model flag, chi2, and the invalidity flag.
struct smf_fit 
{
  double params[NUM_PARAMS+2];
};

void shutdown_clients(void); //Deprecated.

// For the non-deprecated function below, please see
// the comments in all_smf.c.
void init_frac_below8(void);
float all_smf_chi2_err(struct smf_fit test);
float all_smf_chi2_err_write(struct smf_fit test); //Deprecated.
float chi2_type(struct smf_fit test);
void clear_stats(void);
void add_to_stats(struct smf_fit *a);
void stats_to_step(int64_t num_stats);
void random_step(struct smf_fit *a);
void init_orth_matrix(void);
void print_orth_matrix(int type, int syst, float scatter, float z);
void init_mcmc_from_args(int argc, char **argv);
int mcmc_smf_fit(int length, int output);
void set_identity_portion(void);
//void read_params(char *buffer, double *data, int max_n);
struct smf single_smf_cb(float z, void *extra_data);
struct smf linear_smf_cb(float z, void *extra_data);
void assert_model(struct smf_fit *a);

extern int nonlinear_luminosity, vel_dispersion;

#endif /* ALL_SMF_H */
