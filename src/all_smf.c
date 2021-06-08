#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <sys/time.h>
#include <sys/wait.h>
#include "smf.h"
#include "all_smf.h"
#include "observations.h"
#include "jacobi.h"
#include "mt_rand.h"
#include "calc_sfh.h"
#include "mlist.h"
#include "psf_cache.h"

struct timestep *steps = NULL;
int64_t num_outputs = 0;
int64_t num_points = 0;
struct obs_smf_point *obs_smfs = NULL;
int64_t num_obs_smfs = 0;
int no_matching_scatter = 0;
int no_obs_scatter = 0;
int no_systematics = 0;
int nonlinear_luminosity = 1;
int no_bh_mergers = 0;
int vel_dispersion = 0;
int z_scaling_only = 0;
int no_z_scaling = 0;
int no_sfr_constraint = 0;
int dplaw_only = 0;



double bh_scatter_obs = 0.4;
double area_sdss = 14500.0 / 41252.96124941928; //fractional survey area of sdss, ***relative to*** the whole sky.
double frac_below8[MBH_BINS] = {0};
double frac_below11[MBH_BINS] = {0};

// Calculate the frac_below8 and frac_below11, i.e., 
// the fractional probability that BH in each mass bin
// get scattered below 10^8 and 10^11 Msun, respectively
void init_frac_below8(void)
{
  for (int i=0; i<MBH_BINS; i++)
  {
    double mbh = MBH_MIN + (i + 0.5) * MBH_INV_BPDEX;
    frac_below8[i] = 0.5 * (1 + erf((8.0 - mbh) / bh_scatter_obs * M_SQRT1_2));
    frac_below11[i] = 0.5 * (1 + erf((11.0 - mbh) / bh_scatter_obs * M_SQRT1_2));
  }
}



float z_limit = 0;
float z_max = 0;
float z_min = 100;
float inv_temperature = 1;
struct smf_fit initial_smf_fit, smf_sum;
double identity_portion = 3e-5;
// Covariance matrix for all the model parameters.
double cov_matrix[NUM_PARAMS][NUM_PARAMS];
// orth_matrix is the correlation coefficient matrix for all parameters.
double orth_matrix[NUM_PARAMS][NUM_PARAMS];

// Eigenvalue vector. This eigenvalue vector determines the Gaussian
// width in the proposal of new samples.
double eigenvalues[NUM_PARAMS] = { EFF_0_STEP, EFF_0_A_STEP, EFF_0_A2_STEP, 
				   EFF_0_A3_STEP, EXPSCALE_STEP, M_1_STEP,
				   M_1_A_STEP, M_1_A2_STEP, M_1_A3_STEP,
				   ALPHA_STEP, ALPHA_A_STEP, ALPHA_A2_STEP,
				   ALPHA_A3_STEP, 
				   BETA_STEP, BETA_A_STEP, BETA_A2_STEP,
				   DELTA_STEP, DELTA_A_STEP, DELTA_A2_STEP,
				   GAMMA_STEP, GAMMA_A_STEP, GAMMA_A2_STEP,
				   LAMBDA_STEP, LAMBDA_A_STEP, LAMBDA_A2_STEP,
				   MU_STEP, KAPPA_STEP,
				   SCATTER_STEP, SIGMA_Z_STEP, ICL_STEP,
				   MU_A_STEP, KAPPA_A_STEP, SCATTER_A_STEP,
				   ICL_FRAC_E_STEP, BURST_DUST_STEP,
				   BURST_DUST_AMP_STEP, BURST_DUST_Z_STEP,
				   RHO_05_STEP,
           BH_BETA_0_STEP, BH_BETA_1_STEP, BH_BETA_2_STEP,
           BH_GAMMA_0_STEP, BH_GAMMA_1_STEP, BH_GAMMA_2_STEP,
           BH_ALPHA_0_STEP, BH_ALPHA_1_STEP, BH_DELTA_0_STEP,
           BH_DUTY_0_STEP, BH_DUTY_1_STEP,
           BH_MERGE_F_0_STEP, BH_MERGE_F_1_STEP, 
           BH_MERGE_W_0_STEP, BH_MERGE_W_1_STEP, 
           BH_EFFICIENCY_0_STEP,
           BH_EFFICIENCY_1_STEP,
           BH_SCATTER_0_STEP, BH_SCATTER_1_STEP,
           RHO_BH_STEP, 
           DC_MBH_0_STEP, DC_MBH_W_0_STEP,
          ETA_MU_0_STEP, ETA_MU_1_STEP, BH_DELTA_1_STEP, 
          BH_BETA_0_STEP, BH_BETA_1_STEP, BH_BETA_2_STEP,
           BH_GAMMA_0_STEP, BH_GAMMA_1_STEP, BH_GAMMA_2_STEP,
           1e-9 };

// A backup eigenvalue vector, in case we need to use them
// when the vector above has been updated.
double eigenvalues_original[NUM_PARAMS] = { EFF_0_STEP, EFF_0_A_STEP, EFF_0_A2_STEP, 
           EFF_0_A3_STEP, EXPSCALE_STEP, M_1_STEP,
           M_1_A_STEP, M_1_A2_STEP, M_1_A3_STEP,
           ALPHA_STEP, ALPHA_A_STEP, ALPHA_A2_STEP,
           ALPHA_A3_STEP, 
           BETA_STEP, BETA_A_STEP, BETA_A2_STEP,
           DELTA_STEP, DELTA_A_STEP, DELTA_A2_STEP,
           GAMMA_STEP, GAMMA_A_STEP, GAMMA_A2_STEP,
           LAMBDA_STEP, LAMBDA_A_STEP, LAMBDA_A2_STEP,
           MU_STEP, KAPPA_STEP,
           SCATTER_STEP, SIGMA_Z_STEP, ICL_STEP,
           MU_A_STEP, KAPPA_A_STEP, SCATTER_A_STEP,
           ICL_FRAC_E_STEP, BURST_DUST_STEP,
           BURST_DUST_AMP_STEP, BURST_DUST_Z_STEP,
           RHO_05_STEP,
           BH_BETA_0_STEP, BH_BETA_1_STEP, BH_BETA_2_STEP,
           BH_GAMMA_0_STEP, BH_GAMMA_1_STEP, BH_GAMMA_2_STEP,
           BH_ALPHA_0_STEP, BH_ALPHA_1_STEP, BH_DELTA_0_STEP,
           BH_DUTY_0_STEP, BH_DUTY_1_STEP,
           BH_MERGE_F_0_STEP, BH_MERGE_F_1_STEP, 
           BH_MERGE_W_0_STEP, BH_MERGE_W_1_STEP, 
           BH_EFFICIENCY_0_STEP,
           BH_EFFICIENCY_1_STEP,
           BH_SCATTER_0_STEP, BH_SCATTER_1_STEP,
           RHO_BH_STEP, 
           DC_MBH_0_STEP, DC_MBH_W_0_STEP,
          ETA_MU_0_STEP, ETA_MU_1_STEP, BH_DELTA_1_STEP, 
          BH_BETA_0_STEP, BH_BETA_1_STEP, BH_BETA_2_STEP,
           BH_GAMMA_0_STEP, BH_GAMMA_1_STEP, BH_GAMMA_2_STEP,
           1e-9 };

// Read model parameters from buffer and store them
// in data. Repeat for at most max_n times.
// static void read_params(char *buffer, double *data, int max_n) 
void read_params(char *buffer, double *data, int max_n) 
{
  int num_entries = 0;
  char *cur_pos = buffer, *end_pos;
  float val = strtod(cur_pos, &end_pos);
  while (cur_pos != end_pos && num_entries < max_n) 
  {
    data[num_entries] = val;
    num_entries++;
    cur_pos=end_pos;
    while (*cur_pos==' ' || *cur_pos=='\t' || *cur_pos=='\n') cur_pos++;
    val = strtod(cur_pos, &end_pos);
  }
}

// Read the initial model parameters and calculate the initial chi2.
void initial_conditions(void) 
{
  char buffer[2048];
  fgets(buffer, 2048, stdin);
  read_params(buffer, initial_smf_fit.params, NUM_PARAMS);

  // set some parameters to predefined values if needed.
  assert_model(&initial_smf_fit);

  fprintf(stderr, "Initial params:\n");
  int i;
  for (i=0; i<NUM_PARAMS-1; i++) fprintf(stderr, "%.12f ", initial_smf_fit.params[i]);
  printf("%.12f\n", initial_smf_fit.params[NUM_PARAMS-1]);

  CHI2(initial_smf_fit) = all_smf_chi2_err(initial_smf_fit);

  // Calculate the chi2's contributed by each type of observational data.
  chi2_type(initial_smf_fit);
  // Initialize the orth_matrix that specifies how well the model parameters
  // are correlated to each other.
  init_orth_matrix();
  fprintf(stderr, "Initial Chi2 fit: %f\n", CHI2(initial_smf_fit));
}

// The function to calculate chi2 given the parameter set.
float all_smf_chi2_err(struct smf_fit test) 
{
  float chi2 = 0;
  double chi2_bhmf = 0;
  int i;

  // Two limiting cases for the smhm model.
  z_limit = 8.5;
  struct smf smf_limit = smhm_at_z(z_limit, test);
  struct smf smf_mid = smhm_at_z(2, test);
  // struct smf smf_z4 = smhm_at_z(4.0, test);
  INVALID(test) = 0;
  if (no_z_scaling) 
  {
    if (M_1(test) < 9.5) return -1;
    if (ALPHA(test) < -4) return -1;
    if (GAMMA(test) > 1.05) return -1;
  }

  // exclude model parameters that lie outside of physical ranges.
  if (fabs(KAPPA(test))>KAPPA_PRIOR*2.5 || 
      ALPHA(test) >= -1 || smf_limit.alpha >= -1 || smf_mid.alpha >= -1 ||
      ALPHA(test) <= -15 || smf_limit.alpha <= -15 || smf_mid.alpha <= -15 ||
      BETA(test) < ALPHA(test) || smf_limit.beta < smf_limit.alpha || smf_mid.beta < smf_mid.alpha ||
      BETA(test) > 15 || smf_limit.beta > 15 || smf_mid.beta > 15 ||
      DELTA(test) < 0.01 || smf_limit.delta < 0.01 || smf_mid.delta < 0.01 ||
      fabs(DELTA(test)) > 10 || fabs(smf_limit.delta) > 10 || fabs(smf_mid.delta) > 10 ||
      GAMMA(test) < -5 || 
      SIGMA_Z(test) < 0 || EFF_0(test) > 1 ||
      SIGMA_Z(test) < 0 || EFF_0(test) > 3 ||
      EFF_0_A(test) < -2 ||
      smf_limit.epsilon > 1e5 || smf_mid.epsilon > 1e4 ||
      SIGMA_A(test) > 0 ||
      SCATTER(test) > 0.3 ||
      SCATTER(test) < 0 || smf_limit.obs_scatter > 1.0 ||
      BURST_DUST_Z(test) < 0.8 || BURST_DUST_Z(test)>6 ||
      BURST_DUST(test) < 0 || BURST_DUST(test) > 1 ||
      BURST_DUST_AMP(test) < 0 ||
      BURST_DUST_AMP(test) > 5 || fabs(MU_A(test)) > 0.8 ||
      fabs(KAPPA_A(test)) > KAPPA_LIMIT || SCATTER_A(test) > SCATTER(test) ||
      RHO_05(test) < 0.23 || RHO_05(test) > 1.0 ||

      smf_limit.icl_frac < 1e-20 || smf_limit.icl_frac > 1 ||
      ICL_FRAC(test) < -20 || ICL_FRAC(test) > 0 ||
      QM_0(test) < 0 || QM_0(test) > 8 ||
      smf_limit.qm < 0 ||
      QWIDTH_0(test) < 0.001 || QWIDTH_0(test) > 2 ||
      smf_limit.qwidth < 0.001 || smf_limit.qwidth > 2 ||

      smf_limit.bh_beta < 4 || smf_limit.bh_beta > 12 ||
      smf_limit.bh_gamma < 0.5 || smf_limit.bh_gamma > 8.5 ||
      BH_ALPHA_0(test) < -2.0 || BH_ALPHA_0(test) >= 1.0 || 
      BH_DELTA_0(test) < -2.0 || BH_DELTA_0(test) > 3.0 ||
      BH_ALPHA_0(test) > BH_DELTA_0(test) ||
      smf_limit.bh_alpha > smf_limit.bh_delta || 
      fabs(smf_limit.bh_alpha) > 4 ||
      fabs(smf_limit.bh_alpha) >= 1.0 ||
      fabs(smf_limit.bh_delta) > 4 ||

      BH_MERGE_F_0(test) < -5 || BH_MERGE_F_0(test) > 5 ||
      smf_limit.f_merge_bh < 1e-5 || smf_limit.f_merge_bh > 1e5 || 
      BH_ETA_CRIT_0(test) < -2.5 || BH_ETA_CRIT_0(test) > -0.3 ||
      smf_limit.bh_eta_crit < -2.5 || smf_limit.bh_eta_crit > -0.3 ||
      log10(smf_limit.bh_efficiency_rad) < -2.5 || log10(smf_limit.bh_efficiency_rad) > -0.3 ||
      BH_DUTY_0(test) < 0.001 || BH_DUTY_0(test) > 1 ||
      BH_DUTY_1(test) > BH_DUTY_0(test) ||
      smf_limit.bh_duty < 0.001 || smf_limit.bh_duty > 1 ||
      BH_SCATTER_0(test) < 0 || BH_SCATTER_0(test) > 2 ||
      smf_limit.bh_scatter < 0 || smf_limit.bh_scatter > 2 ||
      DC_MBH_0(test) < 4 || DC_MBH_0(test) > 10 ||
      DC_MBH_W_0(test) < 0.001 || DC_MBH_W_0(test) > 4 ||
      ABHMF_SHIFT(test) < -1 || ABHMF_SHIFT(test) > 1)

  {
    //fprintf(stderr, "Is NOT OK!\n");
  		return (-1);
  }
  #pragma omp parallel
    {
      calc_sfh(&test);

      // The number density of massive BHs should not have decreased from z=1 to z=0.
      double bhz0 = calc_bhmf(BH_MASS_TO_REQUIRE_ND_GROWTH, 0);
      double bhz1 = calc_bhmf(BH_MASS_TO_REQUIRE_ND_GROWTH, 1);
      if ((bhz1 < 0) || (bhz0 < 0) || (bhz0 < bhz1)) 
      {  
        //fprintf(stderr, "BH ND down. bhz0: %.6e, bhz1: %.6e\n", bhz0, bhz1);
        INVALIDATE(&test, "BH Number Density Decreased from z=1 to z=0!");
      }

  #pragma omp for schedule(simd:dynamic, 8) //private(m,epsilon,smf_val)
      for (i=0; i<num_obs_smfs; i++)
        if (!(INVALID(test)))
        {
          obs_smfs[i].chi2 = calc_single_chi2_err_point(&(obs_smfs[i]));          
        }
  	  
    }

  if (INVALID(test)) 
  {
    for (i = 0; i < num_outputs; i++)
    {
      if (steps[i].flag_alloc)
      {
        gsl_spline_free(steps[i].spline); //free the spline object since we are allocating it every time.
        gsl_spline_free(steps[i].spline_sfr);
        if (steps[i].alloc2smf)
        {
          gsl_spline_free(steps[i].spline2); //same for the second segment of the SMHM relation.
          gsl_spline_free(steps[i].spline_sfr2);
        }
      }

      if (steps[i].flag_alloc == 2)
      {
        gsl_spline_free(steps[i].spline_uv); //free the spline object since we are allocating it every time.
        if (steps[i].alloc2uvlf) gsl_spline_free(steps[i].spline_uv2); //same for the second segment of the UV-HM relation.
      }
    }
    return -1;
  }

  for (i=0; i<num_obs_smfs; i++) 
  {
    chi2 += obs_smfs[i].chi2;
  }

  //fprintf(stderr, "Before adding priors, chi2=%f\n", chi2);

  // Apply the priors on several model parameters. See the table of priors
  // in Zhang et al. (2021)
  chi2 += pow(MU(test)/MU_PRIOR, 2) + pow(KAPPA(test)/KAPPA_PRIOR, 2)
    + pow(MU_A(test)/MU_PRIOR, 2) + pow(KAPPA_A(test)/KAPPA_PRIOR, 2)
    + pow(SCATTER_A(test)/MU_PRIOR, 2)
    + pow((SCATTER(test)-SCATTER_CENTER)/SCATTER_PRIOR, 2)
    + pow(((-0.5*SIGMA_A(test)+SIGMA_Z(test))-SIGMA_Z_CENTER)/SIGMA_Z_PRIOR, 2);

  if (!vel_dispersion) 
  {
    chi2 += pow((BH_BETA_0(test)-BH_BETA_CENTER)/BH_BETA_PRIOR, 2)
         + pow((BH_GAMMA_0(test)-BH_GAMMA_CENTER)/BH_GAMMA_PRIOR, 2);
  } 
  else 
  {
    chi2 += pow((BH_BETA_0(test)-BH_BETA_CENTER_VD)/BH_BETA_PRIOR, 2)
         + pow((BH_GAMMA_0(test)-BH_GAMMA_CENTER_VD)/BH_GAMMA_PRIOR, 2);
  }

  //fprintf(stderr, "after adding priors on parameters, chi2=%f\n", chi2);

  // Priors on stellar mass--halo mass relations and against a hole in BHARs.
  double chi2_smhm = 0;
  double chi2_bhar = 0;
   for (i = 1; i < num_outputs; i++)
   {
      for (int j=1; j < M_BINS; j++)
      {
        if (steps[i].t[j] >= 2e-8 && steps[i].log_sm[j] > 0 && steps[i].log_sm[j] < steps[i].log_sm[j - 1])
        {
          double delta = steps[i].log_sm[j - 1] - steps[i].log_sm[j];
          chi2_smhm += delta * delta;
        }
    		if (steps[i].scale > 0.5 && steps[i].med_hm_at_a[j] >= 11.5 && steps[i].med_hm_at_a[j] <= 12.5 && steps[i].bh_acc_rate[j] <= steps[i].bh_acc_rate[j-1] && steps[i].bh_acc_rate[j] <= steps[i].bh_acc_rate[j+1])
    		{
    			double dd = 0.5 * (steps[i].bh_acc_rate[j+1] + steps[i].bh_acc_rate[j-1]) / steps[i].bh_acc_rate[j];
    			chi2_bhar += dd * dd;
    		}
      }
   }

   chi2 += chi2_smhm;
   chi2 += chi2_bhar;

   if (!no_z_scaling && !no_sfr_constraint) 
   {
    double chi2_sfr = 0;
    double chi2_icl = 0;
    double chi2_kin = 0;
    double chi2_rad = 0;
    double chi2_qso = 0;
    
    double penalty_rising_sfh = rising_sfh_penalty();
    // fprintf(stderr, "The chi2 penalty for rising SFH in massive halos at low-z: %f\n", penalty_rising_sfh); 

    // priors on high-z bright quasars
    double ratio_z6_qso = ratio_high_z_qso(5.7, 6.5, 47, area_sdss);
    fprintf(stderr, "The ratio of quasars with Lbol>10^47, Mbh<10^8 vs. Mbh>10^8, and 5.7 < z < 6.5 is %e\n", ratio_z6_qso);

    double n_z6_qso = number_high_z_low_mass_qso(5.7, 6.5, 47, area_sdss);
    fprintf(stderr, "The number of quasars with Lbol>10^47, Mbh<10^8, and 5.7 < z < 6.5 is %e\n", n_z6_qso);

    if (n_z6_qso > 0) chi2_qso = 2 * n_z6_qso;
    chi2_qso += log10(ratio_z6_qso) * 10;

    // Prior on recent SFR histories among massive halos
    double recent_sfr = recent_sfh_in_massive_halos();

    if (recent_sfr > SFR_CONSTRAINT)
      chi2_sfr = pow((recent_sfr-SFR_CONSTRAINT)/(SFR_CONSTRAINT_WIDTH), 2);

    // Prior on intra-cluster light to stellar mass ratios
    double recent_icl_star_ratio = recent_Micl_Mstar_ratio_in_massive_halos();
    // Prior on recent (low redshift) radiative luminosities among massive halos.
    double recent_radiative_power = recent_radiative_power_in_massive_halos();

    if (recent_icl_star_ratio > ICL_RATIO_CONSTRAINT)
      chi2_icl = pow((recent_icl_star_ratio - ICL_RATIO_CONSTRAINT) / (ICL_RATIO_CONSTRAINT_WIDTH), 2);

    if (recent_radiative_power > RAD_POWER_CONSTRAINT_HIGH)
      chi2_rad = pow((recent_radiative_power - RAD_POWER_CONSTRAINT_HIGH) / (RAD_POWER_CONSTRAINT_WIDTH), 2);


    chi2 += chi2_sfr + chi2_icl + chi2_qso + chi2_rad;
    chi2 += penalty_rising_sfh;
    fprintf(stderr, "recent rad power: %f, chi2_rad:%f\n", recent_radiative_power, chi2_rad);
  }


  // Note that since we re-allocate space for these spline objects every time due to
  // their varying sizes, we should free them every time after using them to calculate
  // the observables, to avoid out-of-memory kills.
  for (i = 0; i < num_outputs; i++)
  {
    if (steps[i].flag_alloc)
    {
      gsl_spline_free(steps[i].spline); //free the spline object since we are allocating it every time.
      gsl_spline_free(steps[i].spline_sfr);
      if (steps[i].alloc2smf)
      {
        gsl_spline_free(steps[i].spline2); //same for the second segment of the SMHM relation.
        gsl_spline_free(steps[i].spline_sfr2);
      }
    }

    if (steps[i].flag_alloc == 2)
    {
      gsl_spline_free(steps[i].spline_uv); //free the spline object since we are allocating it every time.
      if (steps[i].alloc2uvlf) gsl_spline_free(steps[i].spline_uv2); //same for the second segment of the UV-HM relation.
    }
  }
  return (chi2);
}

// Calculate the chi2's for each different type of observations.
// This function is very similar to all_smf_chi2_err(), except 
// that we gather chi2's for different types of data.
float chi2_type(struct smf_fit test) 
{
  float chi2 = 0;
  double chi2_bhmf = 0;
  int i;
  z_limit = 8.5;
  struct smf smf_limit = smhm_at_z(z_limit, test);
  struct smf smf_mid = smhm_at_z(2, test);

  INVALID(test) = 0;
  if (no_z_scaling) 
  {
    if (M_1(test) < 9.5) return -1;
    if (ALPHA(test) < -4) return -1;
    if (GAMMA(test) > 1.05) return -1;
  }


  if (fabs(KAPPA(test))>KAPPA_PRIOR*2.5 || 
      ALPHA(test) >= -1 || smf_limit.alpha >= -1 || smf_mid.alpha >= -1 ||
      ALPHA(test) <= -15 || smf_limit.alpha <= -15 || smf_mid.alpha <= -15 ||
      BETA(test) < ALPHA(test) || smf_limit.beta < smf_limit.alpha || smf_mid.beta < smf_mid.alpha ||
      BETA(test) > 15 || smf_limit.beta > 15 || smf_mid.beta > 15 ||
      DELTA(test) < 0.01 || smf_limit.delta < 0.01 || smf_mid.delta < 0.01 ||
      fabs(DELTA(test)) > 10 || fabs(smf_limit.delta) > 10 || fabs(smf_mid.delta) > 10 ||
      LAMBDA(test) > 10 || smf_limit.lambda > 10 || 
      LAMBDA(test) < -10 || smf_limit.lambda < -10 || 
      SIGMA_Z(test) < 0 || EFF_0(test) > 3 ||
      smf_limit.epsilon > 1e5 || smf_mid.epsilon > 1e4 ||
      SIGMA_A(test) > 0 ||
      SCATTER(test) > 0.3 ||
      SCATTER(test) < 0 || smf_limit.obs_scatter > 1.0 ||
      BURST_DUST_Z(test) < 0.8 || BURST_DUST_Z(test)>6 ||
      BURST_DUST(test) < 0 || BURST_DUST(test) > 1 ||
      BURST_DUST_AMP(test) < 0 ||
      BURST_DUST_AMP(test) > 5 || fabs(MU_A(test)) > 0.8 ||
      fabs(KAPPA_A(test)) > KAPPA_LIMIT || SCATTER_A(test) > SCATTER(test) ||
      RHO_05(test) < 0.23 || RHO_05(test) > 1.0 ||

      smf_limit.icl_frac < 1e-20 || smf_limit.icl_frac > 1 ||
      ICL_FRAC(test) < -20 || ICL_FRAC(test) > 0 ||
      QM_0(test) < 9 || QM_0(test) > 14 ||
      smf_limit.qm < 9 ||
      QWIDTH_0(test) < 0.001 || QWIDTH_0(test) > 2 ||
      smf_limit.qwidth < 0.001 || smf_limit.qwidth > 2||

      smf_limit.bh_beta < 4 || smf_limit.bh_beta > 12 ||
      smf_limit.bh_gamma < 0.5 || smf_limit.bh_gamma > 8.5 ||
      BH_ALPHA_0(test) < -0.5 || BH_ALPHA_0(test) > 3.0 || 
      BH_DELTA_0(test) < -0.5 || BH_DELTA_0(test) > 3.0 ||
      BH_ALPHA_0(test) > BH_DELTA_0(test) ||
      smf_limit.bh_alpha > smf_limit.bh_delta || 
      fabs(smf_limit.bh_alpha) > 4 ||
      
      BH_ETA_CRIT_0(test) < -2.5 || BH_ETA_CRIT_0(test) > -0.3 ||
      smf_limit.bh_eta_crit < -2.5 || smf_limit.bh_eta_crit > -0.3 ||
      log10(smf_limit.bh_efficiency_rad) < -2.5 || log10(smf_limit.bh_efficiency_rad) > -0.3 ||
      BH_DUTY_0(test) < 0.001 || BH_DUTY_0(test) > 1 ||
      BH_DUTY_1(test) > BH_DUTY_0(test) ||
      smf_limit.bh_duty < 0.001 || smf_limit.bh_duty > 1 ||
      BH_SCATTER_0(test) < 0 || BH_SCATTER_0(test) > 0.4 ||
      smf_limit.bh_scatter < 0 || smf_limit.bh_scatter > 2 ||
      DC_MBH_0(test) < 4 || DC_MBH_0(test) > 10 ||
      DC_MBH_W_0(test) < 0.001 || DC_MBH_W_0(test) > 4 ||
      ABHMF_SHIFT(test) < -1 || ABHMF_SHIFT(test) > 1)
      
      {
        fprintf(stderr, "Is NOT OK!\n");
      }
  #pragma omp parallel
    {
      calc_sfh(&test);      
  #pragma omp for schedule(simd:dynamic, 8) //private(m,epsilon,smf_val)
      for (i=0; i<num_obs_smfs; i++)
        // if (!(INVALID(test)))
        {
          obs_smfs[i].chi2 = calc_single_chi2_err_point(&(obs_smfs[i]));
        }
      
    }

  // the array to store chi2's of different observations.
  double chi2_type[11] = {0};

  for (i=0; i<num_obs_smfs; i++) 
  {
    chi2 += obs_smfs[i].chi2;
    // Add up the contribution from each type.
    chi2_type[obs_smfs[i].type] += obs_smfs[i].chi2;
  }

  for (i=0; i < 11; i++)
    fprintf(stderr, "The total chi2 for the %d-th type data is: %f\n", i, chi2_type[i]);


  if (INVALID(test)) 
  {
    for (i = 0; i < num_outputs; i++)
    {
      if (steps[i].flag_alloc)
      {
        gsl_spline_free(steps[i].spline); //free the spline object since we are allocating it every time.
        gsl_spline_free(steps[i].spline_sfr);
        if (steps[i].alloc2smf)
        {
          gsl_spline_free(steps[i].spline2); //same for the second segment of the SMHM relation.
          gsl_spline_free(steps[i].spline_sfr2);
        }
      }

      if (steps[i].flag_alloc == 2)
      {
        gsl_spline_free(steps[i].spline_uv); //free the spline object since we are allocating it every time.
        if (steps[i].alloc2uvlf) gsl_spline_free(steps[i].spline_uv2); //same for the second segment of the UV-HM relation.
        // gsl_spline_free(stseps[i].spline_std_uv);
      }
    }
    return -1;
  }

  for (i=0; i<num_obs_smfs; i++) 
  {
    chi2 += obs_smfs[i].chi2;
  }
  chi2 += pow(MU(test)/MU_PRIOR, 2) + pow(KAPPA(test)/KAPPA_PRIOR, 2)
    + pow(MU_A(test)/MU_PRIOR, 2) + pow(KAPPA_A(test)/KAPPA_PRIOR, 2)
    + pow(SCATTER_A(test)/MU_PRIOR, 2)
    + pow((SCATTER(test)-SCATTER_CENTER)/SCATTER_PRIOR, 2)
    + pow(((-0.5*SIGMA_A(test)+SIGMA_Z(test))-SIGMA_Z_CENTER)/SIGMA_Z_PRIOR, 2)
    + pow((BH_SCATTER_0(test)-BH_SCATTER_CENTER)/BH_SCATTER_PRIOR, 2);
 
  if (!vel_dispersion) 
  {
  chi2 += pow((BH_BETA_0(test)-BH_BETA_CENTER)/BH_BETA_PRIOR, 2)
    + pow((BH_GAMMA_0(test)-BH_GAMMA_CENTER)/BH_GAMMA_PRIOR, 2);
  } 
  else 
  {
    chi2 += pow((BH_BETA_0(test)-BH_BETA_CENTER_VD)/BH_BETA_PRIOR, 2)
         + pow((BH_GAMMA_0(test)-BH_GAMMA_CENTER_VD)/BH_GAMMA_PRIOR, 2);
  }


   double chi2_smhm = 0;
   for (i = 1; i < num_outputs; i++)
   {
        for (int j=1; j < M_BINS; j++)
        {
                if (steps[i].t[j] >= 2e-8 && steps[i].log_sm[j] > 0 && steps[i].log_sm[j] < steps[i].log_sm[j - 1])
                {
                        double delta = steps[i].log_sm[j - 1] - steps[i].log_sm[j];
                        chi2_smhm += delta * delta;

                }
        }
   }
   chi2 += chi2_smhm;


   if (!no_z_scaling && !no_sfr_constraint) {
    double chi2_sfr = 0;
    double chi2_icl = 0;
    double chi2_kin = 0;
    double chi2_rad = 0;
    double recent_sfr = recent_sfh_in_massive_halos_nocorr();
 
    if (recent_sfr > SFR_CONSTRAINT)
      chi2_sfr = pow((recent_sfr-SFR_CONSTRAINT)/(SFR_CONSTRAINT_WIDTH), 2);

    double recent_icl_star_ratio = recent_Micl_Mstar_ratio_in_massive_halos();
    double recent_radiative_power = recent_radiative_power_in_massive_halos();

    if (recent_icl_star_ratio > ICL_RATIO_CONSTRAINT)
      chi2_icl = pow((recent_icl_star_ratio - ICL_RATIO_CONSTRAINT) / (ICL_RATIO_CONSTRAINT_WIDTH), 2);

    if (recent_radiative_power > RAD_POWER_CONSTRAINT_HIGH)
      chi2_rad = pow((recent_radiative_power - RAD_POWER_CONSTRAINT_HIGH) / (RAD_POWER_CONSTRAINT_WIDTH), 2);


    chi2 += chi2_sfr + chi2_icl + chi2_rad;

  }

  // Note that since we re-allocate space for these spline objects every time due to
  // their varying sizes, we should free them every time after using them to calculate
  // the observables, to avoid out-of-memory kills.
  for (i = 0; i < num_outputs; i++)
  {
    if (steps[i].flag_alloc)
    {
      gsl_spline_free(steps[i].spline); //free the spline object since we are allocating it every time.
      gsl_spline_free(steps[i].spline_sfr);
      if (steps[i].alloc2smf)
      {
        gsl_spline_free(steps[i].spline2); //same for the second segment of the SMHM relation.
        gsl_spline_free(steps[i].spline_sfr2);
      }
    }

    if (steps[i].flag_alloc == 2)
    {
      gsl_spline_free(steps[i].spline_uv); //free the spline object since we are allocating it every time.
      if (steps[i].alloc2uvlf) gsl_spline_free(steps[i].spline_uv2); //same for the second segment of the UV-HM relation.
    }
  }
  return (chi2);
}

// MCMC execution function.
// length is the number of steps for MCMC,
// and run is an identifier to tell if 
// the MCMC is still in the burn-in phase,
// which determines where the MCMC steps are
// output to.
// run == 0 for burn-in phase, otherwise for
// final MCMC.
int mcmc_smf_fit(int length, int run) 
{
  int i, j, dups=0;
  float last_chi2, chi2, first_chi2;
  struct smf_fit last_smf_fit=initial_smf_fit, cur_smf_fit=initial_smf_fit;
  // Determine which file stream to save data, based on the running mode.
  FILE *out_fh = (run > 1) ? stdout : stderr;
  first_chi2 = last_chi2 = all_smf_chi2_err(initial_smf_fit);
  double chi2_repeat = last_chi2;

  // the counter and flag to monitor the repetition of model parameters in MCMC
  int64_t count_repeat = 0;
  int64_t repeated = 0;

  for (i=1; i<length; i++) 
  {
    repeated = 0;
    // Restart the MCMC process if we're still in the burn-in phase
    // and there are already more than 100 repeated steps in the chain.
    if (((count_repeat >= 100) && !run)) 
    {
      fprintf(stderr, "before resetting, i=%d\n", i);
      for (j=0; j<NUM_PARAMS; j++) eigenvalues[j] = eigenvalues_original[j];
      clear_stats();
      count_repeat =0;
      if (i == 100) 
      {
        num_points = dups = i=0;
      }
      set_identity_portion();
      continue;
    }
    random_step(&cur_smf_fit); // Propose new parameter set
    CHI2(cur_smf_fit) = chi2 = all_smf_chi2_err(cur_smf_fit); // Calculate the new chi2
    // Determine if the new step should be accepted.
    if (!isfinite(chi2) || chi2<0 ||
	  (chi2>last_chi2 && drand48()>exp(0.5*inv_temperature*(last_chi2-chi2))))
    {
      chi2 = last_chi2;
      cur_smf_fit = last_smf_fit;
      dups++;
      if (chi2_repeat == last_chi2) repeated = 1;
    }
    chi2_repeat = last_chi2 = chi2;
    last_smf_fit = cur_smf_fit;

    // Add the step (either repeated or not) to the chain and recalculate the
    // statistics, including the mean parameter values and covariance matrix.
    add_to_stats(&cur_smf_fit);
    num_points++;
    // If repeated, add 1 to count_repeat.
    count_repeat = repeated ? count_repeat + 1 : 0;
    // If not in the burn-in phase, we update the input_matrix, eigenvalues, and 
    // orth_matrix after adding the new step.
    if (run > 0) stats_to_step(num_points);
    for (j=0; j<NUM_PARAMS+1; j++) fprintf(out_fh, "%.10f ", cur_smf_fit.params[j]);
    fprintf(out_fh, "%f\n", CHI2(cur_smf_fit));
  }
  initial_smf_fit = cur_smf_fit;
  return dups;
}

// Initialize the MCMC process with input arguments.
void init_mcmc_from_args(int argc, char **argv)
{
  int i;
  r250_init(87L);

  if (argc < 10 || argc > 10000) 
  {
    printf("Usage: %s model matching_scatter obs_scatter systematics nonlinear_luminosity? no_bh_mergers? vel_dispersion? mf_cache file1 file2 ...\n", argv[0]);
    exit(1);
  }


  i = atol(argv[1]);
  if (i == 0) no_z_scaling = 1; // No redshift evolution of any parameters.
  if (i == 2) no_sfr_constraint = 1; // No constraint on the recent SFR histories of massive halos.
  if (i == 3) dplaw_only = 1; // No Gaussian boost in the SFR--Vmax relation.
                              // Note that even in the fiducial model we no longer
                              // have this boost, so this argument is effectively
                              // deprecated.

  if (!atol(argv[2])) { no_matching_scatter = 1; } // No intrinsic scatter in the stellar
                                                    // mass--halo mass relation
  setup_psf(atoi(argv[3])); // Set up the cache for Gaussian distribution.
  if (!atol(argv[3])) no_obs_scatter = 1; // No random scatter in observed vs. true
                                          // stellar mass
  if (!atol(argv[4])) { no_systematics = 1; } // No galaxy systematic effects (like
                                              // mu, kappa)
  nonlinear_luminosity = atol(argv[5]); // Adopt the non-linear scaling relation between
                                        // the radiative and the total Eddington ratio.
  if (atol(argv[6])) { no_bh_mergers = 1; } // No BH mergers.
  if (atol(argv[7])) // Use M-sigma relation instead of the black hole mass--bulge mass
  {                  // relation. This is deprecated for now. See the fprintf() below.
    fprintf(stderr, "M-sigma option is deprecated for now. Automatically switch to black hole mass--bulge mass relation.\n");
    vel_dispersion = 0; 
  }
  
  // Load the cached halo mass functions
  load_mf_cache(argv[8]);

  // Initialize the time steps. See gen_mlists.c
  init_timesteps();

  // Read in the observational data points.
  num_obs_smfs = 0;
  for (i=9; i<argc; i++) 
  {
    if (!strcmp(argv[i], "nopoints")) continue;
    // See observations.c
    load_real_smf(&obs_smfs, &num_obs_smfs, argv[i]);
    if (obs_smfs[num_obs_smfs-1].z_high > z_limit) z_limit = obs_smfs[num_obs_smfs-1].z_high;
  }
}

// Set the scale factor for the identity matrix added to input_matrix 
// each time.
void set_identity_portion(void) 
{
  int64_t i;
  for (i=0; i<NUM_PARAMS; i++) 
  {
    if ((identity_portion > 0.05*fabs(eigenvalues[i])) &&
	  fabs(eigenvalues[i]>0))
      identity_portion = fabs(eigenvalues[i]*0.05);
  }
  fprintf(stderr, "#Set Identity coefficient to %.3e\n", identity_portion);
}

#ifndef GEN_SMF
int main(int argc, char **argv)
{
  int64_t i;
  float temps[5] = {0.1, 0.2, 0.3, 0.5, 0.8};
  int64_t num_temps = 5;
  init_frac_below8();
  gsl_set_error_handler_off(); // Turn off the built-in error handler of GSL.
                               // This is because it kills the program every
                               // time it gets triggered.
  
  init_mcmc_from_args(argc, argv);
  initial_conditions();
  fprintf(stderr, "Total number of data points: %ld", num_obs_smfs);
  // return 0;

  // // Users can turn on the randomness if needed.
  // srand(time(0));
  // srand48(time(0));
  // init_genrand64((unsigned long long) time(0));

  clear_stats(); //Clean all the statistics before starting.
  set_identity_portion();

  // Burn-in.
  for (i=0; i<num_temps; i++) {
    inv_temperature = temps[i];
    fprintf(stderr, "#Temp now: %f\n", 1.0/inv_temperature);
    mcmc_smf_fit(BURN_IN_LENGTH*temps[i], 0);
    stats_to_step(num_points);
    if (i<num_temps-1) clear_stats();
  }
  return 0;
  inv_temperature = 1;

  // Final MCMC
  fprintf(stderr, "#Temp now: %f\n", 1.0/inv_temperature);
  set_identity_portion();
  mcmc_smf_fit(MCMC_LENGTH/4.0, 1);
  set_identity_portion();
  //clear_stats();
  mcmc_smf_fit(MCMC_LENGTH/8.0, 1);
  print_orth_matrix(atol(argv[1]), atol(argv[3]), atof(argv[2]),
		    obs_smfs[num_obs_smfs-1].z_high);
  mcmc_smf_fit(MCMC_LENGTH, 2);
  return 0;
}

#endif

// Initialize the orth_matrix, which specifies how 
// different parameters are correlated with each other.
void init_orth_matrix(void) 
{
  int i,j;
  struct smf_fit address;
  for (i=0; i<NUM_PARAMS; i++)
    for (j=0; j<NUM_PARAMS; j++)
      orth_matrix[i][j] = (i==j) ? 1 : 0;

  if (!no_systematics) 
  {
    i = &(MU(address)) - address.params;
    j = &(EFF_0(address)) - address.params;
    orth_matrix[i][j] = -1; //Build in anti-correlation between MU and EFF
  }
}

// Print the orth_matrix.
void print_orth_matrix(int type, int syst, float scatter, float z) 
{
  FILE *output;
  char buffer[1024];
  int i,j;
  sprintf(buffer, "mcmc_runs/orth_matrix_%d_%d_%f_%f.dat", type, syst, scatter, z);
  output = fopen(buffer, "w");
  for (i=0; i<NUM_PARAMS; i++) 
  {
    fprintf(output, "Eig[%d]: %f\n", i, eigenvalues[i]);
  }
  fprintf(output, "\n");
  for (i=0; i<NUM_PARAMS; i++) 
  {
    for (j=0; j<NUM_PARAMS; j++)
      fprintf(output, "%+.4f ", orth_matrix[i][j]);
    fprintf(output, "\n");
  }
  fclose(output);
}

// Clear the statistics saved in the mean parameter values
// and the covariant matrix.
void clear_stats(void) 
{
  int i, j;
  memset(&smf_sum, 0, sizeof(struct smf_fit));
  for (i=0; i<NUM_PARAMS; i++) 
  {
    smf_sum.params[i] = 0;
    for (j=0; j<NUM_PARAMS; j++) 
    {
      cov_matrix[i][j] = 0;
    }
  }
  num_points = 0;
}

// Add new step to the mean parameter values
// and the covariant matrix.
void add_to_stats(struct smf_fit *a) 
{
  int i,j;
  for (i=0; i<NUM_PARAMS; i++) 
  {
    smf_sum.params[i]+=a->params[i];
    for (j=0; j<NUM_PARAMS; j++) 
    {
      cov_matrix[i][j] += a->params[i]*a->params[j];
    }
  }
}

// Update the input_matrix and eigenvalues
// after adding a new step.
void stats_to_step(int64_t num_stats) 
{
  double mul = 1.0/(double)num_stats;
  double sd_mul = sqrt(2.4*2.4/((double)NUM_PARAMS));
  int i,j;
  double input_matrix[NUM_PARAMS][NUM_PARAMS];

  for (i=0; i<NUM_PARAMS; i++)
    for (j=0; j<NUM_PARAMS; j++)
      input_matrix[i][j] = mul*cov_matrix[i][j] -
	                         mul*mul*smf_sum.params[i]*smf_sum.params[j];

  for (i=0; i<NUM_PARAMS; i++) input_matrix[i][i] += identity_portion;

  jacobi_decompose(input_matrix, eigenvalues, orth_matrix);
  for (i=0; i<NUM_PARAMS; i++) 
  {
    if (eigenvalues[i] <=0 ) eigenvalues[i] = 0;
    eigenvalues[i] = sqrt(eigenvalues[i])*sd_mul;
  }
}

// fix some parameters to certain values as indicated
// by various flags.
void assert_model(struct smf_fit *a) 
{
  LAMBDA(*a) = LAMBDA_A(*a) = 0;
  DELTA_A2(*a) = 0;
  EXPSCALE(*a) = 4;

  if (no_bh_mergers) 
  {
    BH_MERGE_F_0(*a) = -100;
    BH_MERGE_F_1(*a) = 0;
    BH_MERGE_W_0(*a) = 0;
    BH_MERGE_W_1(*a) = 0;
  }

  if (no_z_scaling) 
  {
    M_1_A(*a) = 0;
    M_1_A2(*a) = 0;
    EFF_0_A(*a) = 0;
    EFF_0_A2(*a) = 0;
    EFF_0_A3(*a) = 0;
    EXPSCALE(*a) = 0;
    DELTA_A(*a) = 0;
    DELTA_A2(*a) = 0;
    ALPHA_A(*a) = 0;
    ALPHA_A2(*a) = 0;
    BETA_A(*a) = 0;
    BETA_A2(*a) = 0;
    GAMMA_A(*a) = 0;
    GAMMA_A2(*a) = 0;
    LAMBDA_A(*a) = 0;
    LAMBDA_A2(*a) = 0;
    SCATTER_A(*a) = 0;
    MU_A(*a) = 0;
    KAPPA_A(*a) = 0;
    BURST_DUST(*a) = 0;
    SIGMA_Z(*a) = 0.04;
    SCATTER(*a) = 0.16;
    ALPHA(*a) = -1.490795596100;
    DELTA(*a) = 3.7203;
    GAMMA(*a) = 0.39400;
  }

  if (z_scaling_only) 
  {
    EFF_0(*a) = EFF_0(initial_smf_fit);
    M_1(*a) = M_1(initial_smf_fit);
    DELTA(*a) = DELTA(initial_smf_fit);
    BETA(*a) = BETA(initial_smf_fit);
    GAMMA(*a) = GAMMA(initial_smf_fit);
    LAMBDA(*a) = LAMBDA(initial_smf_fit);
    KAPPA(*a) = KAPPA(initial_smf_fit);
    MU(*a) = MU(initial_smf_fit);
    SCATTER(*a) = SCATTER(initial_smf_fit);
    SIGMA_Z(*a) = SIGMA_Z(initial_smf_fit);    
  }

  if (no_matching_scatter) 
  {
    SCATTER(*a) = 0;
  }
  if (no_obs_scatter) 
  {
    SIGMA_Z(*a) = 0;
  }
  if (no_systematics) 
  {

    KAPPA(*a) = 0;
    MU(*a) = 0;
    KAPPA_A(*a) = 0;
    MU_A(*a) = 0;
    BURST_DUST_AMP(*a) = 0;
    BURST_DUST_Z(*a) = 1;
    BURST_DUST(*a) = 0;
  }

  BH_SCATTER_1(*a) = 0;
  BH_ETA_CRIT_0(*a) = -1.5;
  ABHMF_SHIFT(*a) = 0.2;

  DELTA(*a) = 0.055;
  KAPPA_A(*a) = 0.0;

  BURST_DUST_AMP(*a) = 0;
  BURST_DUST_Z(*a) = 1;
  BURST_DUST(*a) = 0;
  SCATTER_A(*a) = 0;
  SIGMA_A(*a) = 0;
  
  ICL_FRAC_Z(*a) = 0;
  ICL_FRAC_E(*a) = 0;
}

// The function to propose new steps in MCMC.
void random_step(struct smf_fit *a) 
{
  int i,j,num_params=6; // To enhance the sampling efficiency, we randomly
                        // choose 6 parameters to vary.
  for (i=0; i<num_params; i++) 
  {
    j = rand()%NUM_PARAMS;
    // Generate a Gaussian random vector in the eigen vector space, and
    // transform it back into the parameter space.
    vector_madd(a->params, normal_random(0,eigenvalues[j]), orth_matrix[j]);
  }
  assert_model(a);
}
