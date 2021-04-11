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

void init_frac_below8(void)
// Calculate the frac_below8, i.e., the fractional probability that BH in each mass bin
// get scattered below 10^8 Msun.
{
  for (int i=0; i<MBH_BINS; i++)
  {
    double mbh = MBH_MIN + (i + 0.5) * MBH_INV_BPDEX;
    frac_below8[i] = 0.5 * (1 + erf((8.0 - mbh) / bh_scatter_obs * M_SQRT1_2));
    frac_below11[i] = 0.5 * (1 + erf((11.0 - mbh) / bh_scatter_obs * M_SQRT1_2));
    //fprintf(stderr, "i_Mbh=%d, mbh=%f, frac_below8=%e\n", i, mbh, frac_below8[i]);
  }
}



float z_limit = 0;
float z_max = 0;
float z_min = 100;
float inv_temperature = 1;
struct smf_fit initial_smf_fit, smf_sum;
double identity_portion = 3e-5;
double cov_matrix[NUM_PARAMS][NUM_PARAMS];
double orth_matrix[NUM_PARAMS][NUM_PARAMS];
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
           // BH_NORM_STEP, 
           DC_MBH_0_STEP, DC_MBH_W_0_STEP,
          ETA_MU_0_STEP, ETA_MU_1_STEP, BH_DELTA_1_STEP, 
          BH_BETA_0_STEP, BH_BETA_1_STEP, BH_BETA_2_STEP,
           BH_GAMMA_0_STEP, BH_GAMMA_1_STEP, BH_GAMMA_2_STEP,
           1e-9 };

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
           // BH_NORM_STEP, 
           DC_MBH_0_STEP, DC_MBH_W_0_STEP,
          ETA_MU_0_STEP, ETA_MU_1_STEP, BH_DELTA_1_STEP, 
          BH_BETA_0_STEP, BH_BETA_1_STEP, BH_BETA_2_STEP,
           BH_GAMMA_0_STEP, BH_GAMMA_1_STEP, BH_GAMMA_2_STEP,
           1e-9 };

int ind_var[4] = {17, 29, 33, 37};

//void read_params(char *buffer, double *data, int max_n) {
static void read_params(char *buffer, double *data, int max_n) {
  int num_entries = 0;
  char *cur_pos = buffer, *end_pos;
  float val = strtod(cur_pos, &end_pos);
  while (cur_pos != end_pos && num_entries < max_n) {
    data[num_entries] = val;
    num_entries++;
    cur_pos=end_pos;
    while (*cur_pos==' ' || *cur_pos=='\t' || *cur_pos=='\n') cur_pos++;
    val = strtod(cur_pos, &end_pos);
  }
}


void initial_conditions(void) {
  char buffer[2048];
  fgets(buffer, 2048, stdin);
  read_params(buffer, initial_smf_fit.params, NUM_PARAMS);
  if (no_matching_scatter) SCATTER(initial_smf_fit) = 0;
  assert_model(&initial_smf_fit);

  //fprintf(stderr, "Initial Chi2: %e\n", calc_chi2(initial_smf_fit.params));
  fprintf(stderr, "Initial params:\n");
  int i;
  for (i=0; i<NUM_PARAMS-1; i++) fprintf(stderr, "%.12f ", initial_smf_fit.params[i]);
  printf("%.12f\n", initial_smf_fit.params[NUM_PARAMS-1]);

  CHI2(initial_smf_fit) = all_smf_chi2_err(initial_smf_fit);
  double chi22;
  //for (i=0; i<20; i++) 
 //{  
 	//CHI2(initial_smf_fit) = all_smf_chi2_err(initial_smf_fit);
  	chi2_type(initial_smf_fit);
  //}
  init_orth_matrix();
  fprintf(stderr, "Initial Chi2 fit: %f\n", CHI2(initial_smf_fit));
}

float all_smf_chi2_err(struct smf_fit test) {
  float chi2 = 0;
  double chi2_bhmf = 0;
  int i;
  z_limit = 8.5;
  struct smf smf_limit = smhm_at_z(z_limit, test);
  struct smf smf_mid = smhm_at_z(2, test);
  struct smf smf_z4 = smhm_at_z(4.0, test);
  INVALID(test) = 0;
  if (no_z_scaling) {
    if (M_1(test) < 9.5) return -1;
    if (ALPHA(test) < -4) return -1;
    if (GAMMA(test) > 1.05) return -1;
  }
  //  fprintf(stderr, "%f %f\n", smf_limit.alpha, smf_limit.beta);
  //fprintf(stderr, "%f %f %f\n", ALPHA(test), ALPHA_A(test), ALPHA_A2(test));

  if (fabs(KAPPA(test))>KAPPA_PRIOR*2.5 || 
      ALPHA(test) >= -1 || smf_limit.alpha >= -1 || smf_mid.alpha >= -1 ||
      ALPHA(test) <= -15 || smf_limit.alpha <= -15 || smf_mid.alpha <= -15 ||
      BETA(test) < ALPHA(test) || smf_limit.beta < smf_limit.alpha || smf_mid.beta < smf_mid.alpha ||
      BETA(test) > 15 || smf_limit.beta > 15 || smf_mid.beta > 15 ||
      DELTA(test) < 0.01 || smf_limit.delta < 0.01 || smf_mid.delta < 0.01 ||
      fabs(DELTA(test)) > 10 || fabs(smf_limit.delta) > 10 || fabs(smf_mid.delta) > 10 ||
      //fabs(GAMMA(test)) > 50 || fabs(log10(smf_limit.gamma)) > 50 || fabs(log10(smf_mid.gamma)) > 50 ||

      GAMMA(test) < -5 || 
      //LAMBDA(test) > 10 || smf_limit.lambda > 10 || 
      //LAMBDA(test) < -10 || smf_limit.lambda < -10 || 
      SIGMA_Z(test) < 0 || EFF_0(test) > 1 ||
      SIGMA_Z(test) < 0 || EFF_0(test) > 3 ||
      EFF_0_A(test) < -2 ||
      smf_limit.epsilon > 1e5 || smf_mid.epsilon > 1e4 ||
      SIGMA_A(test) > 0 ||
      //      smf_limit.sm_0 > smf_limit.m_1+1 ||
      // fabs(EXPSCALE(test))>20 || EXPSCALE(test) < 0 ||
      //Z_CEIL(test) < 10 ||
      SCATTER(test) > 0.3 ||
      SCATTER(test) < 0 || smf_limit.obs_scatter > 1.0 ||
      // ICL_FRAC(test) > 15 || smf_limit.icl_m > 15 ||
      //ICL_FRAC(test) > 1.5 || ICL_FRAC(test) < 0 ||
      BURST_DUST_Z(test) < 0.8 || BURST_DUST_Z(test)>6 ||
      BURST_DUST(test) < 0 || BURST_DUST(test) > 1 ||
      BURST_DUST_AMP(test) < 0 ||
      BURST_DUST_AMP(test) > 5 || fabs(MU_A(test)) > 0.8 ||
      fabs(KAPPA_A(test)) > KAPPA_LIMIT || SCATTER_A(test) > SCATTER(test) ||
      //ICL_FRAC_E(test) < 0 || ICL_FRAC_E(test)>3 ||
      RHO_05(test) < 0.23 || RHO_05(test) > 1.0 ||

      smf_limit.bh_beta < 4 || smf_limit.bh_beta > 12 ||
      smf_limit.bh_gamma < 0.5 || smf_limit.bh_gamma > 8.5 ||
      BH_ALPHA_0(test) < -2.0 || BH_ALPHA_0(test) >= 1.0 || 
      BH_DELTA_0(test) < -2.0 || BH_DELTA_0(test) > 3.0 ||
      BH_ALPHA_0(test) > BH_DELTA_0(test) ||
      smf_limit.bh_alpha > smf_limit.bh_delta || 
      fabs(smf_limit.bh_alpha) > 4 ||
      fabs(smf_limit.bh_alpha) >= 1.0 ||
      fabs(smf_limit.bh_delta) > 4 ||

      smf_limit.icl_frac < 1e-20 || smf_limit.icl_frac > 1 ||
      ICL_FRAC(test) < -20 || ICL_FRAC(test) > 0 ||
      QM_0(test) < 0 || QM_0(test) > 8 ||
      smf_limit.qm < 0 ||
      QWIDTH_0(test) < 0.001 || QWIDTH_0(test) > 2 ||
      smf_limit.qwidth < 0.001 || smf_limit.qwidth > 2 ||

      BH_MERGE_F_0(test) < -5 || BH_MERGE_F_0(test) > 5 ||
      smf_limit.f_merge_bh < 1e-5 || smf_limit.f_merge_bh > 1e5 || 
      //smf_limit.bh_merge < 12 || smf_limit.bh_merge > 17 ||
      //BH_MERGE_W_0(test) < 0.001 || BH_MERGE_W_0(test) > 4 ||
      //smf_limit.bh_merge_width < 0.001 || smf_limit.bh_merge_width > 4 ||
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
   //fprintf(stderr, "Is ok!\n");
#pragma omp parallel
  {
    calc_sfh(&test);

    double bhz0 = calc_bhmf(BH_MASS_TO_REQUIRE_ND_GROWTH, 0);
    double bhz1 = calc_bhmf(BH_MASS_TO_REQUIRE_ND_GROWTH, 1);
    // if (bhz0 < bhz1)
    //  fprintf(stderr, "bhz0: %.6e, bhz1: %.6e\n", bhz0, bhz1); 
    if ((bhz1 < 0) || (bhz0 < 0) || (bhz0 < bhz1)) {
      //chi2_bhmf = 200 * pow(bhz1 / bhz0, 2);
      //chi2 += chi2_bhmf;
      
      //fprintf(stderr, "BH ND down. bhz0: %.6e, bhz1: %.6e\n", bhz0, bhz1);
      INVALIDATE(&test, "BH Number Density Decreased from z=1 to z=0!");
    }
    //fprintf(stderr, "chi2_bhmf: %f\n", chi2_bhmf);
    //fprintf(stderr, "After calc_sfh(), invalid=%f\n", INVALID(test));
    //for (int j=0; j<NUM_PARAMS+1; j++) fprintf(stderr, "%.10f ", test.params[j]);
    //fprintf(stderr, "%f\n", CHI2(test));
#pragma omp for schedule(simd:dynamic, 8) //private(m,epsilon,smf_val)
    for (i=0; i<num_obs_smfs; i++)
      if (!(INVALID(test)))
      {
        obs_smfs[i].chi2 = calc_single_chi2_err_point(&(obs_smfs[i]));
        
	 // fprintf(stderr, "z_low=%f, z_high=%f, mass=%f, type=%d, chi2=%f\n", obs_smfs[i].z_low, obs_smfs[i].z_high, obs_smfs[i].mass, obs_smfs[i].type, obs_smfs[i].chi2);
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
          // gsl_spline_free(steps[i].spline_uv); //free the spline object since we are allocating it every time.
          // if (steps[i].alloc2uvlf) gsl_spline_free(steps[i].spline_uv2); //same for the second segment of the UV-HM relation.
          // // gsl_spline_free(stseps[i].spline_std_uv);

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
    //if (obs_smfs[i].type == UVLF_TYPE)
    //	fprintf(stderr, "z_low=%f, z_high=%f, val=%e, chi2=%f\n", obs_smfs[i].z_low, obs_smfs[i].z_high, obs_smfs[i].val, obs_smfs[i].chi2);
  }

  //fprintf(stderr, "Before adding priors, chi2=%f\n", chi2);


  chi2 += pow(MU(test)/MU_PRIOR, 2) + pow(KAPPA(test)/KAPPA_PRIOR, 2)
    + pow(MU_A(test)/MU_PRIOR, 2) + pow(KAPPA_A(test)/KAPPA_PRIOR, 2)
    + pow(SCATTER_A(test)/MU_PRIOR, 2)
    + pow((SCATTER(test)-SCATTER_CENTER)/SCATTER_PRIOR, 2)
    + pow(((-0.5*SIGMA_A(test)+SIGMA_Z(test))-SIGMA_Z_CENTER)/SIGMA_Z_PRIOR, 2);
    // + pow((BH_SCATTER_0(test)-BH_SCATTER_CENTER)/BH_SCATTER_PRIOR, 2);
  //if (Z_CEIL(test) < Z_CEIL_CENTER) 
  //{
//	  chi2 += pow((Z_CEIL(test) - Z_CEIL_CENTER) / Z_CEIL_WIDTH, 2);
    // + pow((BH_SCATTER_0(test)-BH_SCATTER_CENTER)/BH_SCATTER_PRIOR, 2);
  //fprintf(stderr, "chi2_z_ceil: %f\n", pow((Z_CEIL(test) - Z_CEIL_CENTER) / Z_CEIL_WIDTH, 2));
  //} 
   if (!vel_dispersion) {
   chi2 += pow((BH_BETA_0(test)-BH_BETA_CENTER)/BH_BETA_PRIOR, 2)
     + pow((BH_GAMMA_0(test)-BH_GAMMA_CENTER)/BH_GAMMA_PRIOR, 2);
 } else {
   chi2 += pow((BH_BETA_0(test)-BH_BETA_CENTER_VD)/BH_BETA_PRIOR, 2)
     + pow((BH_GAMMA_0(test)-BH_GAMMA_CENTER_VD)/BH_GAMMA_PRIOR, 2);
 }

  //fprintf(stderr, "after adding priors on parameters, chi2=%f\n", chi2);

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
   //fprintf(stderr, "chi2_bhar=%f, chi2_smhm=%f\n", chi2_bhar, chi2_smhm);

   if (!no_z_scaling && !no_sfr_constraint) {
    double chi2_sfr = 0;
    double chi2_icl = 0;
    double chi2_kin = 0;
    double chi2_rad = 0;
    double chi2_qso = 0;
    
    double penalty_rising_sfh = rising_sfh_penalty();
    fprintf(stderr, "The chi2 penalty for rising SFH in massive halos at low-z: %f\n", penalty_rising_sfh); 

    double ratio_z6_qso = ratio_high_z_qso(5.7, 6.5, 47, area_sdss);
    fprintf(stderr, "The ratio of quasars with Lbol>10^47, Mbh<10^8 vs. Mbh>10^8, and 5.7 < z < 6.5 is %e\n", ratio_z6_qso);


    double n_z6_qso = number_high_z_low_mass_qso(5.7, 6.5, 47, area_sdss);
    fprintf(stderr, "The number of quasars with Lbol>10^47, Mbh<10^8, and 5.7 < z < 6.5 is %e\n", n_z6_qso);

    if (n_z6_qso > 0) chi2_qso = 2 * n_z6_qso;
    chi2_qso += log10(ratio_z6_qso) * 10;

    double recent_sfr = recent_sfh_in_massive_halos();
    //double recent_sfr = recent_sfh_in_massive_halos_nocorr();
    //fprintf(stderr, "recent sfr: %.6e\n", recent_sfr);
    if (recent_sfr > SFR_CONSTRAINT)
      chi2_sfr = pow((recent_sfr-SFR_CONSTRAINT)/(SFR_CONSTRAINT_WIDTH), 2);
      //chi2 += pow((recent_sfr-SFR_CONSTRAINT)/(SFR_CONSTRAINT_WIDTH), 2);
    //fprintf(stderr, "recent sfr: %.6e\n", recent_sfr);
    double recent_icl_star_ratio = recent_Micl_Mstar_ratio_in_massive_halos();
    double recent_radiative_power = recent_radiative_power_in_massive_halos();
    // double recent_kinetic_frac = recent_kinetic_frac_in_massive_halos();
    //fprintf(stderr, "recent_icl_star_ratio: %.3f\n", recent_icl_star_ratio);
    if (recent_icl_star_ratio > ICL_RATIO_CONSTRAINT)
      chi2_icl = pow((recent_icl_star_ratio - ICL_RATIO_CONSTRAINT) / (ICL_RATIO_CONSTRAINT_WIDTH), 2);
      //chi2 += pow((recent_icl_star_ratio - ICL_RATIO_CONSTRAINT) / (ICL_RATIO_CONSTRAINT_WIDTH), 2);
    //double recent_sfr_nocorr = recent_sfh_in_massive_halos_nocorr();
    //if (recent_sfr_nocorr > SFR_CONSTRAINT)
      //chi2 += pow((recent_sfr_nocorr-SFR_CONSTRAINT)/(SFR_CONSTRAINT_WIDTH), 2);

    //if (recent_kinetic_frac < KIN_FRAC_CONSTRAINT_LOW)
      //chi2_kin = pow((recent_kinetic_frac - KIN_FRAC_CONSTRAINT_LOW) / (KIN_POWER_CONSTRAINT_WIDTH), 2);

    // if (recent_kinetic_power < KIN_POWER_CONSTRAINT_LOW)
    //   chi2_kin = pow((recent_kinetic_power - KIN_POWER_CONSTRAINT_LOW) / (KIN_POWER_CONSTRAINT_WIDTH), 2);
    // else if (recent_kinetic_power > KIN_POWER_CONSTRAINT_HIGH)
    //   chi2_kin = pow((recent_kinetic_power - KIN_POWER_CONSTRAINT_HIGH) / (KIN_POWER_CONSTRAINT_WIDTH), 2);


     // if (recent_radiative_power < RAD_POWER_CONSTRAINT_LOW)
     //   chi2_rad = pow((recent_radiative_power - RAD_POWER_CONSTRAINT_LOW) / (RAD_POWER_CONSTRAINT_WIDTH), 2);
    // else 
    if (recent_radiative_power > RAD_POWER_CONSTRAINT_HIGH)
      chi2_rad = pow((recent_radiative_power - RAD_POWER_CONSTRAINT_HIGH) / (RAD_POWER_CONSTRAINT_WIDTH), 2);


    chi2 += chi2_sfr + chi2_icl + chi2_qso + chi2_rad;
    chi2 += penalty_rising_sfh;
    fprintf(stderr, "recent rad power: %f, chi2_rad:%f\n", recent_radiative_power, chi2_rad);
    //fprintf(stderr, "chi2_sfr: %.6f, chi2_icl=%f, chi2_rad: %.6f\n", chi2_sfr, chi2_icl, chi2_rad);
  }
  //fprintf(stderr, "After adding all priors, chi2=%f\n", chi2);

  //double chi2_merge = 0;
  //for (i = 0; i < num_outputs; i++)
  //{
//	for (int j = 0; j < M_BINS; j++)
//	{
//		if (steps[i].new_bh_mass[j] > 0 && steps[i].new_bh_mass[j] < steps[i].bh_merged[j]) chi2_merge += fabs(steps[i].bh_merged[j] / steps[i].new_bh_mass[j]);
//	}
  //}
  //fprintf(stderr, "chi2_merge=%f\n", chi2_merge); 
  //chi2 += chi2_merge;
  //if (isfinite(chi2) && chi2 > 0) chi2 = chi2_merge;

  //double chi2_csfr = 0;
  //chi2_csfr = pow(log10(steps[24].observed_cosmic_sfr / steps[4].observed_cosmic_sfr) / 0.5, 2);
  //chi2 += chi2_csfr;
  //fprintf(stderr, "chi2_csfr=%f\n", chi2_csfr);
  
  //double chi2_sm = 0;
  //chi2_sm += steps[24].log_sm[19] > 9? 0 : 100 * exp10(9 - steps[24].log_sm[19]);
  //chi2_sm += steps[24].log_sm[24] > 10? 0 : 100 * exp10(10 - steps[24].log_sm[24]);
  //chi2 += chi2_sm;
  //fprintf(stderr, "chi2_sm=%f\n", chi2_sm);

  //fprintf(stderr, "chi2_bhmf=%f\n", chi2_bhmf);

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
          // gsl_spline_free(steps[i].spline_uv); //free the spline object since we are allocating it every time.
          // if (steps[i].alloc2uvlf) gsl_spline_free(steps[i].spline_uv2); //same for the second segment of the UV-HM relation.
          // // gsl_spline_free(stseps[i].spline_std_uv);

        }

        if (steps[i].flag_alloc == 2)
        {
          gsl_spline_free(steps[i].spline_uv); //free the spline object since we are allocating it every time.
          if (steps[i].alloc2uvlf) gsl_spline_free(steps[i].spline_uv2); //same for the second segment of the UV-HM relation.
          // gsl_spline_free(stseps[i].spline_std_uv);
        }
    //if ((1 / steps[i].scale - 1 > 6) && (steps[i].cosmic_sfr <= steps[i-1].cosmic_sfr))
     //{
       //double ratio = steps[i-1].cosmic_sfr / steps[i].cosmic_sfr;
       //chi2 += ratio * ratio;
       //chi2_csfr += ratio * ratio;
     //}
  }
   //fprintf(stderr, "chi2_csfr: %.6f\n", chi2_csfr);

  // printf("The chi2 calculated by all_smf_chi2_err is: %f\n", chi2);
  //if (chi2_bhmf == 0 && isfinite(chi2))
//{
//	for (i = 0; i < NUM_PARAMS; i++) fprintf(stderr, "%.12f ", test.params[i]);
//	fprintf(stderr, "\n");
//}

  return (chi2);
}

float all_smf_chi2_err_write(struct smf_fit test) {
  float chi2 = 0;
  double chi2_bhmf = 0;
  //double chi2_type[10] = {0};
  int i;
  z_limit = 8.5;
  struct smf smf_limit = smhm_at_z(z_limit, test);
  struct smf smf_mid = smhm_at_z(2, test);
  struct smf smf_z4 = smhm_at_z(4.0, test);
  INVALID(test) = 0;
  if (no_z_scaling) {
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

      GAMMA(test) < -5 || 

      SIGMA_Z(test) < 0 || EFF_0(test) > 1 ||
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

      smf_limit.bh_beta < 4 || smf_limit.bh_beta > 12 ||
      smf_limit.bh_gamma < 0.5 || smf_limit.bh_gamma > 8.5 ||
      BH_ALPHA_0(test) < -0.5 || BH_ALPHA_0(test) > 3.0 || 
      BH_DELTA_0(test) < -0.5 || BH_DELTA_0(test) > 3.0 ||
      BH_ALPHA_0(test) > BH_DELTA_0(test) ||
      smf_limit.bh_alpha > smf_limit.bh_delta || 
      fabs(smf_limit.bh_alpha) > 4 ||
      fabs(smf_limit.bh_delta) > 4 ||

      smf_limit.icl_frac < 1e-20 || smf_limit.icl_frac > 1 ||
      ICL_FRAC(test) < -20 || ICL_FRAC(test) > 0 ||
      QM_0(test) < 9 || QM_0(test) > 14 ||
      smf_limit.qm < 9 ||
      QWIDTH_0(test) < 0.001 || QWIDTH_0(test) > 2 ||
      smf_limit.qwidth < 0.001 || smf_limit.qwidth > 2 ||

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
          for (i=0; i<17; i++) fprintf(stderr, "-1 ");
          fprintf(stderr, "\n");
          return (-1);
      }
#pragma omp parallel
  {
    calc_sfh(&test);

    double bhz0 = calc_bhmf(BH_MASS_TO_REQUIRE_ND_GROWTH, 0);
    double bhz1 = calc_bhmf(BH_MASS_TO_REQUIRE_ND_GROWTH, 1);

    if ((bhz1 < 0) || (bhz0 < 0) || (bhz0 < bhz1)) {
      for (i=0; i<17; i++) fprintf(stderr, "-1 ");
      fprintf(stderr, "\n");
      INVALIDATE(&test, "BH Number Density Decreased from z=1 to z=0!");
    }

#pragma omp for schedule(simd:dynamic, 8) //private(m,epsilon,smf_val)
    for (i=0; i<num_obs_smfs; i++)
      if (!(INVALID(test)))
      {
        obs_smfs[i].chi2 = calc_single_chi2_err_point(&(obs_smfs[i]));
        
   // fprintf(stderr, "z_low=%f, z_high=%f, mass=%f, type=%d, chi2=%f\n", obs_smfs[i].z_low, obs_smfs[i].z_high, obs_smfs[i].mass, obs_smfs[i].type, obs_smfs[i].chi2);
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
          // gsl_spline_free(stseps[i].spline_std_uv);
        }
      }
      for (i=0; i<17; i++) fprintf(stderr, "-1 ");
      fprintf(stderr, "\n");
      return -1;
    }
  double chi2_type[11] = {0};
  for (i=0; i<num_obs_smfs; i++) 
  {
    chi2 += obs_smfs[i].chi2;
    chi2_type[obs_smfs[i].type] += obs_smfs[i].chi2;
    //if (obs_smfs[i].type == UVLF_TYPE)
    //fprintf(stderr, "z_low=%f, z_high=%f, val=%e, chi2=%f\n", obs_smfs[i].z_low, obs_smfs[i].z_high, obs_smfs[i].val, obs_smfs[i].chi2);
  }

  //fprintf(stderr, "Before adding priors, chi2=%f\n", chi2);


  chi2 += pow(MU(test)/MU_PRIOR, 2) + pow(KAPPA(test)/KAPPA_PRIOR, 2)
    + pow(MU_A(test)/MU_PRIOR, 2) + pow(KAPPA_A(test)/KAPPA_PRIOR, 2)
    + pow(SCATTER_A(test)/MU_PRIOR, 2)
    + pow((SCATTER(test)-SCATTER_CENTER)/SCATTER_PRIOR, 2)
    + pow(((-0.5*SIGMA_A(test)+SIGMA_Z(test))-SIGMA_Z_CENTER)/SIGMA_Z_PRIOR, 2);

   if (!vel_dispersion) {
   chi2 += pow((BH_BETA_0(test)-BH_BETA_CENTER)/BH_BETA_PRIOR, 2)
     + pow((BH_GAMMA_0(test)-BH_GAMMA_CENTER)/BH_GAMMA_PRIOR, 2);
 } else {
   chi2 += pow((BH_BETA_0(test)-BH_BETA_CENTER_VD)/BH_BETA_PRIOR, 2)
     + pow((BH_GAMMA_0(test)-BH_GAMMA_CENTER_VD)/BH_GAMMA_PRIOR, 2);
 }

  //fprintf(stderr, "after adding priors on parameters, chi2=%f\n", chi2);

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
   //fprintf(stderr, "chi2_bhar=%f, chi2_smhm=%f\n", chi2_bhar, chi2_smhm);
   double chi2_sfr, chi2_icl, chi2_rad;
   if (!no_z_scaling && !no_sfr_constraint) {
    //double chi2_sfr = 0;
    //double chi2_icl = 0;
    double chi2_kin = 0;
    //double chi2_rad = 0;
    double recent_sfr = recent_sfh_in_massive_halos();
    //double recent_sfr = recent_sfh_in_massive_halos_nocorr();
    //fprintf(stderr, "recent sfr: %.6e\n", recent_sfr);
    if (recent_sfr > SFR_CONSTRAINT)
      chi2_sfr = pow((recent_sfr-SFR_CONSTRAINT)/(SFR_CONSTRAINT_WIDTH), 2);
      //chi2 += pow((recent_sfr-SFR_CONSTRAINT)/(SFR_CONSTRAINT_WIDTH), 2);
    //fprintf(stderr, "recent sfr: %.6e\n", recent_sfr);
    double recent_icl_star_ratio = recent_Micl_Mstar_ratio_in_massive_halos();
    double recent_radiative_power = recent_radiative_power_in_massive_halos();
    // double recent_kinetic_frac = recent_kinetic_frac_in_massive_halos();
    //fprintf(stderr, "recent_icl_star_ratio: %.3f\n", recent_icl_star_ratio);
    if (recent_icl_star_ratio > ICL_RATIO_CONSTRAINT)
      chi2_icl = pow((recent_icl_star_ratio - ICL_RATIO_CONSTRAINT) / (ICL_RATIO_CONSTRAINT_WIDTH), 2);

    if (recent_radiative_power > RAD_POWER_CONSTRAINT_HIGH)
      chi2_rad = pow((recent_radiative_power - RAD_POWER_CONSTRAINT_HIGH) / (RAD_POWER_CONSTRAINT_WIDTH), 2);


    chi2 += chi2_sfr + chi2_icl + chi2_rad;

  }
  double bhz0 = calc_bhmf(BH_MASS_TO_REQUIRE_ND_GROWTH, 0);
    double bhz1 = calc_bhmf(BH_MASS_TO_REQUIRE_ND_GROWTH, 1);
  for (i=0; i < 11; i++) fprintf(stderr, "%.6f ", chi2_type[i]);
  fprintf(stderr, "%.6f %.6f %.6f %.6f %.6f %.6e %.6e\n", chi2_smhm, chi2_bhar, chi2_sfr, chi2_icl, chi2_rad, bhz0, bhz1);

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
          // gsl_spline_free(steps[i].spline_uv); //free the spline object since we are allocating it every time.
          // if (steps[i].alloc2uvlf) gsl_spline_free(steps[i].spline_uv2); //same for the second segment of the UV-HM relation.
          // // gsl_spline_free(stseps[i].spline_std_uv);

        }

        if (steps[i].flag_alloc == 2)
        {
          gsl_spline_free(steps[i].spline_uv); //free the spline object since we are allocating it every time.
          if (steps[i].alloc2uvlf) gsl_spline_free(steps[i].spline_uv2); //same for the second segment of the UV-HM relation.
          // gsl_spline_free(stseps[i].spline_std_uv);
        }
  }

  return (chi2);
}

float chi2_type(struct smf_fit test) {
  float chi2 = 0;
  double chi2_bhmf = 0;
  int i;
  z_limit = 8.5;
  struct smf smf_limit = smhm_at_z(z_limit, test);
  struct smf smf_mid = smhm_at_z(2, test);
  struct smf smf_z4 = smhm_at_z(4.0, test);
  INVALID(test) = 0;
  if (no_z_scaling) {
    if (M_1(test) < 9.5) return -1;
    if (ALPHA(test) < -4) return -1;
    if (GAMMA(test) > 1.05) return -1;
  }
  //  fprintf(stderr, "%f %f\n", smf_limit.alpha, smf_limit.beta);
  //fprintf(stderr, "%f %f %f\n", ALPHA(test), ALPHA_A(test), ALPHA_A2(test));

  if (fabs(KAPPA(test))>KAPPA_PRIOR*2.5 || 
      ALPHA(test) >= -1 || smf_limit.alpha >= -1 || smf_mid.alpha >= -1 ||
      ALPHA(test) <= -15 || smf_limit.alpha <= -15 || smf_mid.alpha <= -15 ||
      BETA(test) < ALPHA(test) || smf_limit.beta < smf_limit.alpha || smf_mid.beta < smf_mid.alpha ||
      BETA(test) > 15 || smf_limit.beta > 15 || smf_mid.beta > 15 ||
      DELTA(test) < 0.01 || smf_limit.delta < 0.01 || smf_mid.delta < 0.01 ||
      fabs(DELTA(test)) > 10 || fabs(smf_limit.delta) > 10 || fabs(smf_mid.delta) > 10 ||
      //fabs(GAMMA(test)) > 50 || fabs(log10(smf_limit.gamma)) > 50 || fabs(log10(smf_mid.gamma)) > 50 ||

      // GAMMA(test) < -1 || 
      LAMBDA(test) > 10 || smf_limit.lambda > 10 || 
      LAMBDA(test) < -10 || smf_limit.lambda < -10 || 
      // SIGMA_Z(test) < 0 || EFF_0(test) > 1 ||
      SIGMA_Z(test) < 0 || EFF_0(test) > 3 ||
      smf_limit.epsilon > 1e5 || smf_mid.epsilon > 1e4 ||
      SIGMA_A(test) > 0 ||
      //      smf_limit.sm_0 > smf_limit.m_1+1 ||
      // fabs(EXPSCALE(test))>20 || EXPSCALE(test) < 0 ||
      //Z_CEIL(test) < 10 ||
      SCATTER(test) > 0.3 ||
      SCATTER(test) < 0 || smf_limit.obs_scatter > 1.0 ||
      // ICL_FRAC(test) > 15 || smf_limit.icl_m > 15 ||
      //ICL_FRAC(test) > 1.5 || ICL_FRAC(test) < 0 ||
      BURST_DUST_Z(test) < 0.8 || BURST_DUST_Z(test)>6 ||
      BURST_DUST(test) < 0 || BURST_DUST(test) > 1 ||
      BURST_DUST_AMP(test) < 0 ||
      BURST_DUST_AMP(test) > 5 || fabs(MU_A(test)) > 0.8 ||
      fabs(KAPPA_A(test)) > KAPPA_LIMIT || SCATTER_A(test) > SCATTER(test) ||
      //ICL_FRAC_E(test) < 0 || ICL_FRAC_E(test)>3 ||
      RHO_05(test) < 0.23 || RHO_05(test) > 1.0 ||

      smf_limit.bh_beta < 4 || smf_limit.bh_beta > 12 ||
      smf_limit.bh_gamma < 0.5 || smf_limit.bh_gamma > 8.5 ||
      BH_ALPHA_0(test) < -0.5 || BH_ALPHA_0(test) > 3.0 || 
      BH_DELTA_0(test) < -0.5 || BH_DELTA_0(test) > 3.0 ||
      BH_ALPHA_0(test) > BH_DELTA_0(test) ||
      smf_limit.bh_alpha > smf_limit.bh_delta || 
      fabs(smf_limit.bh_alpha) > 4 ||


      smf_limit.icl_frac < 1e-20 || smf_limit.icl_frac > 1 ||
      ICL_FRAC(test) < -20 || ICL_FRAC(test) > 0 ||
      QM_0(test) < 9 || QM_0(test) > 14 ||
      smf_limit.qm < 9 ||
      QWIDTH_0(test) < 0.001 || QWIDTH_0(test) > 2 ||
      smf_limit.qwidth < 0.001 || smf_limit.qwidth > 2||


      //BH_MERGE_F_0(test) < 12 || BH_MERGE_F_0(test) > 17 ||
      //smf_limit.bh_merge < 12 || smf_limit.bh_merge > 17 ||
      //BH_MERGE_W_0(test) < 0.001 || BH_MERGE_W_0(test) > 4 ||
      //smf_limit.bh_merge_width < 0.001 || smf_limit.bh_merge_width > 4 ||
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
        //fprintf(stderr, "Is NOT OK!\n");
          //return (-1);
      }
  //fprintf(stderr, "Is ok!\n");
#pragma omp parallel
  {
    calc_sfh(&test);

    // double bhz0 = calc_bhmf(BH_MASS_TO_REQUIRE_ND_GROWTH, 0);
    // double bhz1 = calc_bhmf(BH_MASS_TO_REQUIRE_ND_GROWTH, 1);
    // // if (bhz0 < bhz1)
    //  //fprintf(stderr, "bhz0: %.6e, bhz1: %.6e\n", bhz0, bhz1); 
    // if ((bhz0 < 0) || (bhz0 < bhz1)) {
    //   //chi2_bhmf = pow(bhz1 / bhz0, 2);
    //   //chi2 += chi2_bhmf;
      
    //     fprintf(stderr, "BH ND down. bhz0: %.6e, bhz1: %.6e\n", bhz0, bhz1);
    //   INVALIDATE(&test, "BH Number Density Decreased from z=1 to z=0!");
    // }
    //fprintf(stderr, "chi2_bhmf: %f\n", chi2_bhmf);
    //fprintf(stderr, "After calc_sfh(), invalid=%f\n", INVALID(test));
    //for (int j=0; j<NUM_PARAMS+1; j++) fprintf(stderr, "%.10f ", test.params[j]);
    //fprintf(stderr, "%f\n", CHI2(test));
    
#pragma omp for schedule(simd:dynamic, 8) //private(m,epsilon,smf_val)
    
    

    for (i=0; i<num_obs_smfs; i++)
      // if (!(INVALID(test)))
      {
        obs_smfs[i].chi2 = calc_single_chi2_err_point(&(obs_smfs[i]));
        
   // fprintf(stderr, "z_low=%f, z_high=%f, mass=%f, type=%d, chi2=%f\n", obs_smfs[i].z_low, obs_smfs[i].z_high, obs_smfs[i].mass, obs_smfs[i].type, obs_smfs[i].chi2);
      }
    
  }

  double chi2_type[11] = {0};
  double chi2_qpdf_z[5] = {0};
  double zlows[5] = {0.1, 0.5, 1.0, 1.5, 2.0};
  for (i=0; i<num_obs_smfs; i++) 
    {
      chi2 += obs_smfs[i].chi2;
      chi2_type[obs_smfs[i].type] += obs_smfs[i].chi2;
      //for (int ii=0; ii<5; ii++)
      //  if (obs_smfs[i].type == QPDF_ETA_TYPE && obs_smfs[i].z_low == zlows[ii])
       // {
        //fprintf(stderr, "the %d-th data point is at z_low=%.2f\n", i, obs_smfs[i].z_low);
        //chi2_qpdf_z[ii] += obs_smfs[i].chi2;
        //}
	//printf("z_low=%f, z_high=%f, mass=%f, type=%d, chi2=%f\n", obs_smfs[i].z_low, obs_smfs[i].z_high, obs_smfs[i].mass, obs_smfs[i].type, obs_smfs[i].chi2);

    }

    // fprintf(stderr, "The total chi2 for all the types before the priors are: %f\n", chi2);
  //for (i=0; i<5; i++)
    //    fprintf(stderr, "The total chi2 from QPDF_ETA at zlow=%.2f is: %.6f\n",
//zlows[i], chi2_qpdf_z[i]);
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
          // gsl_spline_free(steps[i].spline_uv); //free the spline object since we are allocating it every time.
          // if (steps[i].alloc2uvlf) gsl_spline_free(steps[i].spline_uv2); //same for the second segment of the UV-HM relation.
          // // gsl_spline_free(stseps[i].spline_std_uv);

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
    //if (obs_smfs[i].type == UVLF_TYPE)
    //  fprintf(stderr, "z_low=%f, z_high=%f, val=%e, chi2=%f\n", obs_smfs[i].z_low, obs_smfs[i].z_high, obs_smfs[i].val, obs_smfs[i].chi2);
  }
  chi2 += pow(MU(test)/MU_PRIOR, 2) + pow(KAPPA(test)/KAPPA_PRIOR, 2)
    + pow(MU_A(test)/MU_PRIOR, 2) + pow(KAPPA_A(test)/KAPPA_PRIOR, 2)
    + pow(SCATTER_A(test)/MU_PRIOR, 2)
    + pow((SCATTER(test)-SCATTER_CENTER)/SCATTER_PRIOR, 2)
    + pow(((-0.5*SIGMA_A(test)+SIGMA_Z(test))-SIGMA_Z_CENTER)/SIGMA_Z_PRIOR, 2)
    + pow((BH_SCATTER_0(test)-BH_SCATTER_CENTER)/BH_SCATTER_PRIOR, 2);
  //if (Z_CEIL(test) < Z_CEIL_CENTER) 
  //{
//    chi2 += pow((Z_CEIL(test) - Z_CEIL_CENTER) / Z_CEIL_WIDTH, 2);
    // + pow((BH_SCATTER_0(test)-BH_SCATTER_CENTER)/BH_SCATTER_PRIOR, 2);
  //fprintf(stderr, "chi2_z_ceil: %f\n", pow((Z_CEIL(test) - Z_CEIL_CENTER) / Z_CEIL_WIDTH, 2));
  //} 
  if (!vel_dispersion) {
  chi2 += pow((BH_BETA_0(test)-BH_BETA_CENTER)/BH_BETA_PRIOR, 2)
    + pow((BH_GAMMA_0(test)-BH_GAMMA_CENTER)/BH_GAMMA_PRIOR, 2);
} else {
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
    //fprintf(stderr, "recent sfr: %.6e\n", recent_sfr);
    if (recent_sfr > SFR_CONSTRAINT)
      chi2_sfr = pow((recent_sfr-SFR_CONSTRAINT)/(SFR_CONSTRAINT_WIDTH), 2);
      //chi2 += pow((recent_sfr-SFR_CONSTRAINT)/(SFR_CONSTRAINT_WIDTH), 2);
    //fprintf(stderr, "recent sfr: %.6e\n", recent_sfr);
    double recent_icl_star_ratio = recent_Micl_Mstar_ratio_in_massive_halos();
    double recent_radiative_power = recent_radiative_power_in_massive_halos();
    //double recent_kinetic_frac = recent_kinetic_frac_in_massive_halos();
    //fprintf(stderr, "recent_icl_star_ratio: %.3f\n", recent_icl_star_ratio);
    if (recent_icl_star_ratio > ICL_RATIO_CONSTRAINT)
      chi2_icl = pow((recent_icl_star_ratio - ICL_RATIO_CONSTRAINT) / (ICL_RATIO_CONSTRAINT_WIDTH), 2);
      //chi2 += pow((recent_icl_star_ratio - ICL_RATIO_CONSTRAINT) / (ICL_RATIO_CONSTRAINT_WIDTH), 2);
    //double recent_sfr_nocorr = recent_sfh_in_massive_halos_nocorr();
    //if (recent_sfr_nocorr > SFR_CONSTRAINT)
      //chi2 += pow((recent_sfr_nocorr-SFR_CONSTRAINT)/(SFR_CONSTRAINT_WIDTH), 2);

    //if (recent_kinetic_frac < KIN_FRAC_CONSTRAINT_LOW)
      //chi2_kin = pow((recent_kinetic_frac - KIN_FRAC_CONSTRAINT_LOW) / (KIN_POWER_CONSTRAINT_WIDTH), 2);

    // if (recent_kinetic_power < KIN_POWER_CONSTRAINT_LOW)
    //   chi2_kin = pow((recent_kinetic_power - KIN_POWER_CONSTRAINT_LOW) / (KIN_POWER_CONSTRAINT_WIDTH), 2);
    // else if (recent_kinetic_power > KIN_POWER_CONSTRAINT_HIGH)
    //   chi2_kin = pow((recent_kinetic_power - KIN_POWER_CONSTRAINT_HIGH) / (KIN_POWER_CONSTRAINT_WIDTH), 2);


    // if (recent_radiative_power < RAD_POWER_CONSTRAINT_LOW)
    //   chi2_rad = pow((recent_radiative_power - RAD_POWER_CONSTRAINT_LOW) / (RAD_POWER_CONSTRAINT_WIDTH), 2);
    // else 
    if (recent_radiative_power > RAD_POWER_CONSTRAINT_HIGH)
      chi2_rad = pow((recent_radiative_power - RAD_POWER_CONSTRAINT_HIGH) / (RAD_POWER_CONSTRAINT_WIDTH), 2);


    chi2 += chi2_sfr + chi2_icl + chi2_rad;
    //fprintf(stderr, "recent rad power: %f, chi2_rad:%f\n", recent_radiative_power, chi2_rad);
    //fprintf(stderr, "chi2_kin: %.6f, chi2_rad: %.6f\n", chi2_kin, chi2_rad);
  }

  //double chi2_merge = 0;
  //for (i = 0; i < num_outputs; i++)
  //{
//  for (int j = 0; j < M_BINS; j++)
//  {
//    if (steps[i].new_bh_mass[j] > 0 && steps[i].new_bh_mass[j] < steps[i].bh_merged[j]) chi2_merge += fabs(steps[i].bh_merged[j] / steps[i].new_bh_mass[j]);
//  }
  //}
  //fprintf(stderr, "chi2_merge=%f\n", chi2_merge); 
  //chi2 += chi2_merge;
  //if (isfinite(chi2) && chi2 > 0) chi2 = chi2_merge;

  //double chi2_csfr = 0;
  //chi2_csfr = pow(log10(steps[24].observed_cosmic_sfr / steps[4].observed_cosmic_sfr) / 0.5, 2);
  //chi2 += chi2_csfr;
  //fprintf(stderr, "chi2_csfr=%f\n", chi2_csfr);
  
  //double chi2_sm = 0;
  //chi2_sm += steps[24].log_sm[19] > 9? 0 : 100 * exp10(9 - steps[24].log_sm[19]);
  //chi2_sm += steps[24].log_sm[24] > 10? 0 : 100 * exp10(10 - steps[24].log_sm[24]);
  //chi2 += chi2_sm;
  //fprintf(stderr, "chi2_sm=%f\n", chi2_sm);

  //fprintf(stderr, "chi2_bhmf=%f\n", chi2_bhmf);

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
          // gsl_spline_free(steps[i].spline_uv); //free the spline object since we are allocating it every time.
          // if (steps[i].alloc2uvlf) gsl_spline_free(steps[i].spline_uv2); //same for the second segment of the UV-HM relation.
          // // gsl_spline_free(stseps[i].spline_std_uv);

        }

        if (steps[i].flag_alloc == 2)
        {
          gsl_spline_free(steps[i].spline_uv); //free the spline object since we are allocating it every time.
          if (steps[i].alloc2uvlf) gsl_spline_free(steps[i].spline_uv2); //same for the second segment of the UV-HM relation.
          // gsl_spline_free(stseps[i].spline_std_uv);
        }
    //if ((1 / steps[i].scale - 1 > 6) && (steps[i].cosmic_sfr <= steps[i-1].cosmic_sfr))
     //{
       //double ratio = steps[i-1].cosmic_sfr / steps[i].cosmic_sfr;
       //chi2 += ratio * ratio;
       //chi2_csfr += ratio * ratio;
     //}
  }
   //fprintf(stderr, "chi2_csfr: %.6f\n", chi2_csfr);

  // printf("The chi2 calculated by all_smf_chi2_err is: %f\n", chi2);
  //if (chi2_bhmf == 0 && isfinite(chi2))
//{
//  for (i = 0; i < NUM_PARAMS; i++) fprintf(stderr, "%.12f ", test.params[i]);
//  fprintf(stderr, "\n");
//}

  return (chi2);
}

int mcmc_smf_fit(int length, int run) {
  int i, j, dups=0;
  float last_chi2, chi2, first_chi2;
  struct smf_fit last_smf_fit=initial_smf_fit, cur_smf_fit=initial_smf_fit;
  FILE *out_fh = (run > 1) ? stdout : stderr;
  first_chi2 = last_chi2 = all_smf_chi2_err(initial_smf_fit);
  double chi2_repeat = last_chi2;
  int64_t count_repeat = 0;
  int64_t repeated = 0;

  for (i=1; i<length; i++) {
    repeated = 0;
    if (((count_repeat >= 100) && !run)) {
      fprintf(stderr, "before resetting, i=%d\n", i);
    // if (((i==100 || count_repeat >= 100) && last_chi2==first_chi2 && !run)) {
      // for (j=0; j<NUM_PARAMS; j++) eigenvalues[j]*=1.01; //it was 0.5, but I think if we got no luck in a certain region of parameter space, we should step out of it and explore more possibilities.
      for (j=0; j<NUM_PARAMS; j++) eigenvalues[j] = eigenvalues_original[j];
      // for (j=0; j<NUM_PARAMS; j++) fprintf(stderr, "eigenvalues[%d]=%e ", j, eigenvalues[j]);
      clear_stats();
      count_repeat =0;
      if (i == 100) 
      {
        num_points = dups = i=0;
      }
      // num_points = dups = i=0;
      //fprintf(stderr, "after resetting, i=%d\n", i);
      set_identity_portion();
      continue;
    }
    random_step(&cur_smf_fit);
    CHI2(cur_smf_fit) = chi2 = all_smf_chi2_err(cur_smf_fit);
    if (!isfinite(chi2) || chi2<0 ||
	(chi2>last_chi2 && drand48()>exp(0.5*inv_temperature*(last_chi2-chi2))))
    {
      //fprintf(stderr, "new chi2: %e, old chi2: %e\n", chi2, last_chi2); //used to see whether we are having reasonable new steps (chi2 > 0), or the new step doesn't make sense at all.
      chi2 = last_chi2;
      cur_smf_fit = last_smf_fit;
      dups++;
      if (chi2_repeat == last_chi2) repeated = 1;
    }
    chi2_repeat = last_chi2 = chi2;
    last_smf_fit = cur_smf_fit;
    add_to_stats(&cur_smf_fit);
    num_points++;
    //fprintf(stderr, "Repeated? %d\n", repeated);
    count_repeat = repeated ? count_repeat + 1 : 0;
    //fprintf(stderr, "Number of repetitions? %d\n", count_repeat);
    if (run > 0) stats_to_step(num_points);
    for (j=0; j<NUM_PARAMS+1; j++) fprintf(out_fh, "%.10f ", cur_smf_fit.params[j]);
    fprintf(out_fh, "%f\n", CHI2(cur_smf_fit));
  }
  initial_smf_fit = cur_smf_fit;
  return dups;
}

void init_mcmc_from_args(int argc, char **argv)
{
  int i;
  r250_init(87L);

  // if (argc < 7 || argc > 10000) {
  //   printf("Usage: %s model matching_scatter obs_scatter systematics mf_cache file1 file2 ...\n", argv[0]);
  //   exit(1);
  // }

  if (argc < 10 || argc > 10000) {
    printf("Usage: %s model matching_scatter obs_scatter systematics nonlinear_luminosity? no_bh_mergers? vel_dispersion? mf_cache file1 file2 ...\n", argv[0]);
    exit(1);
  }


  i = atol(argv[1]);
  if (i == 0) no_z_scaling = 1;
  if (i == 2) no_sfr_constraint = 1;
  if (i == 3) dplaw_only = 1;

  if (!atol(argv[2])) { no_matching_scatter = 1; }
  setup_psf(atoi(argv[3]));
  if (!atol(argv[3])) no_obs_scatter = 1;
  if (!atol(argv[4])) { no_systematics = 1; }
  // load_mf_cache(argv[5]);
  nonlinear_luminosity = atol(argv[5]);
  //  if (atol(argv[5])) { nonlinear_luminosity = 1; }
  if (atol(argv[6])) { no_bh_mergers = 1; }
  if (atol(argv[7])) { vel_dispersion = 1; }
  load_mf_cache(argv[8]);
  init_timesteps();

  num_obs_smfs = 0;
  for (i=9; i<argc; i++) {
    if (!strcmp(argv[i], "nopoints")) continue;
    load_real_smf(&obs_smfs, &num_obs_smfs, argv[i]);
    if (obs_smfs[num_obs_smfs-1].z_high > z_limit) z_limit = obs_smfs[num_obs_smfs-1].z_high;
  }
}

void set_identity_portion(void) {
  int64_t i;
  for (i=0; i<NUM_PARAMS; i++) {
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
  gsl_set_error_handler_off();
  
  init_mcmc_from_args(argc, argv);
  initial_conditions();
  fprintf(stderr, "Total number of data points: %ld", num_obs_smfs);
  return 0;

  srand(time(0));
  srand48(time(0));
  init_genrand64((unsigned long long) time(0));

  //for (i=0; i < num_outputs; i++)
  //{
    //printf("i=%d, aexp=%f, z=%f\n", i, steps[i].scale, 1 / steps[i].scale - 1);
  //}

  //return 0;

  clear_stats();
  set_identity_portion();
  for (i=0; i<num_temps; i++) {
    inv_temperature = temps[i];
    fprintf(stderr, "#Temp now: %f\n", 1.0/inv_temperature);
    mcmc_smf_fit(BURN_IN_LENGTH*temps[i], 0);
    stats_to_step(num_points);
    if (i<num_temps-1) clear_stats();
  }
  //return 0;
  inv_temperature = 1;
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

void init_orth_matrix(void) {
  int i,j;
  struct smf_fit address;
  for (i=0; i<NUM_PARAMS; i++)
    for (j=0; j<NUM_PARAMS; j++)
      orth_matrix[i][j] = (i==j) ? 1 : 0;

  if (!no_systematics) {
    i = &(MU(address)) - address.params;
    j = &(EFF_0(address)) - address.params;
    orth_matrix[i][j] = -1; //Build in anti-correlation between MU and EFF
  }
}

void print_orth_matrix(int type, int syst, float scatter, float z) {
  FILE *output;
  char buffer[1024];
  int i,j;
  sprintf(buffer, "mcmc_runs/orth_matrix_%d_%d_%f_%f.dat", type, syst, scatter, z);
  output = fopen(buffer, "w");
  for (i=0; i<NUM_PARAMS; i++) {
    fprintf(output, "Eig[%d]: %f\n", i, eigenvalues[i]);
  }
  fprintf(output, "\n");
  for (i=0; i<NUM_PARAMS; i++) {
    for (j=0; j<NUM_PARAMS; j++)
      fprintf(output, "%+.4f ", orth_matrix[i][j]);
    fprintf(output, "\n");
  }
  fclose(output);
}

void clear_stats(void) {
  int i, j;
  memset(&smf_sum, 0, sizeof(struct smf_fit));
  for (i=0; i<NUM_PARAMS; i++) {
    smf_sum.params[i] = 0;
    for (j=0; j<NUM_PARAMS; j++) {
      cov_matrix[i][j] = 0;
    }
  }
  num_points = 0;
}

void add_to_stats(struct smf_fit *a) {
  int i,j;
  for (i=0; i<NUM_PARAMS; i++) {
    smf_sum.params[i]+=a->params[i];
    for (j=0; j<NUM_PARAMS; j++) {
      cov_matrix[i][j] += a->params[i]*a->params[j];
    }
  }
}

void stats_to_step(int64_t num_stats) {
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
  for (i=0; i<NUM_PARAMS; i++) {
    if (eigenvalues[i] <=0 ) eigenvalues[i] = 0;
    eigenvalues[i] = sqrt(eigenvalues[i])*sd_mul;
  }
  //clear_stats();
}


void assert_model(struct smf_fit *a) {
  LAMBDA(*a) = LAMBDA_A(*a) = 0;
  //BETA(*a) = BETA_A(*a) = BETA_A2(*a) = 0;
  //DELTA_A(*a) = 
  DELTA_A2(*a) = 0;

  if (dplaw_only) {
    GAMMA(*a) = GAMMA_A(*a) = GAMMA_A2(*a) = 0;
    DELTA(*a) = 0.1;
  }

  //EFF_0_A3(*a) = 0;
  EXPSCALE(*a) = 4;
  //ALPHA_A2(*a) = 0;
  //DELTA_A2(*a) = 0;
  //GAMMA_A2(*a) = 0;

  if (no_bh_mergers) {
    BH_MERGE_F_0(*a) = 19;
    BH_MERGE_F_1(*a) = 0;
    BH_MERGE_W_0(*a) = 0.01;
    BH_MERGE_W_1(*a) = 0;
  }

  if (no_z_scaling) {
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
    //-1.740158237100 0.000000000000 0.000000000000 0.000000000000 0.000000000000 11.367077550400 0.000000000000 0.000000000000 -1.490795596100 0.000000000000 0.000000000000 0.000000000000 0.000000000000 0.000000000000 3.720329751000 0.000000000000 0.000000000000 0.394003775600 0.000000000000 0.000000000000 0.000000000000 0.000000000000 0.000000000000 0.000000000000 0.000000000000 0.160000000000 0.040000000000 1.500000000000 0.000000000000 0.000000000000 0.000000000000 1.394759967100 0.000000000000 0.000000000000 1.000000000000 0.234424754700 0.000000000000 1.765546427700
  }

  if (z_scaling_only) {
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

  if (no_matching_scatter) {
    SCATTER(*a) = 0;
  }
  if (no_obs_scatter) {
    SIGMA_Z(*a) = 0;
  }
  if (no_systematics) {
    //ICL_FRAC(*a) = 16;
    //ICL_FRAC_E(*a) = 0;
    KAPPA(*a) = 0;
    MU(*a) = 0;
    KAPPA_A(*a) = 0;
    MU_A(*a) = 0;
    BURST_DUST_AMP(*a) = 0;
    BURST_DUST_Z(*a) = 1;
    BURST_DUST(*a) = 0;
  }

  //BH_GAMMA_0(*a) = 1.1;
  //BH_SCATTER_0(*a) = 0.4;
  BH_SCATTER_1(*a) = 0;
  BH_ETA_CRIT_0(*a) = -1.5;
  //BH_ETA_CRIT_1(*a) = 0.0;
  //BH_EFFICIENCY_0(*a) = -1.0;
  ABHMF_SHIFT(*a) = 0.2;
  //EFF_0(*a) = 0.109; EFF_0_A(*a) = 3.441; EFF_0_A2(*a) = 5.079; EFF_0_A3(*a) = -0.781;
  //M_1(*a) = 2.151; M_1_A(*a) = 1.658; M_1_A2(*a) = 1.68; M_1_A3(*a) = -0.233;
  //ALPHA(*a) = -5.598; ALPHA_A(*a) = 20.731; ALPHA_A2(*a) = 13.455; ALPHA_A3(*a) = -1.321;
  //BETA(*a) = -1.911; BETA_A(*a) = -0.395; BETA_A2(*a) = 0.747;
  //GAMMA(*a) = -1.699; GAMMA_A(*a) = -4.206; GAMMA_A2(*a) = -0.809;
  DELTA(*a) = 0.055;
  //MU(*a) = 0.041;
  KAPPA_A(*a) = 0.0;
  // KAPPA(*a) = 0;
  // KAPPA_A(*a) = 0.155;
  //SCATTER(*a) = 0.23096395;
  //SCATTER_A(*a) = 0.05583049;

  //QM_0(*a) = 12.214; QM_1(*a) = 0.0592; QM_2(*a) = 0.139;
  //QWIDTH_0(*a) = 0.498; QWIDTH_1(*a) = 0.372; QWIDTH_2(*a) = 0.0124;

  BURST_DUST_AMP(*a) = 0;
  BURST_DUST_Z(*a) = 1;
  BURST_DUST(*a) = 0;
  SCATTER_A(*a) = 0;
  //SIGMA_Z(*a) = 0.071;
  SIGMA_A(*a) = 0;
  
  ICL_FRAC_Z(*a) = 0;
  ICL_FRAC_E(*a) = 0;
  // BURST_DUST(*a) = 0.823;
  // BURST_DUST_AMP(*a) = 0.273;
  // BURST_DUST_Z(*a) = 1.077;
  // RHO_05(*a) = 0.799;
}

void random_step(struct smf_fit *a) {
  int i,j,num_params=6;
  for (i=0; i<num_params; i++) { //NUM_PARAMS/2; i++) {
    //j = ind_var[i];
    //if (i == 5) j = 4;
    j = rand()%NUM_PARAMS;
    //if (j < 38 && j > 62)
    //{
    //    i--;
    //    continue;
    //}
    vector_madd(a->params, normal_random(0,eigenvalues[j]), orth_matrix[j]);
  }
  assert_model(a);
}
