#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "observations.h"
#include "smf.h"
#include "all_smf.h"
#include "distance.h"
#include "integrate.h"
#include "mlist.h"
#include "calc_sfh.h"
#include "expcache2.h"
#include "universe_time.h"
//#include "fitter.h"

#define MSUN_YR_TO_ERG_S 5.6631e46

#define INTERP(y,x) { double s1, s2; s1 = (1.0-f)*steps[i].x[j]+f*steps[i+1].x[j];     \
    s2 = (1.0-f)*steps[i].x[j+1]+f*steps[i+1].x[j+1];			               \
    y = s1+mf*(s2-s1); }

#define LINTERP(y,x) { double s1, s2; s1 = (1.0-f)*log10(steps[i].x[j])+f*log10(steps[i+1].x[j]); \
    s2 = (1.0-f)*log10(steps[i].x[j+1])+f*log10(steps[i+1].x[j+1]);	\
    y = s1+mf*(s2-s1); }

float fitter(float *params) {
  struct smf_fit test;
  int i;
  for (i=0; i<NUM_PARAMS; i++)
    test.params[i] = params[i];

  assert_model(&test);

  for (i=0; i<NUM_PARAMS; i++)
    params[i] = test.params[i];
  //iterations++;
  float err = all_smf_chi2_err(test);
  if (!isfinite(err) || err<0) return 1e30;
  return err;
}

float calc_chi2(float *params) {
  return fitter(params);
}

int main(int argc, char **argv)
{
  int64_t i, j, k;
  struct smf_fit smf;
  if (argc<2+NUM_PARAMS) {
    fprintf(stderr, "Usage: %s mass_cache (mcmc output)\n", argv[0]);
    exit(1);
  }
  for (i=0; i<NUM_PARAMS; i++)
    smf.params[i] = atof(argv[i+2]);
  assert_model(&smf);
  gsl_set_error_handler_off();
  nonlinear_luminosity = 1;
  setup_psf(1);
  load_mf_cache(argv[1]);
  init_timesteps();
  INVALID(smf) = 0;
  //double chi2 = calc_chi2(smf.params);
  //printf("Actual chi2=%e\n", chi2);
  calc_sfh(&smf);
  printf("#Is the model invalid? %e\n", INVALID(smf));
  double t,m;
  // int t;
  printf("#z rho_tot rho_rad rho_rad_6m7, rho_rad_7m8, rho_rad_8m9, rho_rad_9m10, rho_kin\n");

  double mbh_lims[] = {6,7,8,9,10};
  //double rho_rad_split[4] = {0};

  for (i=0; i<num_outputs; i++) 
  {
    double z = 1.0 / steps[i].scale - 1;
    double rho_rad = 0, rho_tot = 0, rho_kin = 0;
    double rho_rad_split[4] = {0};
    double eff_rad = steps[i].smhm.bh_efficiency_rad;
    double mbh_min = steps[i].bh_mass_min, mbh_max = steps[i].bh_mass_max;
    double mbh_inv_bpdex = (mbh_max - mbh_min) / MBH_BINS;
    double scatter_tot = sqrt(steps[i].smhm.bh_scatter*steps[i].smhm.bh_scatter+
			steps[i].smhm.scatter*steps[i].smhm.scatter*steps[i].smhm.bh_gamma*steps[i].smhm.bh_gamma);
    double corr = exp(0.5 * pow((scatter_tot * M_LN10), 2));
    double eta_crit = steps[i].smhm.bh_eta_crit;
    // double mass_real = 13.5351-0.23712*z+2.0187*exp(-z/4.48394);

    // for (j=0; j<M_BINS; j++)
    // {
    //   double bhar = steps[i].bh_acc_rate[j] > 0 ? steps[i].bh_acc_rate[j] : 1e-8;
    //   rho_tot += eff_rad * bhar * corr * MSUN_YR_TO_ERG_S * steps[i].t[j];
    // }

    for (j=0; j<MBH_BINS; j++) 
    {
      double rho_rad_mbh = 0;
      double mbh = mbh_min + (j + 0.5) * mbh_inv_bpdex;
      for (k=0; k<LBOL_BINS; k++)
      {
        double lbol = LBOL_MIN + (k + 0.5) * LBOL_INV_BPDEX;
        rho_rad_mbh += exp10(lbol) * steps[i].lum_func_full[j*LBOL_BINS + k] * LBOL_INV_BPDEX;

        double eta_rad = lbol - 38.1 - mbh;
        if (eta_rad < eta_crit)
        {
          double eta_kin = log10(exp10(0.5*(eta_crit + eta_rad)) - exp10(eta_rad));
          double lkin = eta_kin + mbh + 38.1;
          if (eta_kin < -6 || mbh < 8) continue;
	  // double dlog_rad_dlog_kin = (exp10(0.5*eta_crit) - exp10(0.5*eta_rad)) / 
          //                           (0.5*exp10(0.5*eta_crit) - exp10(0.5*eta_rad));
          rho_kin += exp10(lkin) * steps[i].lum_func_full[j*LBOL_BINS + k] * LBOL_INV_BPDEX;
          fprintf(stderr, "scale=%f, mbh=%f, lbol=%f, eta_rad=%f, eta_kin=%f, lkin=%f, contribution=%e\n", steps[i].scale, mbh, lbol, eta_rad, eta_kin, lkin, exp10(lkin) * steps[i].lum_func_full[j*LBOL_BINS + k] * LBOL_INV_BPDEX);
	}

        

      }
      rho_rad += rho_rad_mbh;
      for (k=0; k<4; k++)
      {
        if ((mbh >= mbh_lims[k]) && (mbh < mbh_lims[k+1])) rho_rad_split[k] += rho_rad_mbh;
      }

    }
    rho_tot = rho_kin + rho_rad;
    // rho_kin = rho_tot - rho_rad;
    fprintf(stdout, "%f %e %e ", z, rho_tot, rho_rad);
    for (j=0; j<4; j++)
    {
      fprintf(stdout, "%e ", rho_rad_split[j]);
    }
    fprintf(stdout, "%e\n", rho_kin);
  }
  return 0;
}
