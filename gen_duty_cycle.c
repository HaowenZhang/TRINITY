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

#define INTERP(y,x) { double s1, s2; s1 = (1.0-f)*steps[i].x[j]+f*steps[i+1].x[j];     \
    s2 = (1.0-f)*steps[i].x[j+1]+f*steps[i+1].x[j+1];			               \
    y = s1+mf*(s2-s1); }

#define LINTERP(y,x) { double s1, s2; s1 = (1.0-f)*log10(steps[i].x[j])+f*log10(steps[i+1].x[j]); \
    s2 = (1.0-f)*log10(steps[i].x[j+1])+f*log10(steps[i+1].x[j+1]);	\
    y = s1+mf*(s2-s1); }


#define LBOL_THRESH 44.0

double _frac_above_thresh(double thresh, double alpha, double delta, double norm);

int main(int argc, char **argv)
{
  int64_t i,j;
  struct smf_fit smf;
  if (argc<2+NUM_PARAMS) {
    fprintf(stderr, "Usage: %s mass_cache (mcmc output)\n", argv[0]);
    exit(1);
  }
  for (i=0; i<NUM_PARAMS; i++)
    smf.params[i] = atof(argv[i+2]);
gsl_set_error_handler_off();
  nonlinear_luminosity = 1;
  setup_psf(1);
  load_mf_cache(argv[1]);
  init_timesteps();
  INVALID(smf) = 0;
  calc_sfh(&smf);
  double t,m;
  printf("#1+z M_h SM M_bh duty_cycle dc_beta dc_mass_factor f(Lbol>%g) f(mdot>0.1) f(mdot>1) f(mdot>10)\n", LBOL_THRESH);
  for (t=0; t<num_outputs-1; t+=1.0/3.0) {
    i = t;
    double f = t-i;
    double zp1 = (1.0-f)/steps[i].scale + f/steps[i+1].scale;
    double dc = (1.0-f)*steps[i].smhm.bh_duty + f*steps[i+1].smhm.bh_duty;
    double dc_mbh = (1.0-f)*steps[i].smhm.dc_mbh + f*steps[i+1].smhm.dc_mbh;
    double dc_mbh_w = (1.0-f)*steps[i].smhm.dc_mbh_w + f*steps[i+1].smhm.dc_mbh_w;
    double bher_dist[BHER_BINS];
    //double efficiency = (1.0-f)*steps[i].smhm.bh_efficiency + f*steps[i+1].smhm.bh_efficiency;
    double alpha = (1.0-f)*steps[i].smhm.bh_alpha + f*steps[i+1].smhm.bh_alpha;
    double delta = (1.0-f)*steps[i].smhm.bh_alpha + f*steps[i+1].smhm.bh_delta;
    for (j=0; j<BHER_BINS; j++) {
      bher_dist[j] = (1.0-f)*steps[i].bher_dist[j] + f*steps[i+1].bher_dist[j];
      bher_dist[j] *= dc*BHER_INV_BPDEX;
    }
    double z = zp1 - 1;
    double mass_real = 13.5351-0.23712*z+2.0187*exp(-z/4.48394);
    for (m=8; m<15; m+=0.05) {
      if (m >= mass_real) continue;
      double mf = (m-M_MIN)*BPDEX+0.5;
      double sm;
      double dc = (1.0-f)*steps[i].smhm.bh_duty + f*steps[i+1].smhm.bh_duty;
      int64_t j = mf;
      mf -= j;
      
      double lbh_mass, log_bh_mass, bh_eta, log_bh_acc_rate, dc_obs;
      INTERP(lbh_mass,log_bh_mass);
      LINTERP(log_bh_mass,bh_mass_avg);
      INTERP(sm, sm_avg);
      LINTERP(log_bh_acc_rate,bh_acc_rate);
      INTERP(bh_eta,bh_eta);
      INTERP(dc_obs, dc_obs);
      if (!isfinite(log_bh_acc_rate)) continue;
      double f_lbol=0;

      double dc_mass_factor = exp((log_bh_mass - dc_mbh) / dc_mbh_w);
      dc_mass_factor = dc_mass_factor / (1 + dc_mass_factor);
      dc *= dc_mass_factor;

      for (j=0; j<BHER_BINS; j++) {
	double edd_frac = BHER_MIN+j*BHER_INV_BPDEX;
	double edd_r_obs = edd_frac + bh_eta;
	//double edd_r = edd_r_obs + log10(0.1/efficiency);
	double lir = -5.26 - 2.5*(edd_r_obs+lbh_mass); //=72.5 - 2.5(Lbol-7)
	double lbol = (lir - 72.5)/-2.5 + 7.0;
	double ff_lbol = (lbol > LBOL_THRESH) ? 1 : 0;
	if (!ff_lbol && (lbol+BHER_INV_BPDEX>LBOL_THRESH))
	  ff_lbol += (lbol+BHER_INV_BPDEX-LBOL_THRESH)*BHER_BPDEX;

	f_lbol += bher_dist[j]*ff_lbol;
      }

      double thresh_01 = -1.0 - bh_eta;
      double thresh_1 = 0 - bh_eta;
      double thresh_10 = 1.0 - bh_eta;
      // double sn = dc/doublePL_norm(alpha, delta, BHER_MIN, BHER_EFF_MAX, NULL);
      double sn = 1.0;
      double f_mdot01 = _frac_above_thresh(thresh_01, alpha, delta, sn);
      double f_mdot1 = _frac_above_thresh(thresh_1, alpha, delta, sn);
      double f_mdot10 = _frac_above_thresh(thresh_10, alpha, delta, sn);
      printf("%f %f %f %f %f %f %f %e %e %e %e %f %e %e\n", zp1, m, log10(sm), log_bh_mass, dc, dc_obs, dc_mass_factor, f_lbol, f_mdot01, f_mdot1, f_mdot10, bh_eta, bher_dist[(int64_t)((-bh_eta-BHER_MIN)*BHER_BPDEX)]*BHER_BPDEX, bher_dist[(int64_t)((-bh_eta-BHER_MIN)*BHER_BPDEX)]*BHER_BPDEX);
    }
  }
  return 0;
}


double _frac_above_thresh(double thresh, double alpha, double delta, double norm) {
  if (thresh < BHER_MIN) thresh = BHER_MIN;
  if (thresh >= BHER_EFF_MAX) return 0;
  //return (doublePL_norm(alpha, delta, thresh, BHER_EFF_MAX, NULL)*norm);
  return 1.0;
}
