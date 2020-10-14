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
#include "mah.h"
#include "check_syscalls.h"
#include "expcache2.h"
#include "sm_limits.h"

#define NUM_ZS 9
extern int64_t num_outputs;
extern struct timestep *steps;

double nd_to_mass(double scale, double nd) {
  double mass = 17;
  double last_tot = 0;
  double tot = 0;
  for (; tot < nd; mass -= 0.01) {
    last_tot = tot;
    tot += pow(10, mf_cache(scale, mass))/100.0;
  }
  mass += 0.015;
  return (mass - 0.01*(nd-last_tot)/(tot-last_tot)) ;
}


float _mar_from_mbins(int64_t n, int64_t j) {
  int64_t i;
  if (!n) return pow(10, M_MIN+(j+0.5)*INV_BPDEX)/steps[n].dt;
  if (n>=num_outputs-1) return _mar_from_mbins(num_outputs-2, j);
  if (j>=M_BINS-1) return 0;
  //Find out average progenitor mass:
  double sum = 0;
  double count = 0;
  for (i=0; i<M_BINS; i++) {
    if (!steps[n].mmp[i*M_BINS + j]) continue;
    sum += pow(10, M_MIN+(i+0.5)*INV_BPDEX)*steps[n].mmp[i*M_BINS+j];
    count += steps[n].mmp[i*M_BINS+j];
  }
  if (!count) return 0;
  sum /= count;
  return ((pow(10, M_MIN+(j+0.5)*INV_BPDEX) - sum)/steps[n].dt);
}

float mar_from_mbins(int64_t n, int64_t j) {
  float mar1 = _mar_from_mbins(n,j);
  float mar2 = _mar_from_mbins(n+1,j);
  return (0.5*(mar1+mar2));
}

float biterp (float a, float b, float c, float d, float f1, float f2) {
  float al = log10(a);
  float bl = log10(b);
  float cl = log10(c);
  float dl = log10(d);
  float e = al+f1*(bl-al);
  float f = cl+f1*(dl-cl);
  return (e+f2*(f-e));
}

int main(int argc, char **argv)
{
  double zs[NUM_ZS] = {0.0, 0.1, 1, 2, 3, 4, 5, 6, 7, 8};
  struct smf_fit the_smf;
  int64_t i, j, zi, k;
  float fm;
  double ft;
  FILE *ssfr_smah_f;
  FILE *sfe_smah_f;
  FILE *smhm_smah_f;
  char buffer[1024];
  double smah_ratio[16];
  double smah_norm[16];
  double smah_norm_sm[16];
  double smah_mass[16];
  double smah_started[16];

  if (argc<3+NUM_PARAMS) {
    fprintf(stderr, "Usage: %s mass_cache (mcmc output)\n", argv[0]);
    exit(1);
  }
  for (i=0; i<NUM_PARAMS; i++)
    the_smf.params[i] = atof(argv[i+2]);
  gsl_set_error_handler_off();
  gen_exp10cache();
  setup_psf(1);
  load_mf_cache(argv[1]);
  init_timesteps();
  INVALID(the_smf) = 0;
  calc_sfh(&the_smf);

  for (zi=0; zi<NUM_ZS; zi++) {
    float z = zs[zi];
    sprintf(buffer, "results/ssfr_smah_z%.1f.dat", z);
    ssfr_smah_f = check_fopen(buffer, "w");
    fprintf(ssfr_smah_f, "#Z SM10 SFR10 SM11 SFR11 SM12 SFR12 SM13 SFR13 SM14 SFR14 SM15 SFR15\n");
    sprintf(buffer, "results/sfe_smah_z%.1f.dat", z);
    sfe_smah_f = check_fopen(buffer, "w");
    fprintf(sfe_smah_f, "#Z M10 10 P10 M10.5 10.5 P10.5 M11 11 P11 M11.5 11.5 P11.5 M12 12 P12 M12.5 12.5 P12.5\n");
    sprintf(buffer, "results/smhm_smah_z%.1f.dat", z);
    smhm_smah_f = check_fopen(buffer, "w");
    fprintf(smhm_smah_f, "#Z M10 10 P10 M10.5 10.5 P10.5 M11 11 P11 M11.5 11.5 P11.5 M12 12 P12 M12.5 12.5 P12.5\n");

    calc_step_at_z(z, &i, &ft);
    float max_m = nd_to_mass(steps[i].scale, 1e-9);
    for (k=i; k>0; k--) {
      float min_sm = pow(10, sm_min_mass(1.0/steps[k].scale-1.0));
      float max_sm = pow(10, sm_max_mass(1.0/steps[k].scale-1.0));
      fprintf(ssfr_smah_f, "%f", 1.0/steps[k].scale-1.0);
      fprintf(sfe_smah_f, "%f", 1.0/steps[k].scale-1.0);
      fprintf(smhm_smah_f, "%f", 1.0/steps[k].scale-1.0);
      for (j=10; j<16; j++) {
	if (k==i) smah_started[j-10] = 0;
	float mnow = m_now(steps[i].scale, 10+(j-10));
	float m = m_at_a(mnow, steps[k].scale);
	fm = (m-M_MIN)*BPDEX - 0.5;
	if (fm < 0 || m>max_m) { fprintf(ssfr_smah_f, " 0 0"); fprintf(smhm_smah_f, " 0 0 0"); fprintf(sfe_smah_f, " 0 0 0"); continue; }
	int64_t mb = fm;
	fm -= mb;
	float mar1 = mar_from_mbins(k, mb);
	float mar2 = mar_from_mbins(k, mb+1);
	float mar = mar1 + fm*(mar2-mar1);
	float sm1 = log10(steps[k].sm_avg[mb]);
	float sm2 = log10(steps[k].sm_avg[mb+1]);
	float sm = exp10(sm1 + fm*(sm2-sm1));
	//sm = calc_sm_at_m(m, steps[k].smhm);
	//double scatter_corr = exp(pow(steps[k].smhm.scatter*log(10), 2)/2.0);
	//sm = pow(10, sm)*scatter_corr;
	float sfr1 = log10(steps[k].sfr[mb]);
	float sfr2 = log10(steps[k].sfr[mb+1]);
	float sfr = exp10(sfr1 + fm*(sfr2-sfr1));
	//sfr = calc_ssfr(log(sm)/log(10), 1.0/steps[k].scale-1.0)*sm;
	//float ratio = sfr/sm/mar*pow(10, m);
	//double scatter_corr = exp(pow(0.3*log(10), 2)/2.0);
	mar = ma_rate_avg(mnow, steps[k].scale);//*scatter_corr;
	float ratio = sfr/sm/mar*pow(10, m);
	if (k==i || !smah_started[j-10]) { smah_ratio[j-10] = ratio; }
	if (!(sm > min_sm && sm < max_sm)) ratio *= -1.0;
	fprintf(ssfr_smah_f, " %.6e %.6e", sm, sfr);
	fprintf(smhm_smah_f, " %.2e", sm/pow(10, m));
	ratio = sfr/(0.17*mar);
	if (k==i || !smah_started[j-10]) { smah_norm[j-10] = ratio; smah_norm_sm[j-10] = sm; }
	if (k==i || !smah_started[j-10]) { smah_mass[j-10] = m; }
	if (!(sm > min_sm && sm < max_sm)) ratio *= -1.0;
	fprintf(sfe_smah_f," %.2f %.2e", m, ratio);
	float pred = smah_norm[j-10]*pow(10, (smah_ratio[j-10]-1.0)*(m-smah_mass[j-10]));
	float pred_sm = log10(smah_norm_sm[j-10] * pow(10, smah_ratio[j-10]*(m-smah_mass[j-10])));
	fprintf(smhm_smah_f," %.2f %.2e", m, pow(10, pred_sm - m));
	if (!(sm > min_sm && sm < max_sm)) pred *= -1.0;
	fprintf(sfe_smah_f, " %.2e", pred);
	if (sm > min_sm && sm < max_sm) smah_started[j-10] = 1;
      }
      fprintf(smhm_smah_f, "\n");
      fprintf(ssfr_smah_f, "\n");
      fprintf(sfe_smah_f, "\n");
    }
    fclose(ssfr_smah_f);
    fclose(sfe_smah_f);
  }
  return 0;
}
