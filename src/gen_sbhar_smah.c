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
  FILE *mh_smah_f;
  FILE *ssfr_smah_f;
  FILE *sfe_smah_f;
  FILE *smhm_smah_f;
  char buffer[1024];
  double smah_ratio[16];
  double smah_norm[16];
  double smah_norm_bhm[16];
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
    sprintf(buffer, "results/mh_smah_z%.1f.dat", z);
    mh_smah_f = check_fopen(buffer, "w");
    fprintf(mh_smah_f, "#Z Mh10 fm10 mb10 Mh11 fm11 mb11 Mh12 fm12 mb12 Mh13 fm13 mb13 Mh14 fm14 mb14 Mh15 fm15 mb15\n");
    sprintf(buffer, "results/sbhar_smah_z%.1f.dat", z);
    ssfr_smah_f = check_fopen(buffer, "w");
    fprintf(ssfr_smah_f, "#Z Mbh10 BHAR10 BHMR10 BHER10 Mbh11 BHAR11 BHMR11 BHER11 Mbh12 BHAR12 BHMR12 BHER12 Mbh13 BHAR13 BHMR13 BHER13 Mbh14 BHAR14 BHMR14 BHER14 Mbh15 BHAR15 BHMR15 BHER15 \n");
    sprintf(buffer, "results/bhfe_smah_z%.1f.dat", z);
    sfe_smah_f = check_fopen(buffer, "w");
    fprintf(sfe_smah_f, "#Z M10 10 P10 M11 11 P11 M12 12 P12 M13 13 P13 M14 14 P14 M15 15 P15\n");
    sprintf(buffer, "results/bhmhm_smah_z%.1f.dat", z);
    smhm_smah_f = check_fopen(buffer, "w");
    fprintf(smhm_smah_f, "#Z M10 10 P10 M11 11 P11 M12 12 P12 M13 13 P13 M14 14 P14 M15 15 P15\n");

    calc_step_at_z(z, &i, &ft);
    float max_m = nd_to_mass(steps[i].scale, 1e-9);
    for (k=i; k>0; k--) {
      float min_bhm = pow(10, 5);
      float max_bhm = pow(10, 10);
      fprintf(ssfr_smah_f, "%f", 1.0/steps[k].scale-1.0);
      fprintf(sfe_smah_f, "%f", 1.0/steps[k].scale-1.0);
      fprintf(smhm_smah_f, "%f", 1.0/steps[k].scale-1.0);
      fprintf(mh_smah_f, "%f", 1.0/steps[k].scale-1.0);
      for (j=10; j<16; j++) {// j labels the halo mass
	if (k==i) smah_started[j-10] = 0;
	float mnow = m_now(steps[i].scale, 10+(j-10)); //m_now is the descendant halo mass at z=0, of the halo that has the halo mass of 10 + (j - 10) at scale==steps[i].scale.
	float m = m_at_a(mnow, steps[k].scale); //m is the progenitor halo mass of the halos (with Mh=10+(j-10) at the i-th snapshot) at the k-th snapshot. We don't have a formula to get the mass mapping between two arbitrary snapshots, so we must use the halo mass at z=0 as a intermediate.
  // printf("input=%f, output=%f\n", 10+(j-10)*0.5, m_now(0.2, 10+(j-10)*0.5));
	fm = (m-M_MIN)*BPDEX - 0.5;
	if (fm < 0 || m>max_m) { fprintf(ssfr_smah_f, " %.6e %.6e %.6e %.6e", 0, 0, 0, 0); fprintf(smhm_smah_f, " %.6e %.6e %.6e", 0, 0, 0); fprintf(sfe_smah_f, " %.6e %.6e %.6e", 0, 0, 0); continue; }
	int64_t mb = fm;
	fm -= mb;
        fprintf(mh_smah_f, " %.6f %.6f %d", m, fm, mb);
	float mar1 = mar_from_mbins(k, mb);
	float mar2 = mar_from_mbins(k, mb+1);
	float mar = mar1 + fm*(mar2-mar1);
 
  //float bhm1 = steps[k].bh_mass_avg[mb];
  //float bhm2 = steps[k].bh_mass_avg[mb+1];
  //float bhm = bhm1 + fm*(bhm2-bhm1);
 
  float bhm1 = log10(steps[k].bh_mass_avg[mb]);
  float bhm2 = log10(steps[k].bh_mass_avg[mb+1]);
  float bhm = exp10(bhm1 + fm*(bhm2-bhm1));

  //float bhar1 = steps[k].bh_acc_rate[mb];
  //float bhar2 = steps[k].bh_acc_rate[mb+1];
  //float bhar = bhar1 + fm*(bhar2-bhar1);

  float bhar1 = log10(steps[k].bh_acc_rate[mb]);
  float bhar2 = log10(steps[k].bh_acc_rate[mb+1]);
  float bhar = exp10(bhar1 + fm*(bhar2-bhar1));

  float eff_rad = steps[k].smhm.bh_efficiency_rad;
  
  float eta_avg1 = bhar1 - bhm1 + log10(4.5e8 * eff_rad);
  float eta_avg2 = bhar2 - bhm2 + log10(4.5e8 * eff_rad);
  float eta_avg = exp10(eta_avg1 + fm*(eta_avg2 - eta_avg1));
  //float eta_avg = bhar / bhm * 4.5e8 * eff_rad;

  //float bhmr1 = steps[k].bh_merge_rate[mb];
  //float bhmr2 = steps[k].bh_merge_rate[mb+1];
  //float bhmr = bhmr1 + fm*(bhmr2-bhmr1);

  float bhmr1 = log10(steps[k].bh_merge_rate[mb]);
  float bhmr2 = log10(steps[k].bh_merge_rate[mb+1]);
  float bhmr = exp10(bhmr1 + fm*(bhmr2-bhmr1));


 
	mar = ma_rate_avg(mnow, steps[k].scale);//*scatter_corr;

  float ratio = bhar/bhm/mar*pow(10, m);
	if (k==i || !smah_started[j-10]) { smah_ratio[j-10] = ratio; }
  if (!(bhm > min_bhm && bhm < max_bhm)) ratio *= -1.0;
	fprintf(ssfr_smah_f, " %.6e %.6e %.6e %.6e", bhm, bhar, bhmr, eta_avg);
  fprintf(smhm_smah_f, " %.2e", bhm/pow(10, m));
	ratio = bhar/(0.17*mar);
  if (k==i || !smah_started[j-10]) { smah_norm[j-10] = ratio; smah_norm_bhm[j-10] = bhm; }
	if (k==i || !smah_started[j-10]) { smah_mass[j-10] = m; }
	//if (!(bhm > min_bhm && bhm < max_bhm)) ratio *= -1.0;
	fprintf(sfe_smah_f," %.2f %.2e", m, ratio);
	float pred = smah_norm[j-10]*pow(10, (smah_ratio[j-10]-1.0)*(m-smah_mass[j-10]));
  float pred_bhm = log10(smah_norm_bhm[j-10] * pow(10, smah_ratio[j-10]*(m-smah_mass[j-10])));
	fprintf(smhm_smah_f," %.2f %.2e", m, pow(10, pred_bhm - m));
  if (!(bhm > min_bhm && bhm < max_bhm)) pred *= -1.0;
	fprintf(sfe_smah_f, " %.2e", pred);
  if (bhm > min_bhm && bhm < max_bhm) smah_started[j-10] = 1;
      }
      fprintf(smhm_smah_f, "\n");
      fprintf(ssfr_smah_f, "\n");
      fprintf(sfe_smah_f, "\n");
      fprintf(mh_smah_f, "\n");
    }
    fclose(ssfr_smah_f);
    fclose(sfe_smah_f);
    fclose(mh_smah_f);
  }
  return 0;
}
