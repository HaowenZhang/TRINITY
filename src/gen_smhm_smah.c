#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
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
#include "universe_time.h"
#include "all_luminosities.h"

#define NUM_ZS 16
extern int64_t num_outputs;
extern struct timestep *steps;

#define NUM_LUMS 7

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


double _total_sm(int64_t n, int64_t j) {
  int64_t i;
  double sm=0;
  for (i=0; i<=n; i++) sm += steps[n].sm_hist[j*num_outputs + i];
  return sm;
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

void calc_lums(double lums[NUM_LUMS], float a, float smah_norm, float smah_ratio, float smah_mass, float m_start, float z_start, float sm_now) {
  float z0 = 1.0/a - 1.0;
  int64_t i,j;
  memset(lums, 0, sizeof(double)*NUM_LUMS);
  enum lum_types t_lums[NUM_LUMS] = { L_m1500, L_f125w, L_f160w, L_jwst_f200w, L_jwst_f277w, L_jwst_f356w, L_jwst_f444w };
  float t0 = scale_to_years(a);
  float m0 = m_evolution_avg(m_start, 1.0/(1.0+z_start), 1.0/(1.0+z0));
  for (i=0; i<70; i++) {
    float z1 = z0 + i/10.0;
    float z2 = z1 + 0.1;
    float t1 = scale_to_years(1.0/(1.0+z1));
    float t2 = scale_to_years(1.0/(1.0+z2));
    float m1 = m0;
    float m2 = m_evolution_avg(m0, 1.0/(1.0+z1), 1.0/(1.0+z2));
    m0 = m2;
    float sm1 = smah_norm * pow(10, smah_ratio*(m1-smah_mass));
    float sm2 = smah_norm * pow(10, smah_ratio*(m2-smah_mass));
    float Z = metallicity(log10(0.5*(sm1+sm2)*0.7), 0.5*(z1+z2));
    for (j=0; j<NUM_LUMS; j++) {
      float z_obs = z0;
      if (!j) z_obs = 0;
      lums[j] += (sm1-sm2)*lum_at_age_Z(&tlum, t_lums[j], t0-t1, t0-t2, Z, z_obs);
    }
  }
  for (j=0; j<NUM_LUMS; j++) {
    float df = (!j) ? 1e-18 : distance_factor(z0, 0);
    lums[j] = -2.5*log10(lums[j]*df);
    fprintf(stderr, "Lum for %f at %f (sm=%f): %f (%s)\n", m_start, z0, sm_now, lums[j], lum_names[t_lums[j]]);
  }
}

#define LF_START 20
#define LF_STOP 32
#define LF_BPMAG 10
#define LF_BINS ((LF_STOP-LF_START)*LF_BPMAG+1)

int main(int argc, char **argv)
{
//   double zs[NUM_ZS] = {0.1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
//   struct smf_fit the_smf;
//   int64_t i, j, zi;
//   float fm;
//   double ft;
//   FILE *smhm_smah_f, *direct_ssfr_f, *csfr_f;
//   char buffer[1024];
//   double smah_ratio[M_BINS];
//   double smah_norm[M_BINS];
//   double smah_mass[M_BINS];
//   double smah_started[M_BINS];
//   double csfrs[M_BINS];
//   double masses[M_BINS];
// #define SMF_MIN 6
// #define SMF_MAX 12
// #define SMF_BPDEX 10
// #define SMF_BINS ((SMF_MAX-SMF_MIN)*SMF_BPDEX+2)
//   double smf[SMF_BINS];
//   double smf_nevol[SMF_BINS];
//   double lf[NUM_LUMS][LF_BINS];
//   if (argc<4+NUM_PARAMS) {
//     fprintf(stderr, "Usage: %s mass_cache starting_z (mcmc output)\n", argv[0]);
//     exit(1);
//   }

//   double starting_z = atof(argv[2]);
//   for (i=0; i<NUM_PARAMS; i++)
//     the_smf.params[i] = atof(argv[i+3]);

//   gen_exp10cache();
//   setup_psf(1);
//   load_mf_cache(argv[1]);
//   init_timesteps();
//   INVALID(the_smf) = 0;
//   calc_sfh(&the_smf);
//   load_luminosities("fsps/mags_ab_zall_jwst.dat");

//   sprintf(buffer, "results/direct_ssfr_z%.1f.dat", starting_z);
//   direct_ssfr_f = check_fopen(buffer, "w");
//   fprintf(direct_ssfr_f, "#SM SSFR O_SSFR SHMAR\n");
//   calc_step_at_z(starting_z, &i, &ft);
//   for (j=0; j<M_BINS; j++) {
//     float sfr = steps[i].sfr[j];
//     if (!sfr) continue;
//     float sm = steps[i].sm_avg[j];
//     float tsm = _total_sm(i,j);
//     float mnow = m_now(steps[i].scale, M_MIN+(j+0.5)*INV_BPDEX);
//     float mar = ma_rate_avg(mnow, steps[i].scale);
//     fprintf(direct_ssfr_f, "%f %f %f %f %f\n", log10(sm), log10(sfr/tsm), log10(calc_ssfr(log10(sm), starting_z)), log10(mar)-(M_MIN+(j+0.5)*INV_BPDEX), M_MIN+(j+0.5)*INV_BPDEX);
//   }
//   fclose(direct_ssfr_f);

//   sprintf(buffer, "results/csfr_guess_z%.1f.dat", starting_z);
//   csfr_f = check_fopen(buffer, "w");
//   fprintf(csfr_f, "#Z Total_CSFR Total_Obs_CSFR CSFR_NoSM<4e7 Median_M M_up M_dn\n");
//   int64_t initial_i;
//   calc_step_at_z(starting_z, &initial_i, &ft);
//   sprintf(buffer, "results/smhm_21_z1_%.1f.dat", starting_z);
//   FILE *smhm_21_f = check_fopen(buffer, "w");
//   fprintf(smhm_21_f, "#Z SM-HM SM HM ND M@ND Steve_SMHM\n");

//   for (j=0; j<M_BINS; j++) smah_started[j] = 0;
//   for (zi=0; zi<NUM_ZS; zi++) {
//     float z = zs[zi];
//     if (z < starting_z) continue;
//     sprintf(buffer, "results/smhm_smah_z1_%.1f_z2_%.1f.dat", starting_z, z);
//     smhm_smah_f = check_fopen(buffer, "w");
//     fprintf(smhm_smah_f, "#M SM/HM\n");
//     sprintf(buffer, "results/smf_smah_z1_%.1f_z2_%.1f.dat", starting_z, z);
//     FILE *smf_smah_f = check_fopen(buffer, "w");
//     fprintf(smf_smah_f, "#SM ND(pred) ND(non-evol)\n");
//     sprintf(buffer, "results/lf_smah_z1_%.1f_z2_%.1f.dat", starting_z, z);
//     FILE *lf_smah_f = check_fopen(buffer, "w");
//     fprintf(lf_smah_f, "#L ND(ABS_1500+45) ND(125W) ND(160W) ND(200W) ND(277W) ND(356W) ND(444W)\n");

//     memset(smf, 0, sizeof(double)*SMF_BINS);
//     memset(smf_nevol, 0, sizeof(double)*SMF_BINS);
//     for (j=0; j<NUM_LUMS; j++)
//       memset(lf[j], 0, sizeof(double)*LF_BINS);
    
//     calc_step_at_z(z, &i, &ft);
//     float max_m = nd_to_mass(steps[i].scale, 1e-9);
//     //float m_thresh = nd_to_mass(steps[i].scale, 1e-5);
//     float min_sm = sm_min_mass(1.0/steps[i].scale-1.0);
//     float max_sm = sm_max_mass(1.0/steps[i].scale-1.0);
//     fprintf(smhm_smah_f, "#Min: %f; Max: %f\n", min_sm, max_sm);
//     double csfr = 0, csfr_lim = 0, csfr_llim = 0, csfr_ulim = 0, csfr_lulim=0, nd_21=0, lavg_smhm=0, lavg_m=0, lavg_sm=0;
//     for (j=0; j<M_BINS; j++) {
//       float mnow = m_now(1.0/(1.0+starting_z), M_MIN+j*INV_BPDEX);
//       float mnow_u = m_now(1.0/(1.0+starting_z), M_MIN+(j+0.5)*INV_BPDEX);
//       float mnow_d = m_now(1.0/(1.0+starting_z), M_MIN+(j-0.5)*INV_BPDEX);
//       float midscale = steps[i].scale;// + 0.2*(1.0/(1.0+starting_z)-steps[i].scale);
//       float m = m_at_a(mnow, midscale);
//       float dm = m_at_a(mnow_u, midscale) - m_at_a(mnow_d, midscale);
//       m = m_evolution_avg(M_MIN+j*INV_BPDEX, 1.0/(1.0+starting_z), steps[i].scale);
//       dm = m_evolution_avg(M_MIN+(j+0.5)*INV_BPDEX, 1.0/(1.0+starting_z), steps[i].scale) - m_evolution_avg(M_MIN+(j-0.5)*INV_BPDEX, 1.0/(1.0+starting_z), steps[i].scale);

//       dm *= BPDEX;
//       fm = (m-M_MIN)*BPDEX-0.5;
//       if (fm < 0 || m>max_m) { continue; }
//       int64_t mb = fm;
//       fm -= mb;
//       float mar1 = mar_from_mbins(i, mb);
//       float mar2 = mar_from_mbins(i, mb+1);
//       float mar = mar1 + fm*(mar2-mar1);
//       float sm1 = _total_sm(i,mb); //steps[i].sm_avg[mb];
//       float sm2 = _total_sm(i,mb+1); //steps[i].sm_avg[mb+1];
//       float tsm = sm1 + fm*(sm2-sm1);
//       sm1 = steps[i].sm_avg[mb];
//       sm2 = steps[i].sm_avg[mb+1];
//       float sm = sm1 + fm*(sm2-sm1);
//       /*float icl1 = steps[i].sm_icl[mb];
//       float icl2 = steps[i].sm_icl[mb+1];
//       float icl = icl1 + fm*(icl2-icl1);
//       if (icl > 0.2*sm) icl = 0.2*sm;*/
//       float sfr1 = steps[i].sfr[mb];
//       float sfr2 = steps[i].sfr[mb+1];
//       float sfr = sfr1 + fm*(sfr2-sfr1);
//       //sfr = calc_ssfr(log10(sm), 1.0/(steps[i].scale)-1.0)*sm;
//       float n1 = steps[i].t[mb];
//       float n2 = steps[i].t[mb+1];
//       float n = (n1+(n2-n1)*fm)*dm;
//       //sfr = calc_ssfr(log10(sm), 1.0/steps[i].scale-1.0)*sm;
//       //float ssfr = 5.96e-11*pow(10, -0.35*(log10(sm)-11.03))*exp(-pow(10, log10(sm)-11.03));
//       //sfr = ssfr*sm;
//       //double scatter_corr = exp(pow(0.3*log(10), 2)/2.0);
    
//       //mar = ma_rate_avg(mnow, steps[i].scale);//*scatter_corr;
//       mar = ma_rate_avg_mnow(m, steps[i].scale);//*scatter_corr;

//       float ratio = sfr/tsm/mar*pow(10, m);
//       if (!smah_started[j]) { smah_ratio[j] = ratio; }
//       double sm_scatter_corr = exp(pow(steps[i].smhm.scatter*log(10), 2)/2.0);

//       if (!smah_started[j]) { smah_norm[j] = sm/sm_scatter_corr; }
//       if (!smah_started[j]) { smah_mass[j] = m; }
//       double pred = log10(smah_norm[j] * pow(10, smah_ratio[j]*(m-smah_mass[j])));
//       double pred_sfr = smah_ratio[j]*pow(10, pred)*(mar / pow(10, m))*tsm/sm*sm_scatter_corr;
//       double lums[NUM_LUMS];
//       calc_lums(lums, steps[i].scale, smah_norm[j]*tsm/sm, smah_ratio[j], smah_mass[j], M_MIN+j*INV_BPDEX, starting_z, pred);
//       csfr += pred_sfr*n;
//       int64_t k, l;
//       struct smf c = smhm_at_z(8, the_smf);
//       double pred_nevol = calc_sm_at_m(m, c);
//       for (k=0; k<SMF_BINS; k++) {
// 	double sbin = SMF_MIN+(double)k/(double)SMF_BPDEX;
// 	double dm1 = (pred-sbin)/steps[i].smhm.scatter;
// 	double w1 = exp(-0.5*dm1*dm1)/(sqrt(2.0*M_PI)*steps[i].smhm.scatter);
// 	smf[k] += w1*n;
// 	double dm2 = (pred_nevol-sbin)/steps[i].smhm.scatter;
// 	double w2 = exp(-0.5*dm2*dm2)/(sqrt(2.0*M_PI)*steps[i].smhm.scatter);
// 	smf_nevol[k] += w2*n;
//       }

//       for (l=0; l<NUM_LUMS; l++) {
// 	for (k=0; k<LF_BINS; k++) {
// 	  double lbin = LF_START + (double)k/(double)LF_BPMAG;
// 	  double dm = (lums[l]-lbin)/(steps[i].smhm.scatter*2.5);
// 	  if (!isfinite(dm)) continue;
// 	  double w1 = exp(-0.5*dm*dm)/(sqrt(2.0*M_PI)*steps[i].smhm.scatter*2.5);
// 	  lf[l][k] += n*w1;
// 	}
//       }

//       csfrs[j] = pred_sfr*n;
//       masses[j] = m;
//       /*if (pred_sfr > 0) {
// 	float obs_fraction = 0.5*(1.0-erf((log10(0.4)-log10(pred_sfr/sm_scatter_corr))/(sqrt(2)*steps[i].smhm.scatter)));
// 	csfr_lim += pred_sfr*n*obs_fraction;
// 	if (sm>4e7) csfr_llim += pred_sfr*n*obs_fraction;
// 	obs_fraction = 0.5*(1.0-erf((log10(2.0)-log10(pred_sfr/sm_scatter_corr))/(sqrt(2)*steps[i].smhm.scatter)));
// 	if (sm>1e8) csfr_ulim += pred_sfr*n*obs_fraction;
// 	obs_fraction = 0.5*(1.0-erf((log10(0.1)-log10(pred_sfr/sm_scatter_corr))/(sqrt(2)*steps[i].smhm.scatter)));
// 	if (sm>1e7) csfr_lulim += pred_sfr*n*obs_fraction;
// 	}*/
      
//       if (pred_sfr > 0) {
// 	float uv_lum = lums[0]-45+2.5*log10(sm_scatter_corr);
// 	float scatter = sqrt(2)*steps[i].smhm.scatter*2.5;
// 	float obs_fraction = 0.5*(1.0-erf((uv_lum+18)/scatter));
// 	float obs_21 = 0.5*(1.0-erf((uv_lum+21)/scatter));
// 	float pred_21 = pred;
// 	if (uv_lum + scatter > -21) pred_21 += (uv_lum + scatter + 21.0)/2.5/2.0;
// 	nd_21 += n*obs_21;
// 	lavg_smhm += n*obs_21*(pred_21-m);
// 	lavg_m += n*obs_21*m;
// 	lavg_sm += n*obs_21*(pred_21);
// 	csfr_lim += pred_sfr*n*obs_fraction;
// 	obs_fraction = 0.5*(1.0-erf((uv_lum+17.7)/scatter));
//         csfr_llim += pred_sfr*n*obs_fraction;
// 	obs_fraction = 0.5*(1.0-erf((uv_lum+19)/scatter));
// 	csfr_ulim += pred_sfr*n*obs_fraction;
// 	obs_fraction = 0.5*(1.0-erf((uv_lum+17)/scatter));
// 	csfr_lulim += pred_sfr*n*obs_fraction;
//       }
//       if ((pred > min_sm-1.5 && pred < max_sm))
// 	fprintf(smhm_smah_f,"%.2f %.2f # %e\n", m, pred-m, smah_ratio[j]);
//       smah_started[j] = 1;
//     }
//     double m_med=0, csfr2=0, m_up=0, m_dn=0;

//     for (j=0; j<M_BINS; j++) {
//       csfr2 += csfrs[j];
//       if (csfr2 > csfr/2.0 && !m_med) {
// 	m_med = masses[j] + (csfr/2.0-csfr2)*(masses[j]-masses[j-1])/csfrs[j];
//       }
//       if (csfr2 > csfr/4.0 && !m_dn) {
// 	m_dn = masses[j] + (csfr/4.0-csfr2)*(masses[j]-masses[j-1])/csfrs[j];
//       }
//       if (csfr2 > csfr*0.75 && !m_up) {
// 	m_up = masses[j] + (csfr*0.75-csfr2)*(masses[j]-masses[j-1])/csfrs[j];
//       }
//     }
//     m_up -= m_med;
//     m_dn = m_med - m_dn;
//     fprintf(csfr_f, "%f %f %f %f %f %f %f %f %f\n", 1.0/steps[i].scale-1.0, log10(csfr), log10(csfr_lim), log10(csfr_llim), m_med, m_up, m_dn, log10(csfr_ulim), log10(csfr_lulim));
//     if (nd_21 > 0) {
//       fprintf(smhm_21_f, "%f %f %f %f %e %f %f\n", 1.0/steps[i].scale - 1.0, lavg_smhm/nd_21, lavg_sm/nd_21, lavg_m/nd_21, nd_21, nd_to_mass(steps[i].scale, nd_21), lavg_sm/nd_21 - nd_to_mass(steps[i].scale, nd_21));
//     }
//     fclose(smhm_smah_f);
//     for (j=0; j<SMF_BINS; j++) {
//       double sbin = SMF_MIN+(double)j/(double)SMF_BPDEX;
//       fprintf(smf_smah_f, "%f %e %e\n", sbin, smf[j], smf_nevol[j]);
//     }
//     int64_t k;
//     for (k=0; k<LF_BINS; k++) {
//       double sbin = LF_START + (double)k/(double)LF_BPMAG;
//       fprintf(lf_smah_f, "%f", sbin);
//       for (j=0; j<NUM_LUMS; j++) fprintf(lf_smah_f, " %e", lf[j][k]);
//       fprintf(lf_smah_f, "\n");
//     }
//     fclose(smf_smah_f);
//     fclose(lf_smah_f);
//   }
//   fclose (smhm_21_f);
//   fclose(csfr_f);
  return 0;
}
