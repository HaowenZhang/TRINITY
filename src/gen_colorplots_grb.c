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
#include "all_luminosities.h"
#include "mt_rand.h"
#include "universe_time.h"
#include <string.h>
#include <gsl/gsl_sf.h>

#define NUM_ZS 9
extern int64_t num_outputs;
extern struct timestep *steps;

#define NUM_PLAWS 3
#define NUM_VELS 3
#define NUM_RADS 5
#define TMIN 50e6
#define NUM_L_THRESHES 100
#define FONG_THRESH 26
#define SMF_MIN 7
#define SMF_MAX 12
#define SMF_BPDEX 24
#define NUM_SMS ((SMF_MAX-SMF_MIN)*SMF_BPDEX+1)
#define SMF_MAX_Z 0.2
#define SFR_MIN -13
#define SFR_MAX -8
#define SFR_BPDEX 24
#define NUM_SFRS ((SFR_MAX-SFR_MIN)*SFR_BPDEX+1)
#define GRB_LF_MIN -1
#define GRB_LF_MAX 3
#define GRB_BPDEX 24
#define NUM_GRBS ((GRB_LF_MAX-GRB_LF_MIN)*GRB_BPDEX+1)

#define VELF_MIN 1
#define VELF_MAX 4
#define VELF_BPDEX 10
#define NUM_VELF ((VELF_MAX-VELF_MIN)*VELF_BPDEX+1)

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

float TOTAL_SM(int x, int y, int z) {
  float sm = steps[x].sm_hist[y*num_outputs+z];
  float icl = steps[x].icl_stars[y*num_outputs+z];
  if (!isfinite(sm)) sm = 0;
  if (!isfinite(icl)) icl = 0;
  return (sm+icl+1e-10);
}

float sfh_m_index(int j, float fm, int t, int k) {
  float ft = (float)(((int)t)%3) / 3.0;
  float fk = (float)(((int)k)%3) / 3.0;
  int i = t/3;
  int i2 = k/3;
  if (i2 >= i) fk=0;
  float smh1 = biterp(TOTAL_SM(i,j,i2), TOTAL_SM(i+1,j,i2),
		      TOTAL_SM(i,j+1,i2), TOTAL_SM(i+1,j+1,i2), 
		      ft, fm);
  float smh2 = biterp(TOTAL_SM(i,j,i2+1), TOTAL_SM(i+1,j,i2+1),
		      TOTAL_SM(i,j+1,i2+1), TOTAL_SM(i+1,j+1,i2+1), 
		      ft, fm);
  if (!isfinite(smh2)) smh2=smh1;
  if (!isfinite(smh1)) smh1=-100;
  float sfr = pow(10, smh1+fk*(smh2-smh1))/steps[i].dt;
  //  if (j>25 && j < 30) printf("SFR: %e\n", sfr);
  return sfr;
}

float instrument_completeness(float z, int weighting) {
  if (weighting > 3) return exp(-1.0*z); //(0.5 - 0.5*erf(z-3));
  //return (0.5-0.5*erf(1.0*(z-0.3)));
  //return (exp(-6.0*z)*pow(1.0+z,2));
  //return exp(-4.3*z);
  float c = exp(-4.3*z);
  if (c < 0.005) c = 0.005;
  return c;
  return 1;
}

float missed_hosts_fraction(float *lums, float m, float z, int64_t i, float thresh) {
  float lum = lum_at_hm_z(lums, m, z);
  return 0.5*(1.0+erf((lum-thresh)/(M_SQRT2*2.5*steps[i].smhm.scatter)));
}

float concentration(float hm, float z) {
  return (pow(10, (2.2358 - 0.10*hm))/(1.0+z));
}

float Delta_vir(float z) {
  float om = 0.27;
  float omega = om*(pow(1.0+z,3)/(om*pow(1.0+z, 3) + 1.0-om));
  float x = omega - 1.0;
  return ((1.0/(1.0+x)) * (18*M_PI*M_PI + 82*x - 39*x*x));
}

double x_to_d(double x) {
  double xp1 = 1.0+x;
  return ((log1p(x) - x/xp1)/(x*x*x));
}

double d_to_x(double d) {
  double x = d;
  double td = x_to_d(x);
  double new_x;
  while (fabs((d-td) / d) > 1e-7) {
    double td2 = x_to_d(x+0.1);
    double slope = (td2 - td) / 0.1;
    new_x = x + (d-td)/slope;
    if (new_x < 0) x/=2;
    else x = new_x;
    td = x_to_d(x);
  }
  return x;
}

float ratio_r200c_rvir(float cvir, float vir_dens, float dens_200c) {
  float d = x_to_d(cvir);
  return (d_to_x(d*dens_200c/vir_dens)/cvir);
}

float passive_frac(float sm, float z) {
  return 1.0/(pow(10, -1.3*(sm-(10.2+0.5*z)))+1.0);
}

float escape_velocity(float m, float z, float re_thresh) {
  float dvir = Delta_vir(z);
  float scale = 1.0/(1.0+z);
  float uvol = scale*scale*scale;
  float b_dens =  0.27*2.77519737e11/uvol;
  float g = 4.30117902e-9;
  float c = concentration(m, z);
  float rvir = pow((pow(10,m)/(4.0*M_PI*dvir*b_dens/3.0)),1.0/3.0);
  float c_dens = b_dens / (0.27/uvol / (0.27/uvol + 0.73));
  float dens_200c = 200*c_dens;
  float vir_dens = b_dens*dvir;
  float fr_start = (0.015)*ratio_r200c_rvir(c, vir_dens, dens_200c);
  float fr_end = re_thresh * fr_start;
  //float rs = rvir/c;
  //float rho = pow(10, m) / (4.0*M_PI*(pow(rs, 3))*(log(1.0+c)-1.0/(1.0+1.0/c)));
  float d = 1.0/(log(1.0+c)-c/(1.0+c));
  float phi = -d*g*pow(10, m)*(log(1.0+c*fr_start)/(rvir*fr_start)-log(1.0+c*fr_end)/(rvir*fr_end)); // - g*pow(10,m)/(rvir);
  float vesc = sqrt(2.0*fabs(phi));
  return vesc;
}

float re_obs(float m, float z, float vel) {
  float v0 = escape_velocity(m, z, 7);
  int64_t i = 0;
  float dre0 = 2.5;
  float re0 = 7;
  float vmax = escape_velocity(m,z,20);
  if (vel > vmax) return -1;
  while (fabs(v0-vel)>0.01 && (++i)<20) {
    float v1 = escape_velocity(m, z, re0+dre0);
    float re1 = re0 + (vel-v0)*dre0/(v1-v0);
    if (re1 > 20) { re0 += (20-re0)/2.0; }
    else if (re1 < 0) { re0 /= 2.0; }
    else { re0 = re1; }
    dre0 /= 2.0;
    v0 =  escape_velocity(m, z, re0);
  }
  return re0;
}

float maxwell_boltzmann_cdf(float mu, float v) {
  float a = mu*(sqrt(M_PI/2.0)/2.0);
  float cdf = erf(v/(sqrt(2)*a)) - sqrt(2.0/M_PI)*v*exp(-v*v/(2.0*a*a))/a;
  return cdf;
}

float exponential_cdf(float mu, float v) {
  return (1.0-exp(-v/mu));
}

float escape_fraction(float vel, float hm, float z, float re_thresh) {
  float ev = escape_velocity(hm, z, re_thresh);
  float ef = 1.0-exponential_cdf(vel, ev);
  return ef;
}

float plaw_weight(float dt, float ttot, float index) {
  float weight=1;
  if (index==-1.0) return log(1.0+dt/ttot);
  if (index == -1.5) weight = 1.383012e+02/7.309003e-03;
  else if (index == -0.5) weight = 1.383012e+02/4.348812e+06;
  return (weight*(1.0/(index+1.0))*(pow(ttot+dt, index+1.0)-pow(ttot, index+1.0)));
}

float band_fn(float e, float e0, float a, float b) {
  if (e < (a-b)*e0) return (pow(e,a)*exp(-e/e0));
  return (pow((a-b)*e0,a-b) * pow(e,b) * exp(b-a));
}

float int_band_fn(float e1, float e2) {
  double s = 0;
  int64_t i;
  for (i=0; i<100; i++) {
    double v = e1 + (i+0.5)*(e2-e1)/100.0;
    s += band_fn(v, 511, -1, -2.3);
  }
  return s*(e2-e1)/100.0;
}

#define BAND_E0 (511.0)
#define BAND_POWER_B (-2.3)

double exact_int_band_fn(double e1, double e2) {
  double s1 = 0, s2 = 0;
  if (e1 < BAND_E0) {
    double emax = (e2 > BAND_E0) ? BAND_E0 : e2;
    s1 = gsl_sf_expint_Ei(-emax/BAND_E0) - 
      gsl_sf_expint_Ei(-e1/BAND_E0);
  }
  if (e2 > BAND_E0) {
    double emin = (e1 < BAND_E0) ? BAND_E0 : e1;
    const double amb = -1.0 - BAND_POWER_B;
    s2 = (pow(amb*BAND_E0,amb)*exp(-amb)/(BAND_POWER_B+1.0))*
      (pow(e2,BAND_POWER_B+1.0)-pow(emin,BAND_POWER_B+1.0));
  }
  return s1+s2;
}

double grb_kcorrect(double z) {
  double e1 = 15.0*(1.0+z);
  double e2 = 150.0*(1.0+z);
  return int_band_fn(e1, e2)/int_band_fn(15, 150);
}

double grb_kcorrect_exact(double z) {
  double e1 = 15.0*(1.0+z);
  double e2 = 150.0*(1.0+z);
  return exact_int_band_fn(e1, e2)/exact_int_band_fn(15, 150);
}

double ld_factor(double z) {
  if (!(z > 1e-9)) z = 1e-9;
  double ld = Dl(z);
  return (ld*ld/(1.0+z));
}

double intrinsic_luminosity_fn(double l) {
  double z = 4.0;
  double ld_norm = ld_factor(0.1)/grb_kcorrect_exact(0.1);
  double lg = ld_factor(z)/grb_kcorrect_exact(z)/ld_norm;
  double dz = 0.05;
  while (fabs(lg-l)>0.001 && dz>0.001) {
    double lg2 = ld_factor(z+dz)/grb_kcorrect_exact(z+dz)/ld_norm;
    double nz = z + dz*(l-lg)/(lg2-lg);
    if (nz < z/2.0) nz = z/2.0;
    if (nz > 2.0*z) nz = 2.0*z;
    dz /= 2.0;
    z = nz;
    lg = ld_factor(z)/grb_kcorrect_exact(z)/ld_norm;
  }
  return instrument_completeness(z, 1);
}

double obs_luminosity_fn(double l) {
  double s = 0;
  int64_t i;
  double ld_norm = ld_factor(0.1)/grb_kcorrect_exact(0.1);
  for (i=0; i<1000; i++) {
    double z = (i+0.5)/120.0;
    double v = Vc((i+1.0)/120.0) - Vc((i)/120.0);
    double lz = l*ld_factor(z)/grb_kcorrect_exact(z)/ld_norm;
    //    double dt = scale_to_years(1.0/((i)/120.0+1.0))-scale_to_years(1.0/((i+1.0)/120.0+1.0));
    double rate = (287.0/373.0)*(82170.0/5037.0)/(pow(10,-1*(z-1.243)) + pow(10,0.241*(z-1.243)));
    s += intrinsic_luminosity_fn(lz)*v*rate/(1.0+z);
  }
  return s*1e-8;
}



float *grb_rates[NUM_PLAWS], *grb_dens_cache[NUM_PLAWS];
float smf[NUM_PLAWS][NUM_VELS][NUM_RADS][NUM_SMS+1];
float smf_all[NUM_PLAWS][NUM_VELS][NUM_RADS][NUM_SMS+1];
float smf_sfr_sm[2][NUM_SMS+1];
float ssfrf[NUM_PLAWS][NUM_VELS][NUM_RADS][NUM_SFRS+1];
float ssfrf_all[NUM_PLAWS][NUM_VELS][NUM_RADS][NUM_SFRS+1];
float ssfrf_sfr_sm[2][NUM_SFRS+1];
float grb_lf[NUM_PLAWS][NUM_GRBS+1];
float grb_escape_f[NUM_PLAWS][NUM_VELF+1];

int main(int argc, char **argv)
{
  float m;
  //float *zs = {0.1, 1, 2, 3, 4, 5, 6, 7, 8};
  struct smf_fit the_smf;
  int64_t i, j, k, lt;
  float fm, ft, t;
  FILE *grb_f, *grb_dens_f, *grb_nd_f, *grb_ef, *grb_h;
  char buffer[1024];
  float plaws[NUM_PLAWS] = {-1.5, -1, -0.5};
  float vels[NUM_VELS] = {180.0, 90.0, 45.0};
  float rads[NUM_RADS] = {5.0, 10.0, 15.0, 7.5, 12.5};
  float missed_hosts[NUM_PLAWS][NUM_L_THRESHES+1];
  float missed_hosts_mean_z[NUM_PLAWS][NUM_L_THRESHES+1];
  float missed_hosts_ks[NUM_PLAWS][NUM_L_THRESHES+1];
  float missed_hosts_ks_mean_z[NUM_PLAWS][NUM_L_THRESHES+1];
  float *lums = NULL;
  float *lums_ks = NULL;

  init_cosmology(0.27, 0.73, 0.7);
  init_time_table(0.27, 0.7);
  for (i=-30; i<30; i++) {
    double l = pow(10,i/5.0);
    printf("%f %e\n", l*1.5, intrinsic_luminosity_fn(l));
  }
  return 0;

  if (argc<5+NUM_PARAMS) {
    fprintf(stderr, "Usage: %s mass_cache weighting (mcmc output)\n", argv[0]);
    exit(1);
  }
  int64_t weighting = atol(argv[2]);
  for (i=0; i<NUM_PARAMS; i++)
    the_smf.params[i] = atof(argv[i+3]);

  r250_init(87L);
  gen_exp10cache();
  setup_psf(1);
  load_mf_cache(argv[1]);
  init_timesteps();
  INVALID(the_smf) = 0;
  calc_sfh(&the_smf);
  load_luminosities("../smf_mcmc2/data/mags_ab_zall_dust.dat");
  gen_all_luminosities(L_f160w, &lums, -1, 0);
  gen_all_luminosities(L_2mass_k, &lums_ks, -1, 0);

  sprintf(buffer, "grb.dat");
  grb_f = check_fopen(buffer, "w");
  sprintf(buffer, "grb_escape_f.dat");
  grb_ef = check_fopen(buffer, "w");
  sprintf(buffer, "grb_hostless.dat");
  grb_h = check_fopen(buffer, "w");
  sprintf(buffer, "grb_dens.dat");
  grb_dens_f = check_fopen(buffer, "w");
  grb_nd_f = check_fopen("grb_nd.dat", "w");
  FILE *grb_nd_weighted_f = check_fopen("grb_nd_weighted.dat", "w");
  FILE *grb_smf_f = check_fopen("grb_smf.dat", "w");
  FILE *grb_smf_all_f = check_fopen("grb_smf_all_grb.dat", "w");
  FILE *grb_smf_sfr_sm_f = check_fopen("grb_smf_sfr_sm.dat", "w");
  FILE *grb_sm_offsets_f = check_fopen("grb_sm_offsets.dat", "w");

  float *total_grb_nd[NUM_PLAWS];
  float *total_grb_nd_ef[NUM_PLAWS][NUM_VELS][NUM_RADS];
  float tgrb[NUM_PLAWS];
  float tgrb_weighted[NUM_PLAWS];
  float tgrb_ef[NUM_PLAWS][NUM_VELS][NUM_RADS];

  FILE *output = check_fopen("grb_legends.txt", "w");
  FILE *output2 = check_fopen("grb_plaw_legends.txt", "w");

  FILE *ev = check_fopen("halo_ev.dat", "w");
  fclose(ev);
  int64_t p,v,r;

  for (p=0; p<NUM_PLAWS; p++) {
    total_grb_nd[p] = check_realloc(NULL, sizeof(float)*num_outputs*3, "");
    memset(total_grb_nd[p], 0, sizeof(float)*num_outputs*3);
    fprintf(output2, "P(\\xD\\f{}t) = (\\xD\\f{}t)\\S%g\\N\n", plaws[p]);
    tgrb[p] = tgrb_weighted[p] = 0;
    memset(grb_lf[p], 0, sizeof(float)*(NUM_GRBS+1));
    memset(grb_escape_f[p], 0, sizeof(float)*(NUM_VELF+1));
    for (v=0; v<NUM_VELS; v++) {
      for (r=0; r<NUM_RADS; r++) {
	total_grb_nd_ef[p][v][r] = check_realloc(NULL, sizeof(float)*num_outputs*3, "");
	memset(total_grb_nd_ef[p][v][r], 0, sizeof(float)*num_outputs*3);
	fprintf(output, "P(\\xD\\f{}t) = (\\xD\\f{}t)\\S%g\\N\\n<v> = %g km/s\\nR\\sthresh\\N = %g R\\se\\N\n", 
		plaws[p], vels[v], rads[r]);
	tgrb_ef[p][v][r] = 0;
	memset(smf[p][v][r], 0, sizeof(float)*(NUM_SMS+1));
	memset(smf_all[p][v][r], 0, sizeof(float)*(NUM_SMS+1));
	memset(ssfrf_all[p][v][r], 0, sizeof(float)*(NUM_SFRS+1));
	memset(ssfrf[p][v][r], 0, sizeof(float)*(NUM_SFRS+1));
      }
    }
  }
  fclose(output);
  fclose(output2);
  memset(smf_sfr_sm[0], 0, sizeof(float)*(NUM_SMS+1));
  memset(smf_sfr_sm[1], 0, sizeof(float)*(NUM_SMS+1));
  memset(ssfrf_sfr_sm[0], 0, sizeof(float)*(NUM_SFRS+1));
  memset(ssfrf_sfr_sm[1], 0, sizeof(float)*(NUM_SFRS+1));


  for (p=0; p<NUM_PLAWS; p++) {
    grb_rates[p] = check_realloc(NULL, sizeof(float)*num_outputs*3*(15.1-8.0)*20, "");
    grb_dens_cache[p] = check_realloc(NULL, sizeof(float)*num_outputs*3*(15.1-8.0)*20, "");
    memset(missed_hosts[p], 0, sizeof(float)*(NUM_L_THRESHES+1));
    memset(missed_hosts_ks[p], 0, sizeof(float)*(NUM_L_THRESHES+1));
    memset(missed_hosts_mean_z[p], 0, sizeof(float)*(NUM_L_THRESHES+1));
    memset(missed_hosts_ks_mean_z[p], 0, sizeof(float)*(NUM_L_THRESHES+1));
  }

  for (m=8; m<15.1; m+=0.05) {
    j = (m-M_MIN-INV_BPDEX*0.5)*BPDEX;
    fm = (m-M_MIN)*BPDEX - j;
    for (t=1; t<(num_outputs-1)*3; t++) {
      ft = (float)(((int)t)%3) / 3.0;
      i = t/3;
      float tscale = steps[i].scale + ft*(steps[i+1].scale - steps[i].scale);
      float a1 = tscale - 0.5*(steps[i+1].scale - steps[i].scale)/3.0;
      float a2 = tscale + 0.5*(steps[i+1].scale - steps[i].scale)/3.0;
      float vol = comoving_volume(1.0/a1 - 1.0) - comoving_volume(1.0/a2 - 1.0);
      float scale = 1.0/tscale;
      float grb_rate[NUM_PLAWS] = {0};
      float ttot = TMIN;
      float sfr = biterp(steps[i].sfr[j], steps[i+1].sfr[j],
                         steps[i].sfr[j+1], steps[i+1].sfr[j+1],
                         ft, fm);
      float sm1 = calc_sm_at_m(m, steps[i].smhm)+steps[i].smhm.mu;
      float sm2 = calc_sm_at_m(m, steps[i+1].smhm)+steps[i+1].smhm.mu;
      float sm = sm1 + ft*(sm2-sm1);
      if (!isfinite(sfr)) sfr = -1000;
      for (k=t-3; k>-1; k--) {
	float sfr = sfh_m_index(j, fm, t, k);
	for (p=0; p<NUM_PLAWS; p++) {
	  float weight = plaw_weight(steps[i].dt/3.0, ttot, plaws[p]);
	  grb_rate[p] += weight*sfr;
	}
	ttot += steps[i].dt/3.0;
      }

      float nd = biterp(steps[i].t[j], steps[i+1].t[j],
                        steps[i].t[j+1], steps[i+1].t[j+1],
                        ft, fm);
      nd += log10(BPDEX);
      
      for (p=0; p<NUM_PLAWS; p++) {
	float grb_dens = pow(10, nd)*grb_rate[p];
	grb_rates[p][(int)((m-8.0)*20*3*num_outputs + t)] = grb_rate[p];
	grb_dens_cache[p][(int)((m-8.0)*20*3*num_outputs + t)] = grb_dens;
	grb_rate[p] = (grb_rate[p] > 0) ? log10(grb_rate[p]) : -5;
	if (!isfinite(grb_dens)) continue;
	//Corrections for volume, instrument completeness, and time delay
	if (weighting) grb_dens *= vol*tscale;
	if (weighting>1) grb_dens *= instrument_completeness(1.0/tscale - 1.0, weighting);
	missed_hosts[p][NUM_L_THRESHES] += grb_dens;
	missed_hosts_ks[p][NUM_L_THRESHES] += grb_dens;
	for (lt=0; lt<NUM_L_THRESHES; lt++) {
	  float thresh = 10.0+lt*0.25;
	  float mh = grb_dens*missed_hosts_fraction(lums, m, 1.0/tscale-1.0, 
						    i, thresh);
	  missed_hosts[p][lt] += mh;
	  missed_hosts_mean_z[p][lt] += mh*(1.0/tscale - 1.0);
	  float mh_ks = grb_dens*missed_hosts_fraction(lums_ks, m, 
						       1.0/tscale-1.0, 
						       i, thresh);
	  missed_hosts_ks[p][lt] += mh_ks;
	  missed_hosts_ks_mean_z[p][lt] += mh_ks*(1.0/tscale - 1.0);
	}
      }

      fprintf(grb_f, "%f %f", scale, m);
      fprintf(grb_dens_f, "%f %f", scale, m);
      fprintf(grb_ef, "%f %f", scale, m);
      for (p=0; p<NUM_PLAWS; p++) {
	float grb_dens = grb_rate[p] + nd;
	float grb_dens_vw = grb_dens;
	if (weighting) grb_dens_vw += log10(vol*tscale);
	if (weighting > 1) grb_dens_vw += log10(instrument_completeness(1.0/tscale - 1.0, weighting));
	float escape_v = escape_velocity(m, 1.0/tscale-1.0, 5);
	int64_t escape_bin = (log10(escape_v)-VELF_MIN)*VELF_BPDEX+0.5;
	float contrib = pow(10, grb_dens_vw);
	if ((escape_bin >= 0) && (escape_bin<NUM_VELF) && isfinite(contrib)) {
	  grb_escape_f[p][escape_bin] += contrib;
	  grb_escape_f[p][NUM_VELF] += contrib;
	}

	if (isfinite(grb_dens)) total_grb_nd[p][(int)t]+=pow(10, grb_dens_vw)*0.05;
	float missed_fraction = missed_hosts_fraction(lums, m, 1.0/tscale-1.0, i, FONG_THRESH);
	for (v=0; v<NUM_VELS; v++) {
	  if (v==1 && p==1) {
	    float all_weight = nd+log10(vol*tscale) + log10(instrument_completeness(1.0/tscale - 1.0, weighting));
	    float sfr_weight = pow(10, all_weight + sfr);
	    if (!isfinite(all_weight)) all_weight = -1000;
	    if (!isfinite(sfr_weight)) sfr_weight = 0;
	    smf_sfr_sm[0][NUM_SMS] += sfr_weight;
	    for (lt=0; lt<NUM_SMS; lt++) {
	      float tsm = SMF_MIN + (double)lt/(double)SMF_BPDEX;
	      float dsm = (sm - tsm)/steps[i].smhm.scatter;
	      smf_sfr_sm[0][lt] += sfr_weight*exp(-0.5*dsm*dsm)*M_2_SQRTPI/(2.0*M_SQRT2*steps[i].smhm.scatter);
	    }
	    for (lt=0; lt<NUM_SMS; lt++) {
	      float tsm = SMF_MIN + (double)lt/(double)SMF_BPDEX;
	      float dsm = (sm - tsm)/steps[i].smhm.scatter;
	      float weight = pow(10, all_weight+tsm)*exp(-0.5*dsm*dsm)*M_2_SQRTPI/(2.0*M_SQRT2*steps[i].smhm.scatter);
	      smf_sfr_sm[1][lt] += SMF_BPDEX*weight;
	      smf_sfr_sm[1][NUM_SMS] += weight;
	    }	    

	    /*int64_t trial;
	    for (trial=0; trial<100; trial++) {
	      float tsm = sm + normal_random(0, steps[i].smhm.scatter);
	      float vel = log(dr250())*vels[v];
	      float re = re_obs(m, 1.0/tscale-1.0, vel);
	      if (re > 0) fprintf(grb_sm_offsets_f, "%f %f %e\n", tsm, re, pow(10, grb_dens_vw));
	      }*/
	  }

	  for (r=0; r<NUM_RADS; r++) {
	    float ef = escape_fraction(vels[v], m, 1.0/tscale-1.0, rads[r]);
	    if (weighting > 2) ef = 1.0 - (1.0-ef)*(1.0-missed_fraction);
	    fprintf(grb_ef, " %f", ef);
	    if (isfinite(grb_dens_vw))
	      total_grb_nd_ef[p][v][r][(int)t]+=pow(10, grb_dens_vw)*ef*0.05;
	    else continue;
	    double gd = pow(10, grb_dens_vw)*(1.0-ef);
	    smf_all[p][v][r][NUM_SMS] += gd;
	    for (lt=0; lt<NUM_SMS; lt++) {
	      float tsm = SMF_MIN + (double)lt/(double)SMF_BPDEX;
	      float dsm = (sm - tsm)/steps[i].smhm.scatter;
	      smf_all[p][v][r][lt] += gd*exp(-0.5*dsm*dsm)*M_2_SQRTPI/(2.0*M_SQRT2*steps[i].smhm.scatter);
	    }
	    float pf = passive_frac(sm, 1.0/tscale-1.0);
	    if (pf > 0.99) pf = 0.99;
	    float scatter = sqrt(steps[i].smhm.scatter*steps[i].smhm.scatter
				 + 0.3*0.3);
	    float ssfr1 = (sfr-log10(1.0-pf))-sm;
	    if (!isfinite(ssfr1)) ssfr1 = -1000;
	    float ssfr2 = ssfr1 - 2.0;
	    ssfrf_all[p][v][r][NUM_SFRS] += gd;
	    for (lt=0; lt<NUM_SFRS; lt++) {
	      float tsfr = SFR_MIN + (double)lt/(double)SFR_BPDEX;
	      float dsm1 = (ssfr1 - tsfr)/scatter;
	      float dsm2 = (ssfr2 - tsfr)/scatter;
	      ssfrf_all[p][v][r][lt] += gd*((1.0-pf)*exp(-0.5*dsm1*dsm1)+pf*exp(-0.5*dsm2*dsm2))*M_2_SQRTPI/(2.0*M_SQRT2*scatter);
	    }

	    if (tscale < (1.0/(1.0+SMF_MAX_Z))) continue;
	    smf[p][v][r][NUM_SMS] += gd;
	    ssfrf[p][v][r][NUM_SFRS] += gd;
	    for (lt=0; lt<NUM_SMS; lt++) {
	      float tsm = SMF_MIN + (double)lt/(double)SMF_BPDEX;
	      float dsm = (sm - tsm)/steps[i].smhm.scatter;
	      smf[p][v][r][lt] += gd*exp(-0.5*dsm*dsm)*M_2_SQRTPI/(2.0*M_SQRT2*steps[i].smhm.scatter);
	    }
	    for (lt=0; lt<NUM_SFRS; lt++) {
	      float tsfr = SFR_MIN + (double)lt/(double)SFR_BPDEX;
	      float dsm1 = (ssfr1 - tsfr)/scatter;
	      float dsm2 = (ssfr2 - tsfr)/scatter;
	      ssfrf[p][v][r][lt] += gd*((1.0-pf)*exp(-0.5*dsm1*dsm1)+pf*exp(-0.5*dsm2*dsm2))*M_2_SQRTPI/(2.0*M_SQRT2*scatter);	      
	    }
	  }
	}
	if (grb_rate[p] < -5) grb_rate[p] = -5;
	if (grb_dens < -8) grb_dens = -8;
	fprintf(grb_f, " %f", grb_rate[p]);
	fprintf(grb_dens_f, " %f", grb_dens);
      }
      fprintf(grb_f, " %f %f\n", sm, sfr);
      fprintf(grb_ef, " %f %f\n", sm, sfr);
      fprintf(grb_dens_f, " %f %f\n", sm, sfr);
    }
  }

  fclose(grb_sm_offsets_f);
  fclose(grb_f);
  fclose(grb_ef);
  fclose(grb_dens_f);

  for (t=1; t<(num_outputs-1)*3; t++) {
    ft = (float)(((int)t)%3) / 3.0;
    i = t/3;
    float tscale = steps[i].scale + ft*(steps[i+1].scale - steps[i].scale);
    //float a1 = tscale + 0.5*(steps[i+1].scale - steps[i].scale)/3.0;
    //float a2 = tscale - 0.5*(steps[i+1].scale - steps[i].scale)/3.0;
    //float dz = fabs(1.0/a1-1.0/a2);
    //float dt = steps[i].dt/3.0;
    fprintf(grb_nd_f, "%f", tscale);
    fprintf(grb_nd_weighted_f, "%f", tscale);
    fprintf(grb_h, "%f", tscale);
    for (p=0; p<NUM_PLAWS; p++) {
      fprintf(grb_nd_f, " %f", total_grb_nd[p][(int)t]);
      tgrb[p] += total_grb_nd[p][(int)t];
      float weighted_total = total_grb_nd[p][(int)t]; //*dt/dz;
      fprintf(grb_nd_weighted_f, " %f", weighted_total);
      tgrb_weighted[p] += weighted_total; //*dz;
      if (!total_grb_nd[p][(int)t]) total_grb_nd[p][(int)t]=1;
      for (v=0; v<NUM_VELS; v++) {
	for (r=0; r<NUM_RADS; r++) {
	  fprintf(grb_h, " %f", total_grb_nd_ef[p][v][r][(int)t]/total_grb_nd[p][(int)t]);
	  tgrb_ef[p][v][r]+=total_grb_nd_ef[p][v][r][(int)t];
	}
      }
    }
    fprintf(grb_nd_f, "\n");
    fprintf(grb_nd_weighted_f, "\n");
    fprintf(grb_h, "\n");
  }
  fclose(grb_nd_f);
  fclose(grb_nd_weighted_f);
  fclose(grb_h);

  double ld_norm = ld_factor(0.1)/grb_kcorrect_exact(0.1);	
  for (t=1; t<(num_outputs-1)*3; t++) {
    ft = (float)(((int)t)%3) / 3.0;
    i = t/3;
    float tscale = steps[i].scale + ft*(steps[i+1].scale - steps[i].scale);
    //float dt = steps[i].dt/3.0;
    float z = 1.0/tscale - 1.0;
    double ldf = ld_factor(z)/grb_kcorrect_exact(z)/ld_norm;
    for (p=0; p<NUM_PLAWS; p++) {
      double sgrbs = total_grb_nd[p][(int)t]/tgrb_weighted[p];
      for (lt=0; lt<NUM_GRBS; lt++) {
	double l = pow(10, GRB_LF_MIN+(double)lt/(double)GRB_BPDEX)*ldf;
	double w = intrinsic_luminosity_fn(l)*sgrbs;
	grb_lf[p][lt]+=w;
	grb_lf[p][NUM_GRBS]+=w;
      }
    }
  }

  FILE *grb_lfs = check_fopen("grb_lfs.dat", "w");
  fprintf(grb_lfs, "#L[ph/s/cm^2] ND\n");
  for (lt=0; lt<NUM_GRBS; lt++) {
    fprintf(grb_lfs, "%e", pow(10, GRB_LF_MIN+(double)lt/(double)GRB_BPDEX)*1.4);
    for (p=0; p<NUM_PLAWS; p++)
      fprintf(grb_lfs, " %e", grb_lf[p][lt]/grb_lf[p][NUM_GRBS]);
    fprintf(grb_lfs, "\n");
  }
  fclose(grb_lfs);

  fprintf(grb_smf_f, "#SM SMFS(z<%f)\n", SMF_MAX_Z);
  for (lt=0; lt<NUM_SMS; lt++) {
    fprintf(grb_smf_f, "%e", pow(10, SMF_MIN+(double)lt/(double)SMF_BPDEX));
    for (p=0; p<NUM_PLAWS; p++)
      for (v=0; v<NUM_VELS; v++)
	for (r=0; r<NUM_RADS; r++)
	  fprintf(grb_smf_f, " %e", smf[p][v][r][lt]/smf[p][v][r][NUM_SMS]);
    fprintf(grb_smf_f, "\n");
  }
  fclose(grb_smf_f);

  FILE *grb_ssfrf_f = check_fopen("grb_ssfrf.dat", "w");
  fprintf(grb_ssfrf_f, "#SM SSFRFS(z<%f)\n", SMF_MAX_Z);
  for (lt=0; lt<NUM_SFRS; lt++) {
    fprintf(grb_ssfrf_f, "%e", pow(10, SFR_MIN+(double)lt/(double)SFR_BPDEX));
    for (p=0; p<NUM_PLAWS; p++)
      for (v=0; v<NUM_VELS; v++)
	for (r=0; r<NUM_RADS; r++)
	  fprintf(grb_ssfrf_f, " %e", ssfrf[p][v][r][lt]/ssfrf[p][v][r][NUM_SFRS]);
    fprintf(grb_ssfrf_f, "\n");
  }
  fclose(grb_ssfrf_f);

 FILE *grb_ssfrf_all_f = check_fopen("grb_ssfrf_all.dat", "w");
 fprintf(grb_ssfrf_all_f, "#SM SSFRFS(all z)\n");
  for (lt=0; lt<NUM_SFRS; lt++) {
    fprintf(grb_ssfrf_all_f, "%e", pow(10, SFR_MIN+(double)lt/(double)SFR_BPDEX));
    for (p=0; p<NUM_PLAWS; p++)
      for (v=0; v<NUM_VELS; v++)
	for (r=0; r<NUM_RADS; r++)
	  fprintf(grb_ssfrf_all_f, " %e", ssfrf_all[p][v][r][lt]/ssfrf_all[p][v][r][NUM_SFRS]);
    fprintf(grb_ssfrf_all_f, "\n");
  }
  fclose(grb_ssfrf_all_f);

  fprintf(grb_smf_all_f, "#SM ALL SMFS\n");
  for (lt=0; lt<NUM_SMS; lt++) {
    fprintf(grb_smf_all_f, "%e", pow(10, SMF_MIN+(double)lt/(double)SMF_BPDEX));
    for (p=0; p<NUM_PLAWS; p++)
      for (v=0; v<NUM_VELS; v++)
	for (r=0; r<NUM_RADS; r++)
	  fprintf(grb_smf_all_f, " %e", smf_all[p][v][r][lt]/smf_all[p][v][r][NUM_SMS]);
    fprintf(grb_smf_all_f, "\n");
  }
  fclose(grb_smf_all_f);


  fprintf(grb_smf_sfr_sm_f, "#SFR-weighted SMF, SM-weighted SMF\n");
  for (lt=0; lt<NUM_SMS; lt++) {
    fprintf(grb_smf_sfr_sm_f, "%e", pow(10, SMF_MIN+(double)lt/(double)SMF_BPDEX));
    for (p=0; p<2; p++)
	  fprintf(grb_smf_sfr_sm_f, " %e", smf_sfr_sm[p][lt]/smf_sfr_sm[p][NUM_SMS]);
    fprintf(grb_smf_sfr_sm_f, "\n");
  }
  fclose(grb_smf_sfr_sm_f);

  sprintf(buffer, "grb_missed_hosts_w%"PRId64".dat", weighting);
  output = check_fopen(buffer, "w");
  fprintf(output, "#Lum.Thresh(F160W); Missed hosts; Average Z of missed hosts\n");
  sprintf(buffer, "grb_missed_hosts_ks_w%"PRId64".dat", weighting);
  output2 = check_fopen(buffer, "w");
  fprintf(output2, "#Lum.Thresh(2MASS_Ks); Missed hosts; Average Z of missed hosts\n");

  for (lt=0; lt<NUM_L_THRESHES; lt++) {
    fprintf(output, "%f", 10+lt*0.25);
    fprintf(output2, "%f", 10+lt*0.25);
    for (p=0; p<NUM_PLAWS; p++) {
      fprintf(output, " %f", missed_hosts[p][lt]/missed_hosts[p][NUM_L_THRESHES]);
      fprintf(output2, " %f", missed_hosts_ks[p][lt]/missed_hosts_ks[p][NUM_L_THRESHES]);
    }
    for (p=0; p<NUM_PLAWS; p++) {
      fprintf(output, " %f", missed_hosts_mean_z[p][lt]/missed_hosts[p][lt]);
      fprintf(output2, " %f", missed_hosts_ks_mean_z[p][lt]/missed_hosts_ks[p][lt]);
    }
    fprintf(output, "\n");
    fprintf(output2, "\n");
  }
  fclose(output);
  fclose(output2);

  FILE *grb_escape_energy_f = check_fopen("grb_escape_energy.dat", "w");
  fprintf(grb_escape_energy_f, "#V CUM_SGRB_fraction\n");
  float cum_grb_escape_energy[3] = {0};
  for (lt=0; lt<NUM_VELF; lt++) {
    float v = pow(10, VELF_MIN+(double)lt/(double)VELF_BPDEX);
    fprintf(grb_escape_energy_f, "%e", v);
    for (p=0; p<NUM_PLAWS; p++) {
      cum_grb_escape_energy[p]+=grb_escape_f[p][lt];
      fprintf(grb_escape_energy_f, " %f", cum_grb_escape_energy[p]/grb_escape_f[p][NUM_VELF]);
    }
    fprintf(grb_escape_energy_f, "\n");
  }
  fclose(grb_escape_energy_f);

  output = check_fopen("grb_norm.dat", "w");
  sprintf(buffer, "grb_fracs_w%"PRId64".dat", weighting);
  output2 = check_fopen(buffer, "w");
  FILE *norm_weighted = check_fopen("grb_norm_weighted.dat", "w");
  for (p=0; p<NUM_PLAWS; p++) {
    fprintf(output, "%e\n", tgrb[p]);
    fprintf(norm_weighted, "%e\n", tgrb_weighted[p]);
    for (v=0; v<NUM_VELS; v++) {
      for (r=0; r<NUM_RADS; r++) {
	fprintf(output2, "P=%g; v=%g; r=%g; EF: %g\n", plaws[p], vels[v], rads[r], tgrb_ef[p][v][r]/tgrb[p]);
      }
    }
  }
  fclose(output);
  fclose(norm_weighted);
  fclose(output2);
  return 0;
}
