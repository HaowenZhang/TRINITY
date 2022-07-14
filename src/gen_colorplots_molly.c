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
  return exp(-4.3*z);
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


#define MZ_MIN 5
#define MZ_MAX 10
#define MZ_BPDEX 20
#define NUM_MZ ((MZ_MAX-MZ_MIN)*MZ_BPDEX + 1)

float metallicity_z[NUM_MZ];
float metallicity_all[NUM_MZ];
float *planet_rate = NULL;

void calc_metallicity(float *dist, float *med, float *up, float *down) {
  int64_t i;
  double sum = 0, sum2 = 0;
  *med = *up = *down = 0;
  for (i=0; i<NUM_MZ; i++) sum += dist[i];
  for (i=0; i<NUM_MZ; i++) {
    double mz = MZ_MIN+((double)i)/(double)MZ_BPDEX;
    sum2 += dist[i];
    if (sum2 > sum*(1.0-0.6827)/2.0 && !(*down)) {
      *down = mz + (sum*(1.0-0.6827)/2.0-sum2)/((double)MZ_BPDEX*dist[i]);
    }
    if (sum2 > sum*(1.0)/2.0 && !(*med)) {
      *med = mz + (sum*(1.0)/2.0-sum2)/((double)MZ_BPDEX*dist[i]);
    }
    if (sum2 > sum*(1.0+0.6827)/2.0 && !(*up)) {
      *up = mz + (sum*(1.0+0.6827)/2.0-sum2)/((double)MZ_BPDEX*dist[i]);
    }
  }
}

void calc_median(float *dist, int64_t n, float *med, float *up, float *down) {
  int64_t i;
  double sum = 0, sum2 = 0;
  *med = *up = *down = 0;
  for (i=0; i<n; i++) sum += dist[i];
  for (i=0; i<n; i++) {
    sum2 += dist[i];
    if (sum2 > sum*(1.0-0.6827)/2.0 && !(*down)) {
      *down = i + (sum*(1.0-0.6827)/2.0-sum2)/(dist[i]);
    }
    if (sum2 > sum*(1.0)/2.0 && !(*med)) {
      *med = i + (sum*(1.0)/2.0-sum2)/((double)dist[i]);
    }
    if (sum2 > sum*(1.0+0.6827)/2.0 && !(*up)) {
      *up = i + (sum*(1.0+0.6827)/2.0-sum2)/((double)dist[i]);
    }
  }
}


int main(int argc, char **argv)
{
  float m;
  struct smf_fit the_smf;
  int64_t i, j, k;
  float fm, ft, t;
  char buffer[1024];

  if (argc<6+NUM_PARAMS) {
    fprintf(stderr, "Usage: %s mass_cache normalization power (mcmc output)\n", argv[0]);
    exit(1);
  }
  double fraction = atof(argv[2]);
  fraction /= 0.35; //Chabrier IMF
  double weighting = atof(argv[3]);
  for (i=0; i<NUM_PARAMS; i++)
    the_smf.params[i] = atof(argv[i+4]);

  r250_init(87L);
  gen_exp10cache();
  setup_psf(1);
  load_mf_cache(argv[1]);
  init_timesteps();
  INVALID(the_smf) = 0;
  calc_sfh(&the_smf);

  FILE *mz_hm_f = check_fopen("mz_hm_z.dat", "w");
  sprintf(buffer, "planet_mz_weighted_%g.dat", weighting);
  FILE *output = check_fopen(buffer, "w");
  sprintf(buffer, "sf_metallicity.dat");
  FILE *mz_f = check_fopen(buffer, "w");

  planet_rate = check_realloc(NULL, sizeof(float)*(num_outputs*3+1), "Planet rate");
  memset(planet_rate, 0, sizeof(float)*(num_outputs+1));
  memset(metallicity_all, 0, sizeof(float)*(NUM_MZ));

  float mz_med, mz_up, mz_down;

  for (t=1; t<(num_outputs-1)*3; t++) {
    memset(metallicity_z, 0, sizeof(float)*(NUM_MZ));
    ft = (float)(((int)t)%3) / 3.0;
    i = t/3;
    float tscale = steps[i].scale + ft*(steps[i+1].scale - steps[i].scale);

    for (m=8; m<15.1; m+=0.05) {
      j = (m-M_MIN-INV_BPDEX*0.5)*BPDEX;
      fm = (m-M_MIN)*BPDEX - j;

      float sfr = biterp(steps[i].sfr[j], steps[i+1].sfr[j],
                         steps[i].sfr[j+1], steps[i+1].sfr[j+1],
                         ft, fm);
      float sm1 = calc_sm_at_m(m, steps[i].smhm)+steps[i].smhm.mu;
      float sm2 = calc_sm_at_m(m, steps[i+1].smhm)+steps[i+1].smhm.mu;
      float sm = sm1 + ft*(sm2-sm1);
      double scatter_corr = exp(pow(steps[i].smhm.scatter*log(10), 2)/2.0);
      float median_sfr = sfr - log10(scatter_corr);
      double alpha = 0.238, beta = 0.273, gamma = 6.226;
      float median_mz =  beta*(sm - alpha*median_sfr) + gamma;
      median_mz = metallicity(sm, 1.0/tscale-1.0)+8.69;
      //median_mz = metallicity_moustakas(sm, 1.0/tscale-1.0)+8.69;
      //median_mz = metallicity_mannucci(sm, median_sfr)+8.69;
      if (median_mz < 5) median_mz = 5;
      float scatter_mz = 0.1;
      float nd = biterp(steps[i].t[j], steps[i+1].t[j],
                        steps[i].t[j+1], steps[i+1].t[j+1],
                        ft, fm);
      nd += log10(BPDEX);
      if (!isfinite(median_sfr)) continue;
      fprintf(mz_hm_f, "%f %f %f %f %f\n", 1.0/tscale, m, median_mz, sm, median_sfr);
      for (k=0; k<NUM_MZ; k++) {
	float mz = MZ_MIN + (double)k/(double)MZ_BPDEX;
	double dmz = (mz-median_mz)/scatter_mz;
	double weight = exp(-0.5*dmz*dmz)/scatter_mz;
	double plaw_weight = pow(10, (mz-8.69)*weighting);
	double rnd = pow(10, nd);
	double planet_weight = plaw_weight*weight*rnd;
	planet_rate[(int)t] += planet_weight*pow(10, sfr);
	planet_rate[0] += planet_weight;
	metallicity_z[k] += weight*rnd*pow(10, sfr);
	metallicity_all[k] += weight*rnd*pow(10, sfr);
      }
    }
    calc_metallicity(metallicity_z, &mz_med, &mz_up, &mz_down);
    fprintf(mz_f, "%f %f %f %f\n", tscale, mz_med, mz_up-mz_med, mz_med - mz_down);
  }
  fclose(mz_hm_f);
  calc_metallicity(metallicity_all, &mz_med, &mz_up, &mz_down);
  fprintf(mz_f, "#Overall Median MZ: %f %f %f\n", mz_med, mz_up-mz_med, mz_med - mz_down);
  printf("Overall Median MZ: %f %f %f\n", mz_med, mz_up-mz_med, mz_med - mz_down);
  fclose(mz_f);

  for (t=1; t<(num_outputs-1)*3; t++) {
    ft = (float)(((int)t)%3) / 3.0;
    i = t/3;
    float tscale = steps[i].scale + ft*(steps[i+1].scale - steps[i].scale);
    fprintf(output, "%f %g\n", tscale, planet_rate[(int)t]/planet_rate[0]);
  }
  float planet_med, planet_up, planet_down;
  calc_median(planet_rate+1, (num_outputs-1)*3, &planet_med, &planet_up, &planet_down);
  fprintf(output, "#Median scale: %f\n", steps[(int)(planet_med/3)].scale);
  fclose(output);
  printf("Median Planet Formation Scale: %f (1sig: %f - %f)\n", steps[(int)(planet_med/3)].scale,
	 steps[(int)(planet_down/3)].scale, steps[(int)(planet_up/3)].scale);
  
  return 0;
}

