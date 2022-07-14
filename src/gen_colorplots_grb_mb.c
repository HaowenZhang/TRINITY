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
#include <string.h>

#define NUM_ZS 9
extern int64_t num_outputs;
extern struct timestep *steps;

#define NUM_PLAWS 3
#define NUM_VELS 3
#define NUM_RADS 5
#define TMIN 50e6
#define NUM_L_THRESHES 40
#define FONG_THRESH 26
#define SMF_MIN 7
#define SMF_MAX 12
#define SMF_BPDEX 24
#define NUM_SMS ((SMF_MAX-SMF_MIN)*SMF_BPDEX+1)
#define SMF_MAX_Z 0.2


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
  if (weighting > 3) return exp(-1.0*z);
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
  float ef = 1.0-maxwell_boltzmann_cdf(vel, ev);
  return ef;
}

float plaw_weight(float dt, float ttot, float index) {
  float weight=1;
  if (index==-1.0) return log(1.0+dt/ttot);
  if (index == -1.5) weight = 1.383012e+02/7.309003e-03;
  else if (index == -0.5) weight = 1.383012e+02/4.348812e+06;
  return (weight*(1.0/(index+1.0))*(pow(ttot+dt, index+1.0)-pow(ttot, index+1.0)));
}

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
  //float vels[NUM_VELS] = {360, 180, 90};
  float vels[NUM_VELS] = {180.0, 90.0, 45.0};
  float rads[NUM_RADS] = {5.0, 10.0, 15.0, 7.5, 12.5};
  float *grb_rates[NUM_PLAWS], *grb_dens_cache[NUM_PLAWS];
  float missed_hosts[NUM_PLAWS][NUM_L_THRESHES+1];
  float missed_hosts_mean_z[NUM_PLAWS][NUM_L_THRESHES+1];
  float *lums = NULL;
  float smf[NUM_PLAWS][NUM_VELS][NUM_RADS][NUM_SMS+1];

  if (argc<5+NUM_PARAMS) {
    fprintf(stderr, "Usage: %s mass_cache weighting (mcmc output)\n", argv[0]);
    exit(1);
  }
  int64_t weighting = atol(argv[2]);
  for (i=0; i<NUM_PARAMS; i++)
    the_smf.params[i] = atof(argv[i+3]);

  gen_exp10cache();
  setup_psf(1);
  load_mf_cache(argv[1]);
  init_timesteps();
  INVALID(the_smf) = 0;
  calc_sfh(&the_smf);
  load_luminosities("../smf_mcmc2/data/mags_ab_zall_dust.dat");
  gen_all_luminosities(L_f160w, &lums, -1, 0);

  sprintf(buffer, "grb_mb.dat");
  grb_f = check_fopen(buffer, "w");
  sprintf(buffer, "grb_escape_f_mb.dat");
  grb_ef = check_fopen(buffer, "w");
  sprintf(buffer, "grb_hostless_mb.dat");
  grb_h = check_fopen(buffer, "w");
  sprintf(buffer, "grb_dens_mb.dat");
  grb_dens_f = check_fopen(buffer, "w");
  grb_nd_f = check_fopen("grb_nd_mb.dat", "w");
  FILE *grb_smf_f = check_fopen("grb_smf_mb.dat", "w");

  float *total_grb_nd[NUM_PLAWS];
  float *total_grb_nd_ef[NUM_PLAWS][NUM_VELS][NUM_RADS];
  float tgrb[NUM_PLAWS];
  float tgrb_ef[NUM_PLAWS][NUM_VELS][NUM_RADS];

  FILE *output = check_fopen("grb_legends_mb.txt", "w");
  FILE *output2 = check_fopen("grb_plaw_legends_mb.txt", "w");

  FILE *ev = check_fopen("halo_ev_mb.dat", "w");
  //
  fclose(ev);
  int64_t p,v,r;

  for (p=0; p<NUM_PLAWS; p++) {
    total_grb_nd[p] = check_realloc(NULL, sizeof(float)*num_outputs*3, "");
    memset(total_grb_nd[p], 0, sizeof(float)*num_outputs*3);
    fprintf(output2, "P(\\xt\\f{}) = \\xt\\f{}\\S%g\\N\n", plaws[p]);
    tgrb[p] = 0;
    for (v=0; v<NUM_VELS; v++) {
      for (r=0; r<NUM_RADS; r++) {
	total_grb_nd_ef[p][v][r] = check_realloc(NULL, sizeof(float)*num_outputs*3, "");
	memset(total_grb_nd_ef[p][v][r], 0, sizeof(float)*num_outputs*3);
	fprintf(output, "P(\\xt\\f{}) = \\xt\\f{}\\S%g\\N\\n<v> = %g km/s\\nR\\sthresh\\N = %g R\\se\\N\n", 
		plaws[p], vels[v], rads[r]);
	tgrb_ef[p][v][r] = 0;
	memset(smf[p][v][r], 0, sizeof(float)*(NUM_SMS+1));
      }
    }
  }
  fclose(output);
  fclose(output2);

  for (p=0; p<NUM_PLAWS; p++) {
    grb_rates[p] = check_realloc(NULL, sizeof(float)*num_outputs*3*(15.1-8.0)*20, "");
    grb_dens_cache[p] = check_realloc(NULL, sizeof(float)*num_outputs*3*(15.1-8.0)*20, "");
    memset(missed_hosts[p], 0, sizeof(float)*(NUM_L_THRESHES+1));
    memset(missed_hosts_mean_z[p], 0, sizeof(float)*(NUM_L_THRESHES+1));
  }

  for (m=8; m<15.1; m+=0.05) {
    j = (m-M_MIN-INV_BPDEX*0.5)*BPDEX;
    fm = (m-M_MIN)*BPDEX - j;
    for (t=1; t<(num_outputs-1)*3; t++) {
      ft = (float)(((int)t)%3) / 3.0;
      i = t/3;
      float tscale = steps[i].scale + ft*(steps[i+1].scale - steps[i].scale);
      float a1 = tscale - 0.5*ft*(steps[i+1].scale - steps[i].scale);
      float a2 = tscale + 0.5*ft*(steps[i+1].scale - steps[i].scale);
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
	for (lt=0; lt<NUM_L_THRESHES; lt++) {
	  float thresh = 20.0+lt*0.25;
	  float mh = grb_dens*missed_hosts_fraction(lums, m, 1.0/tscale-1.0, 
						    i, thresh);
	  missed_hosts[p][lt] += mh;
	  missed_hosts_mean_z[p][lt] += mh*(1.0/tscale - 1.0);
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
	if (isfinite(grb_dens)) total_grb_nd[p][(int)t]+=pow(10, grb_dens_vw)*0.05;
	float missed_fraction = missed_hosts_fraction(lums, m, 1.0/tscale-1.0, i, FONG_THRESH);
	for (v=0; v<NUM_VELS; v++) {
	  for (r=0; r<NUM_RADS; r++) {
	    float ef = escape_fraction(vels[v], m, 1.0/tscale-1.0, rads[r]);
	    if (weighting > 2) ef = 1.0 - (1.0-ef)*(1.0-missed_fraction);
	    fprintf(grb_ef, " %f", ef);
	    if (isfinite(grb_dens_vw))
	      total_grb_nd_ef[p][v][r][(int)t]+=pow(10, grb_dens_vw)*ef*0.05;
	    else continue;
	    if (tscale < (1.0/(1.0+SMF_MAX_Z))) continue;
	    double gd = pow(10, grb_dens_vw)*(1.0-ef);
	    smf[p][v][r][NUM_SMS] += gd;
	    for (lt=0; lt<NUM_SMS; lt++) {
	      float tsm = SMF_MIN + (double)lt/(double)SMF_BPDEX;
	      float dsm = (sm - tsm)/steps[i].smhm.scatter;
	      smf[p][v][r][lt] += gd*exp(-0.5*dsm*dsm)*M_2_SQRTPI/(2.0*M_SQRT2*steps[i].smhm.scatter);
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

  fclose(grb_f);
  fclose(grb_ef);
  fclose(grb_dens_f);

  for (t=1; t<(num_outputs-1)*3; t++) {
    ft = (float)(((int)t)%3) / 3.0;
    i = t/3;
    float tscale = steps[i].scale + ft*(steps[i+1].scale - steps[i].scale);
    fprintf(grb_nd_f, "%f", tscale);
    fprintf(grb_h, "%f", tscale);
    for (p=0; p<NUM_PLAWS; p++) {
      fprintf(grb_nd_f, " %f", total_grb_nd[p][(int)t]);
      tgrb[p] += total_grb_nd[p][(int)t];
      if (!total_grb_nd[p][(int)t]) total_grb_nd[p][(int)t]=1;
      for (v=0; v<NUM_VELS; v++) {
	for (r=0; r<NUM_RADS; r++) {
	  fprintf(grb_h, " %f", total_grb_nd_ef[p][v][r][(int)t]/total_grb_nd[p][(int)t]);
	  tgrb_ef[p][v][r]+=total_grb_nd_ef[p][v][r][(int)t];
	}
      }
    }
    fprintf(grb_nd_f, "\n");
    fprintf(grb_h, "\n");
  }
  fclose(grb_nd_f);
  fclose(grb_h);

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

  sprintf(buffer, "grb_missed_hosts_w%"PRId64"_mb.dat", weighting);
  output = check_fopen(buffer, "w");
  fprintf(output, "#Lum.Thresh(F160W); Missed hosts; Average Z of missed hosts\n");
  for (lt=0; lt<NUM_L_THRESHES; lt++) {
    fprintf(output, "%f", 20+lt*0.25);
    for (p=0; p<NUM_PLAWS; p++)
      fprintf(output, " %f", missed_hosts[p][lt]/missed_hosts[p][NUM_L_THRESHES]);
    for (p=0; p<NUM_PLAWS; p++)
      fprintf(output, " %f", missed_hosts_mean_z[p][lt]/missed_hosts[p][lt]);
    fprintf(output, "\n");
  }
  fclose(output);

  output = check_fopen("grb_norm_mb.dat", "w");
  sprintf(buffer, "grb_fracs_w%"PRId64"_mb.dat", weighting);
  output2 = check_fopen(buffer, "w");
  for (p=0; p<NUM_PLAWS; p++) {
    fprintf(output, "%e\n", tgrb[p]);
    for (v=0; v<NUM_VELS; v++) {
      for (r=0; r<NUM_RADS; r++) {
	fprintf(output2, "P=%g; v=%g; r=%g; EF: %g\n", plaws[p], vels[v], rads[r], tgrb_ef[p][v][r]/tgrb[p]);
      }
    }
  }
  fclose(output);
  fclose(output2);
  return 0;
}
