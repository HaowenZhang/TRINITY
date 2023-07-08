#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <inttypes.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <assert.h>
#include <unistd.h>
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
#include "universe_time.h"
#include "smloss.h"
#include "mt_rand.h"
#include "make_sf_catalog.h"

//TODO:
//Rewrite to not copy main progenitor SFHs


struct catalog_halo *old=NULL, *new=NULL;
float *sfh_old = NULL;
float *sfh_new = NULL;
float *icl_old = NULL;
float *icl_new = NULL;
int64_t *id_to_idx = NULL;
int64_t min_id=0, max_id = 0;
int64_t num_old = 0, num_new = 0;
float *scales = NULL;
float *dt = NULL;
float *rem = NULL;
float *tmp_sfh = NULL, *tmp_icl = NULL;
int64_t counts[M_BINS];
int64_t bin_start[M_BINS+1];
double q_frac[M_BINS];
double av_vpeak[M_BINS];
double av_vpeak2[M_BINS];



void calc_geo_av_vpeak(void) {
  int64_t i;
  memset(av_vpeak, 0, sizeof(double)*M_BINS);
  memset(av_vpeak2, 0, sizeof(double)*M_BINS);
  for (i=0; i<num_new; i++) {
    int64_t b = new[i].mf;
    double d = log10(new[i].vp);
    av_vpeak[b] += d;
    av_vpeak2[b] += d*d;
  }
  for (i=0; i<M_BINS; i++) {
    if (!counts[i]) continue;
    av_vpeak[i] /= counts[i];
    av_vpeak2[i] /= counts[i];
    av_vpeak2[i] -= av_vpeak[i]*av_vpeak[i];
    if (av_vpeak2[i]>0) av_vpeak2[i] = sqrt(av_vpeak2[i]);
    else av_vpeak2[i] = 1.0;
  }
}

void calc_ages(float scale) {
  int64_t i;
  //double l04 = log(0.04);
  for (i=0; i<num_new; i++) {
    //new[i].age = scale*new[i].a_4p + (1.0-scale)*new[i].a_half;
    //new[i].age = new[i].a_half;
    //new[i].age = new[i].a_4p;
    //new[i].age = new[i].acc_scale * 4.1 / new[i].c;
    //new[i].age = new[i].acc_q / new[i].mp;
    new[i].age = 1.0*(1.0+(new[i].v - new[i].vp)/new[i].vp);

    //new[i].age = new[i].acc_scale * 4.1 / new[i].c;
    //new[i].age = new[i].peak_scale * 4.1 / new[i].pc;
    //new[i].age = new[i].spin_acc;
    //new[i].age = new[i].a_v05;
    //float frac = 2;
    //if (frac > 1) frac = 1;
    /*float obs_scale = new[i].acc_scale;
      new[i].age = obs_scale * 4.1 / new[i].c;*/

    //if (new[i].a_12 < new[i].age) new[i].age = new[i].a_12;

  }
}

float quenched_sm0(float scale) {
  //return pow(10, 10.2+0.5*(1.0/scale-1.0)); //Default
  return pow(10, 10.4+0.5*(1.0/scale-1.0)); //H&W
}

//Modified!!!
float quenched_fraction(float sm, float sm0) {
  if (!(sm>0)) return 0;
  return 1.0/(pow(sm/sm0, -0.5)+1.0);
}

void calc_quenched_fraction(struct smf_fit *t, float scale) {
  float mf[M_BINS], smm[M_BINS], shortfall[M_BINS], shortfall_norm[M_BINS];
  float z = 1.0/scale - 1.0;
  struct smf c = smhm_at_z(z, *t);
  float sm0 = quenched_sm0(scale);
  int i;
  for (i=0; i<M_BINS; i++) {
    float m = M_MIN + (i+0.5)*INV_BPDEX;
    mf[i] = pow(10, mf_cache(scale, m));
    smm[i] =calc_sm_at_m(m, c);
    q_frac[i] = quenched_fraction(pow(10, smm[i]), sm0);
  }

  float sm = 0;
  float scatter = sqrt(c.scatter*c.scatter + c.obs_scatter+c.obs_scatter);
  int j=0, max_j = 8;

  for (j=0; j<max_j+1; j++) {    
    for (i=0; i<M_BINS; i++) shortfall[i] = shortfall_norm[i] = 0;

    for (sm = 8; sm <12; sm+=0.05) {
      double tot=0, totq=0;
      for (i=0; i<M_BINS; i++) {
	float dm = (smm[i]-sm)/scatter;
	float nd = mf[i]*exp(-0.5*dm*dm);
	tot += nd;
	totq += nd*q_frac[i];
      }
      float model_qf = totq/tot;
      float actual_qf = quenched_fraction(pow(10, sm), sm0);
      //if (j==max_j)
      //printf("%f %f %f\n", sm, totq/tot, quenched_fraction(pow(10, sm), sm0));
	
      for (i=0; i<M_BINS; i++) {
	float dm = (smm[i]-sm)/scatter;
	float nd = exp(-0.5*dm*dm);
	shortfall[i] += nd*(actual_qf-model_qf);
	shortfall_norm[i] += nd;
      }
    }
      
    for (i=0; i<M_BINS; i++) {
      q_frac[i] += shortfall[i]/shortfall_norm[i];
      if (q_frac[i] > 1) q_frac[i] = 1;
      if (q_frac[i] < 0) q_frac[i] = 0;
    }
  }
  /*for (i=0; i<M_BINS*10; i++) {
    float m = M_MIN + i*INV_BPDEX/10.0;
    printf("%f %f\n", m, fq[i]);
    }*/
}


int sort_by_age(const void *a, const void *b) {
  const float *c = a;
  const float *d = b;
  if (*c < *d) return -1;
  if (*c > *d) return 1;
  return 0;
}			      

int main(int argc, char **argv) {
  int64_t i;
  FILE *in;
  char buffer[1024];
  struct smf_fit the_smf;
  if (argc<4) {
    printf("Usage: %s mass_cache.dat mcmc_stats scale_file\n", argv[0]);
    exit(1);
  }

  r250_init(87L);
  in = check_fopen(argv[2], "r");
  fgets(buffer, 1024, in);
  read_params(buffer, the_smf.params, NUM_PARAMS);
  fclose(in);
  
  init_time_table(0.27, H0);
  gen_exp10cache();
  setup_psf(1);
  load_mf_cache(argv[1]);
  init_timesteps();
  INVALID(the_smf) = 0;
  calc_sfh(&the_smf);
  check_realloc_s(tmp_sfh, num_outputs, sizeof(float));
  check_realloc_s(tmp_icl, num_outputs, sizeof(float));

  in = check_fopen(argv[3], "r");
  i = 0;
  while (fgets(buffer, 1024, in)) {
    char filename[1024];
    double scale = atof(buffer);
    snprintf(filename, 1024, "acc_list_%f.bin", scale);
    load_halos(filename, i+1);
    calc_ages(scale);
    calc_losses(atof(buffer), i+1);
    calc_sm_histories(i);
    calc_quenched_fraction(&the_smf, scale);
    quench_galaxies(i);
    calc_geo_av_vpeak();
    calc_sf_from_accretion_rates(i);
    write_catalog(i);
    i++;
    if (argc>4) {
      if (!(i%20)) {
	int64_t pid = fork();
	if (pid>0) exit(0);
      }
    }
  }
  fclose(in);
  return 0;
}


#define V_BPDEX 8
#define V_INV_BPDEX (1.0/((double)V_BPDEX))
#define V_MIN 0
#define V_MAX 4
#define V_BINS ((V_MAX-V_MIN)*V_BPDEX+1)

void calc_sf_from_accretion_rates(int64_t n) {
  int64_t i,j;
  double avg_sfr[M_BINS] = {0};
  double icl_frac[M_BINS] = {0};
  double frac_icl_to_transfer[M_BINS] = {0};

  double v_avg_sfr[V_BINS] = {0};
  int64_t v_counts[V_BINS] = {0};

  for (i=0; i<M_BINS; i++) {
    int64_t step;
    double f;
    calc_step_at_z(1.0/scales[n]-1.0, &step, &f);
    if (step > num_outputs-2) { step = num_outputs-2; f = 1; }
    avg_sfr[i] = steps[step].sfr[i]*(1.0-f) + f*steps[step+1].sfr[i];
    assert(!counts[i] || isfinite(avg_sfr[i]));

    //avg_sfr[i] *= 1.0-qf;
    //avg_sfe[i] = avg_sfr[i]/avg_acc[i];

    icl_frac[i] = steps[step].icl_frac[i]*(1.0-f) + f*steps[step+1].icl_frac[i];
    if (icl_frac[i] > 0.99) icl_frac[i] = 0.99;
    double total_icl = 0;
    for (j=bin_start[i]; j<bin_start[i+1]; j++) total_icl += new[j].icl;
    double avg_icl_transfer = avg_sfr[i]*rem[n]*dt[n]*icl_frac[i]/(1.0-icl_frac[i]);
    double max_icl_transfer = total_icl;
    if (counts[i]) max_icl_transfer /= ((double)counts[i]);
    if (max_icl_transfer < avg_icl_transfer) avg_icl_transfer = max_icl_transfer;
    if (!max_icl_transfer) max_icl_transfer = 1.0;
    frac_icl_to_transfer[i] = avg_icl_transfer / max_icl_transfer;
  }

  for (j=0; j<num_new; j++) {
    if (new[j].vmp < 1) continue;
    if (!(avg_sfr[(int64_t)new[j].mf]>0)) continue;
    int64_t vbin = (log10(new[j].vmp)-V_MIN)*V_BPDEX;
    assert(vbin < V_BINS && vbin>=0);
    v_avg_sfr[vbin]+=avg_sfr[(int64_t)new[j].mf];
    v_counts[vbin]++;
  }

  for (i=0; i<V_BINS; i++) {
    if (!v_counts[i]) continue;
    v_avg_sfr[i] /= (double)v_counts[i];
    //if (n==num_outputs-1) printf("%f %f\n", V_MIN+i*V_INV_BPDEX, v_avg_sfr[i]);
  }

#define SET(x,y) float x = y[mb]; if (counts[mb+1] && (avg_sfr[mb+1]>0)) x+=f*(y[mb+1]-y[mb]);
#define V_SET(x,y) float x = y[vb]; if (v_counts[vb+1] && (v_avg_sfr[vb+1]>0)) x+=vf*(y[vb+1]-y[vb]);
  double scatter_corr = exp(pow(INTR_SCATTER*log(10), 2)/2.0);
  for (j=0; j<num_new; j++) {
    float f = new[j].mf - 0.5;
    int64_t mb = f;
    f -= mb;
    if (!counts[mb]) { mb = new[j].mf; f = 0; }
    //SET(asfr, avg_sfr);

    double vmp = new[j].vmp;
    if (vmp < 1) vmp = 1;
    float vf = (log10(vmp)-V_MIN)*V_BPDEX - 0.5;
    int64_t vb = vf;
    vf -= vb;
    if (!counts[vb]) { vb++; vf=0; }
    V_SET(asfr, v_avg_sfr);

    //SET(avpeak, av_vpeak);
    //SET(svpeak, av_vpeak2);


    //if (!(svpeak > 0)) svpeak = 1.0;
    float obs_scatter = (new[j].q) ? OBS_SIGMA_Q : OBS_SIGMA_SF;
    //float dv = (log10(new[j].vp) - avpeak)/svpeak;

    if (new[j].q)
      new[j].sfr = pow(10, normal_random(Q_SSFR, INTR_SCATTER))*new[j].sm;
    else
      new[j].sfr = asfr*pow(10, normal_random(0, INTR_SCATTER))/scatter_corr;
    //new[j].sfr = asfr*pow(10, dv*INTR_SCATTER)/scatter_corr;

    new[j].obs_sfr = new[j].sfr*pow(10, normal_random(0, obs_scatter));
    sfh_new[j*(n+1)+n] = new[j].sfr*dt[n];
    new[j].sm += sfh_new[j*(n+1)+n]*rem[n];

    SET(ftt, frac_icl_to_transfer);
    new[j].sm += ftt * new[j].icl;
    new[j].icl *= (1.0 - ftt);
    float *sfh = sfh_new+j*(n+1);
    float *icl = icl_new+j*(n+1);
    int64_t k;
    for (k=0; k<n+1; k++) {
      sfh[k] += ftt*icl[k];
      icl[k] *= (1.0 - ftt);
    }
  }
#undef SET
}

/*
void calc_sf_from_accretion_rates(int64_t n) {
  int64_t i,j;
  double avg_sfr[M_BINS] = {0};
  double avg_acc[M_BINS] = {0};
  double avg_sfe[M_BINS] = {0};
  double q_factor[M_BINS] = {0};
  double icl_frac[M_BINS] = {0};
  double frac_icl_to_transfer[M_BINS] = {0};

  for (i=0; i<M_BINS; i++) {
    int64_t step;
    double f;
    calc_step_at_z(1.0/scales[n]-1.0, &step, &f);
    if (step > num_outputs-2) { step = num_outputs-2; f = 1; }
    avg_sfr[i] = steps[step].sfr[i]*(1.0-f) + f*steps[step+1].sfr[i];
    assert(!counts[i] || isfinite(avg_sfr[i]));
    //avg_sm[i] = steps[step].sm_avg[i]*(1.0-f) + f*steps[step+1].sm_avg[i];
    q_factor[i] = 1.0/200.0;
    //if (avg_sfr[i]>0) { q_factor[i] = 0.3e-11 * (avg_sm[i]/avg_sfr[i]); }

    int64_t tq = 0;
    float qf = 0;
    for (j=bin_start[i]; j<bin_start[i+1]; j++) {
      float acc_q = new[j].acc_q;
      if (acc_q < 0) acc_q = 0;
      if (new[j].q) acc_q *= q_factor[i];
      avg_acc[i] += acc_q;
      if (!(acc_q>0) || new[j].q) tq++;
    }
    if (counts[i]) {
      avg_acc[i] /= (double)counts[i];
      qf = (double)tq / (double)counts[i];
    }
    if (!(avg_acc[i]>0)) { avg_sfr[i] = 0; avg_acc[i] = 1.0; }

    //MAJOR HACK!!!
    avg_sfr[i] *= 1.0-qf;

    avg_sfe[i] = avg_sfr[i]/avg_acc[i];

    icl_frac[i] = steps[step].icl_frac[i]*(1.0-f) + f*steps[step+1].icl_frac[i];
    if (icl_frac[i] > 0.99) icl_frac[i] = 0.99;
    double total_icl = 0;
    for (j=bin_start[i]; j<bin_start[i+1]; j++) total_icl += new[j].icl;
    double avg_icl_transfer = avg_sfr[i]*rem[n]*dt[n]*icl_frac[i]/(1.0-icl_frac[i]);
    double max_icl_transfer = total_icl;
    if (counts[i]) max_icl_transfer /= ((double)counts[i]);
    if (max_icl_transfer < avg_icl_transfer) avg_icl_transfer = max_icl_transfer;
    if (!max_icl_transfer) max_icl_transfer = 1.0;
    frac_icl_to_transfer[i] = avg_icl_transfer / max_icl_transfer;
  }

#define SET(x,y) float x = y[mb]; if (counts[mb+1] && (avg_sfe[mb+1]>0)) x+=f*(y[mb+1]-y[mb]);
  for (j=0; j<num_new; j++) {
    int64_t mb = new[j].mf-0.5;
    float f = (new[j].mf-0.5)-mb;
    if (!counts[mb]) { mb = new[j].mf; f = 0; }
    SET(a_acc,avg_acc);
    SET(qf, q_factor);
    SET(asfr, avg_sfr);
    SET(asfe, avg_sfe);

    float acc_q = new[j].acc_q;
    if (acc_q < 0) acc_q = 0;
    if (new[j].q) { acc_q *= qf; new[j].qf = qf; }
    else { new[j].qf = 1.0; }
    new[j].sfr = asfe * acc_q;
    sfh_new[j*(n+1)+n] = new[j].sfr*dt[n];
    new[j].sm += sfh_new[j*(n+1)+n]*rem[n];

    SET(ftt, frac_icl_to_transfer);
    new[j].sm += ftt * new[j].icl;
    new[j].icl *= (1.0 - ftt);
    float *sfh = sfh_new+j*(n+1);
    float *icl = icl_new+j*(n+1);
    int64_t k;
    for (k=0; k<n+1; k++) {
      sfh[k] += ftt*icl[k];
      icl[k] *= (1.0 - ftt);
    }
  }
#undef SET
}
*/


#define SMIN 0
#define SMAX 13
#define SBPDEX 5
#define SBINS ((int64_t)((SMAX-SMIN)*SBPDEX+2))
#define SINV_BPDEX (1.0/(double)SBPDEX)

void quench_galaxies_sm(int64_t n) {
  int64_t i;
  double quenching_age[SBINS] = {0};
  float *ages[SBINS]={NULL};
  int64_t sm_counts[SBINS] = {0};
  
  float sm0 = quenched_sm0(scales[n]);

  for (i=0; i<num_new; i++) {
    if (!(new[i].sm>0)) continue;
    int64_t smb = (log10(new[i].sm)-SMIN)*SBPDEX;
    if (smb < 0 || smb >= SBINS) continue;
    check_realloc_every(ages[smb], sizeof(float), sm_counts[smb], 1000);
    ages[smb][sm_counts[smb]] = new[i].age;
    sm_counts[smb]++;
  }
  for (i=0; i<SBINS; i++) {
    if (!sm_counts[i]) {
      if (i>0 && quenching_age[i-1])
	quenching_age[i] = quenching_age[i-1];
      continue;
    }
    float qfrac = quenched_fraction(pow(10, SMIN + (i+0.5)*SINV_BPDEX), sm0);
    qsort(ages[i], sm_counts[i], sizeof(float), sort_by_age);
    if (sm_counts[i]) sm_counts[i]--;
    quenching_age[i] = ages[i][(int64_t)(qfrac*sm_counts[i])];
    free(ages[i]);
  }

  for (i=0; i<num_new; i++) {
    new[i].q = 0;
    if (!(new[i].sm>0)) continue;
    float f = (log10(new[i].sm)-SMIN)*SBPDEX - 0.5;
    int64_t mb = f;
    f -= mb;
    if (mb < 0) continue;
    if (mb >= SBINS) { new[i].q=1; continue; }
    if (!sm_counts[mb]) { mb++; f=0; }
    float q_age = quenching_age[mb] + f*(quenching_age[mb+1]-quenching_age[mb]);
    new[i].q = (q_age > new[i].age) ? 1 : 0;
  }
}


void quench_galaxies(int64_t n) {
  int64_t i;
  double quenching_age[M_BINS] = {0};
  float *ages=NULL; 

  for (i=0; i<M_BINS; i++) {
    if (!counts[i]) {
      if (i>0 && quenching_age[i-1])
	quenching_age[i] = quenching_age[i-1];
      continue;
    }
    ages = check_realloc(ages, sizeof(float)*counts[i], "Ages");
    int64_t j;
    for (j=0; j<counts[i]; j++) ages[j] = new[bin_start[i]+j].age;
    qsort(ages, counts[i], sizeof(float), sort_by_age);
    quenching_age[i] = ages[(int64_t)(q_frac[i]*counts[i])];
  }
  free(ages);

  for (i=0; i<num_new; i++) {
    int64_t mb = new[i].mf-0.5;
    float f = (new[i].mf-0.5)-mb;
    float q_age = quenching_age[mb] + f*(quenching_age[mb+1]-quenching_age[mb]);
    new[i].q = (q_age > new[i].age) ? 1 : 0;
  }
}

/*void quench_galaxies(int64_t n) {
  int64_t i;
  float scale = scales[n];
  float sm0 = pow(10, 10.2+0.5*(1.0/scale-1.0));
  int64_t scounts[SBINS] = {0};
  int64_t quenched[SBINS] = {0};
  double exp_quenched[SBINS] = {0};

  for (i=0; i<num_new; i++) {
    if (new[i].sm>0) new[i].smf = (log10(new[i].sm)-SMIN)*SBPDEX;
    else { new[i].smf = 0; }
    int64_t mb = new[i].smf;
    scounts[mb]++;
    if (new[i].q || new[i].acc_q <= 0) quenched[mb]++;
  }
  for (i=0; i<SBINS; i++) {
    exp_quenched[i] = quenched_fraction(pow(10, SMIN+(i+0.5)*SINV_BPDEX), sm0)
      *scounts[i];
    exp_quenched[i] -= quenched[i];
    if (exp_quenched[i]<0) exp_quenched[i] = 0;
    if (quenched[i] < scounts[i]) exp_quenched[i] /= (scounts[i]-quenched[i]);
  }

  for (i=0; i<num_new; i++) {
    if (new[i].q || new[i].acc_q <= 0) continue;
    if (new[i].smf < 0 || new[i].smf >= SBINS-1) continue;
    int64_t mb = new[i].smf - 0.5;
    float f = (new[i].smf - 0.5) - mb;
    if (!scounts[mb]) { mb=new[i].smf; f = 0; }
    float qp = exp_quenched[mb];
    if (scounts[mb+1]) qp += f*(exp_quenched[mb+1]-exp_quenched[mb]);
    if (drand48() < qp) new[i].q = 1;
  }
}
*/

void generate_sm_hist(int64_t n, int64_t num_old_inputs) {
  float z = 1.0/scales[num_old_inputs] - 1.0; //Current scale
  int64_t step;
  double f;
  calc_step_at_z(z, &step, &f);
  float mf = new[n].mf - 0.5;
  int64_t mb = mf;
  mf -= mb;
  assert(mb>=0 && mb<M_BINS-1);
  if (step >= num_outputs-1) { step = num_outputs-2; f=1; }
  int64_t i;
#define INTERP(x,s) (steps[s].x[mb*num_outputs+i]*(1.0-mf) +	\
		       mf*steps[s].x[(mb+1)*num_outputs+i])
  memset(tmp_sfh, 0, sizeof(float)*num_outputs);
  memset(tmp_icl, 0, sizeof(float)*num_outputs);
  for (i=0; i<step+1; i++) {
    float sm2 = INTERP(sm_hist, step+1);
    float icl2 = INTERP(icl_stars, step+1);
    float sm1 = INTERP(sm_hist, step);
    float icl1 = INTERP(icl_stars, step);
    tmp_sfh[i] = sm1 + f*(sm2-sm1);
    if (!isfinite(tmp_sfh[i])) tmp_sfh[i] = 0;
    tmp_icl[i] = icl1 + f*(icl2-icl1);
    if (!isfinite(tmp_icl[i])) tmp_icl[i] = 0;
  }
#undef INTERP
  float tused = 0;
  int64_t j=0;
  i=0;
  float *sfh = sfh_new+(n*(num_old_inputs+1));
  float *icl = icl_new+(n*(num_old_inputs+1));
  while (j<num_old_inputs) {
    float dts = dt[j];
    while (dts > steps[i].dt-tused) {
      sfh[j] += tmp_sfh[i]*(steps[i].dt-tused)/(steps[i].dt);
      icl[j] += tmp_icl[i]*(steps[i].dt-tused)/(steps[i].dt);
      dts -= steps[i].dt-tused;
      if (i<num_outputs-1) i++;
      tused = 0;
    }
    tused = dts;
    sfh[j] += tmp_sfh[i]*tused/(steps[i].dt);
    icl[j] += tmp_icl[i]*tused/(steps[i].dt);
    j++;
  }
}

void calc_sm_histories(int num_old_inputs) {
  int64_t i, j;
  for (i=0; i<num_old; i++) {
    int64_t desc = id_to_index(old[i].descid);
    if (desc < 0) continue;
    float *desc_sfh = sfh_new+(desc*(num_old_inputs+1));
    float *prog_sfh = sfh_old+(i*num_old_inputs);
    float *desc_icl = icl_new+(desc*(num_old_inputs+1));
    float *prog_icl = icl_old+(i*num_old_inputs);
    if (old[i].mp < new[desc].prog_mp) //To ICL
      for (j=0; j<num_old_inputs; j++) desc_icl[j]+=prog_sfh[j];
    else {
      for (j=0; j<num_old_inputs; j++) {
	desc_icl[j]+=prog_icl[j];
	desc_sfh[j]+=prog_sfh[j];
      }
      new[desc].q = old[i].q;
    }
  }

  for (i=0; i<num_new; i++) {
    if (!new[i].prog_mp) generate_sm_hist(i, num_old_inputs);
    float *sfh = sfh_new+(i*(num_old_inputs+1));
    float *icl = icl_new+(i*(num_old_inputs+1));
    for (j=0; j<num_old_inputs; j++) {
      new[i].sm += sfh[j]*rem[j];
      new[i].icl += icl[j]*rem[j];
    }
  }
}

void calc_losses(float scale, int64_t num_inputs) {
  int64_t i;
  check_realloc_s(scales, sizeof(float), num_inputs);
  scales[num_inputs-1] = scale;
  check_realloc_s(dt, sizeof(float), num_inputs);
  float prev_scale = 0;
  if (num_inputs>1) prev_scale = scales[num_inputs-2];
  dt[num_inputs-1] = scale_to_years(scale) - scale_to_years(prev_scale);
  check_realloc_s(rem, sizeof(float), num_inputs);
  for (i=0; i<num_inputs; i++) {
    prev_scale = 0;
    if (i>0) prev_scale = scales[i-1];
    rem[i] = calc_sm_loss_int(prev_scale, scales[i], scale);
  }
}


void read_params(char *buffer, double *params, int max_n) {
  int num_entries = 0;
  char *cur_pos = buffer, *end_pos;
  float val = strtod(cur_pos, &end_pos);
  while (cur_pos != end_pos && num_entries < max_n) {
    params[num_entries] = val;
    num_entries++;
    cur_pos=end_pos;
    while (*cur_pos==' ' || *cur_pos=='\t' || *cur_pos=='\n') cur_pos++;
    val = strtod(cur_pos, &end_pos);
  }
}

int sort_by_mf(const void *a, const void *b) {
  const struct catalog_halo *c = a;
  const struct catalog_halo *d = b;
  if (c->mf < d->mf) return -1;
  if (d->mf < c->mf) return 1;
  return 0;
}

void load_halos(char *filename, int num_inputs) {
  int64_t i;
  struct packed_halo *p = NULL;
  if (old) free(old);
  if (sfh_old) free(sfh_old);
  if (icl_old) free(icl_old);
  old = new;
  sfh_old = sfh_new;
  icl_old = icl_new;
  num_old = num_new;
  int64_t new_length=0;
  p = check_mmap_file(filename, 'r', &new_length);
  assert(p && new_length > 0 && !(new_length % sizeof(struct packed_halo)));
  num_new = new_length / sizeof(struct packed_halo);
  new = NULL;
  sfh_new = NULL;
  icl_new = NULL;
  check_realloc_s(new, sizeof(struct catalog_halo), num_new);
  check_realloc_s(sfh_new, sizeof(float), num_new*num_inputs);
  check_realloc_s(icl_new, sizeof(float), num_new*num_inputs);
  memset(new, 0, sizeof(struct catalog_halo)*num_new);
  memset(sfh_new, 0, sizeof(float)*num_new*num_inputs);  
  memset(icl_new, 0, sizeof(float)*num_new*num_inputs);
  memset(counts, 0, sizeof(int64_t)*M_BINS);
  memset(bin_start, 0, sizeof(int64_t)*(M_BINS+1));
  for (i=0; i<num_new; i++) {
    new[i].id = p[i].id;
    new[i].descid = p[i].descid;
    memcpy(new[i].pos, p[i].pos, sizeof(float)*6);
    new[i].m = p[i].m;
    new[i].mp = p[i].mp;
    new[i].v = p[i].v;
    new[i].vp = p[i].vp;
    new[i].vmp = p[i].vmp;
    new[i].r = p[i].r;
    new[i].sat = p[i].sat;
    new[i].acc_q = p[i].acc_q;
    new[i].acc_qne = p[i].acc_qne;
    new[i].a_half = p[i].a_half;
    new[i].a_4p = p[i].a_4p;
    new[i].a_12 = p[i].a_12;
    new[i].c = p[i].c;
    new[i].pc = p[i].pc;
    new[i].acc_scale = p[i].acc_scale;
    new[i].peak_scale = p[i].peak_scale;
    new[i].spin_acc = p[i].spin_acc;
    new[i].va = p[i].va;
    new[i].a_v05 = p[i].a_v05;
    new[i].mf = (log10(new[i].mp/H0)-M_MIN)*BPDEX;
    assert(new[i].mf >= 0 && new[i].mf < M_BINS);
    counts[(int64_t)new[i].mf]++;
    if (!i || new[i].id > max_id) max_id = new[i].id;
    if (!i || new[i].id < min_id) min_id = new[i].id;
  }

  bin_start[0] = 0;
  for (i=1; i<M_BINS+1; i++) bin_start[i] = bin_start[i-1]+counts[i-1];
  memset(counts, 0, sizeof(int64_t)*M_BINS);
  struct catalog_halo *sorted_new = NULL;
  check_realloc_s(sorted_new, sizeof(struct catalog_halo), num_new);
  for (i=0; i<num_new; i++) {
    int64_t mb = new[i].mf;
    sorted_new[bin_start[mb]+counts[mb]] = new[i];
    counts[mb]++;
  }
  free(new);
  new = sorted_new;
  
  check_realloc_s(id_to_idx, sizeof(int64_t), (max_id+1-min_id));
  for (i=0; i<max_id+1-min_id; i++) id_to_idx[i] = -1;
  for (i=0; i<num_new; i++) id_to_idx[new[i].id-min_id] = i;
  munmap(p, new_length);

  for (i=0; i<num_old; i++) {
    if (old[i].descid < 0) continue;
    int64_t desc_idx = id_to_index(old[i].descid);
    if (desc_idx > -1 && old[i].mp > new[desc_idx].prog_mp) new[desc_idx].prog_mp = old[i].mp;
  }
}

int64_t id_to_index(int64_t id) {
  if (id < min_id || id > max_id) return -1;
  //assert(id >= min_id && id <= max_id);
  int64_t idx = id_to_idx[id-min_id];
  //assert(idx > -1);
  return idx;
}


void write_catalog(int64_t n) {
  char buffer[1024];
  snprintf(buffer, 1024, "sfr_catalog_%f.bin", scales[n]);
  FILE *output = check_fopen(buffer, "w");
  check_fwrite(new, sizeof(struct catalog_halo), num_new, output);
  fclose(output);
}
