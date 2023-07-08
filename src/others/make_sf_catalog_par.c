#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <inttypes.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <assert.h>
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
#include <omp.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "make_sf_catalog.h"

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
int64_t counts[M_BINS];
int64_t bin_start[M_BINS+1];
double q_frac[M_BINS];
gsl_rng *rng[MAX_THREADS];

void calc_ages(float scale) {
  int64_t i;
#pragma omp for private(i) 
  for (i=0; i<num_new; i++) {
    new[i].age = scale*new[i].a_4p + (1.0-scale)*new[i].a_half;
    if (new[i].a_12 < new[i].age) new[i].age = new[i].a_12;
  }
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

  init_rngs();
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

  in = check_fopen(argv[3], "r");
  i = 0;
  while (fgets(buffer, 1024, in)) {
    char filename[1024];
    double scale = atof(buffer);
    snprintf(filename, 1024, "acc_list_%f.bin", scale);
    load_halos(filename, i+1);
#pragma omp parallel
    {
    calc_ages(scale);
#pragma omp single
    {
    calc_losses(atof(buffer), i+1);
    }
    calc_sm_histories(i);
#pragma omp single
    {
    calc_quenched_fraction(&the_smf, scale);
    quench_galaxies(i);
    }
    calc_sf_from_accretion_rates(i);
    }
    write_catalog(i);
    i++;
  }
  fclose(in);
  return 0;
}

void init_rngs(void) {
  int64_t i;
  for (i=0; i<MAX_THREADS; i++) {
    rng[i] = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng[i], 3141592653L);
  }
}

double normal_random_gsl(double mean, double stddev) {
  int64_t t = omp_get_thread_num();
  return (mean + gsl_ran_gaussian(rng[t], stddev));
}

void calc_sf_from_accretion_rates(int64_t n) {
  int64_t i,j;
  double avg_sfr[M_BINS] = {0};
  double icl_frac[M_BINS] = {0};
  double frac_icl_to_transfer[M_BINS] = {0};

#pragma omp single
{
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
}

#define SET(x,y) float x = y[mb]; if (counts[mb+1] && (avg_sfr[mb+1]>0)) x+=f*(y[mb+1]-y[mb]);
  double scatter_corr = exp(pow(INTR_SCATTER*log(10), 2)/2.0);
#pragma omp for private(j)
  for (j=0; j<num_new; j++) {
    float f = new[j].mf - 0.5;
    int64_t mb = f;
    f -= mb;
    if (!counts[mb]) { mb = new[j].mf; f = 0; }
    SET(asfr, avg_sfr);
    float obs_scatter = (new[j].q) ? OBS_SIGMA_Q : OBS_SIGMA_SF;

    if (new[j].q)
      new[j].sfr = pow(10, normal_random_gsl(Q_SSFR, INTR_SCATTER))*new[j].sm;
    else
      new[j].sfr = asfr*pow(10, normal_random_gsl(0, INTR_SCATTER))/scatter_corr;

    new[j].obs_sfr = new[j].sfr*pow(10, normal_random_gsl(0, obs_scatter));
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
  float tmp_sfh[MAX_OUTPUTS] = {0};
  float tmp_icl[MAX_OUTPUTS] = {0};
  assert(MAX_OUTPUTS > num_outputs);
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
    sfh[j] = icl[j] = 0;
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
  sfh[j]=icl[j]=0;
}

void calc_sm_histories(int num_old_inputs) {
  int64_t i, j;
#pragma omp for private(i,j)
  for (i=0; i<num_old; i++) {
    if (old[i].descid < 0) continue;
    int64_t desc = id_to_index(old[i].descid);
    float *desc_sfh = sfh_new+(desc*(num_old_inputs+1));
    float *prog_sfh = sfh_old+(i*num_old_inputs);
    float *desc_icl = icl_new+(desc*(num_old_inputs+1));
    float *prog_icl = icl_old+(i*num_old_inputs);
    if (old[i].mp < new[desc].prog_mp) {//To ICL 
      if (new[desc].icl < 0) {
	for (j=0; j<num_old_inputs; j++) desc_icl[j]=prog_sfh[j];
	new[desc].icl = 0;
      } else {
	for (j=0; j<num_old_inputs; j++) desc_icl[j]+=prog_sfh[j];
      }
    }
    else {
      if (new[desc].icl < 0) {
	for (j=0; j<num_old_inputs; j++) desc_icl[j]=prog_icl[j];
	new[desc].icl = 0;
      } else {
	for (j=0; j<num_old_inputs; j++) desc_icl[j]+=prog_icl[j];
      }
      for (j=0; j<num_old_inputs; j++) desc_sfh[j]=prog_sfh[j];
      new[desc].q = old[i].q;
      desc_sfh[j] = desc_icl[j] = 0;
    }
  }

#pragma omp for private(i,j)
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
  memset(sfh_new, 0, sizeof(float)*num_new*num_inputs);
  memset(icl_new, 0, sizeof(float)*num_new*num_inputs);
  memset(counts, 0, sizeof(int64_t)*M_BINS);
  memset(bin_start, 0, sizeof(int64_t)*(M_BINS+1));

  struct catalog_halo n = {0};
  n.icl = -1;
  for (i=0; i<num_new; i++) {
    n.id = p[i].id;
    n.descid = p[i].descid;
    memcpy(n.pos, p[i].pos, sizeof(float)*6);
    n.m = p[i].m;
    n.mp = p[i].mp;
    n.v = p[i].v;
    n.vp = p[i].vp;
    n.r = p[i].r;
    n.sat = p[i].sat;
    n.acc_q = p[i].acc_q;
    n.acc_qne = p[i].acc_qne;
    n.a_half = p[i].a_half;
    n.a_4p = p[i].a_4p;
    n.a_12 = p[i].a_12;
    n.mf = (log10(n.mp/H0)-M_MIN)*BPDEX;
    assert(n.mf >= 0 && n.mf < M_BINS);
    counts[(int64_t)n.mf]++;
    if (!i || n.id > max_id) max_id = n.id;
    if (!i || n.id < min_id) min_id = n.id;
    new[i] = n;
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
    if (old[i].mp > new[desc_idx].prog_mp) new[desc_idx].prog_mp = old[i].mp;
  }
}

int64_t id_to_index(int64_t id) {
  assert(id >= min_id && id <= max_id);
  int64_t idx = id_to_idx[id-min_id];
  assert(idx > -1);
  return idx;
}


void write_catalog(int64_t n) {
  char buffer[1024];
  snprintf(buffer, 1024, "sfr_catalog_%f.bin", scales[n]);
  FILE *output = check_fopen(buffer, "w");
  check_fwrite(new, sizeof(struct catalog_halo), num_new, output);
  fclose(output);
}

float quenched_fraction(float sm, float sm0) {
  if (!(sm>0)) return 0;
  return 1.0/(pow(sm/sm0, -1.3)+1.0);
}

void calc_quenched_fraction(struct smf_fit *t, float scale) {
  int64_t i;
  float sm0 = pow(10, 10.2+0.5*(1.0/scale-1.0));
  struct smf c = smhm_at_z(1.0/scale - 1.0, *t);
  float obs_scatter = sqrt(c.scatter*c.scatter + c.obs_scatter*c.obs_scatter);
  for (i=0; i<M_BINS; i++) {
    float m = M_MIN + ((i+0.5)*INV_BPDEX);
    float sm = calc_sm_at_m(m, c);
    if (m>12)
      sm += obs_scatter*(m-12);
    q_frac[i] = quenched_fraction(pow(10, sm), sm0);
  }
}


int sort_by_age(const void *a, const void *b) {
  const float *c = a;
  const float *d = b;
  if (*c < *d) return -1;
  if (*c > *d) return 1;
  return 0;
}			      
