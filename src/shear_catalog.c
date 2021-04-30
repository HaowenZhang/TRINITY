#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <sys/mman.h>
#include <gsl/gsl_spline.h>
#include "make_sf_catalog.h"
#include "check_syscalls.h"
#include "stringparse.h"

#define BPDEX 8
#define MIN_R -2
#define NUM_BINS 20

struct shear {
  int64_t id;
  int32_t bins[NUM_BINS];
};
int64_t *idx = NULL;
int64_t min_id, max_id;
struct shear *sh = NULL;
int64_t num_sh = 0;

int main(int argc, char **argv) {
  if (argc < 3) {
    fprintf(stderr, "Usage: %s sfr_catalog shear_catalog pmass\n", argv[0]);
    exit(1);
  }

  float pmass = 1.36297e8/0.7*pow(2048, 3) / 100e6;
  if (argc > 3) pmass = atof(argv[3]);

  char buffer[1024];
  FILE *in = check_fopen(argv[2], "r");
  struct shear s;
  int64_t i, j;
  void *data[NUM_BINS+4];
  enum parsetype pt[NUM_BINS+4];
  data[0] = &s.id;
  pt[0] = PARSE_INT64;
  data[1] = data[2] = data[3] = NULL;
  pt[1] = pt[2] = pt[3] = PARSE_SKIP;
  for (i=0; i<NUM_BINS; i++) {
    data[i+4] = s.bins+i;
    pt[i+4] = PARSE_INT32;
  }

  while (fgets(buffer, 1024, in)) {
    if (buffer[0] == '#') continue;
    if (stringparse(buffer, data, pt, NUM_BINS+4)!=(NUM_BINS+4)) continue;
    check_realloc_every(sh, sizeof(struct shear), num_sh, 1000);
    sh[num_sh] = s;
    if (!num_sh) {
      min_id = max_id = s.id;
    } else {
      if (min_id > s.id) min_id = s.id;
      if (max_id < s.id) max_id = s.id;
    }
    num_sh++;
  }
  fclose(in);

  idx = check_realloc(idx, (max_id-min_id+1)*sizeof(int64_t), "ID-index");
  for (i=0; i<max_id-min_id+1; i++) idx[i] = -1;
  for (i=0; i<num_sh; i++) idx[sh[i].id-min_id] = i;


  int64_t length = 0;
  struct catalog_halo *ch = check_mmap_file(argv[1], 'r', &length);
  mlock(ch, length);
  assert(!(length % sizeof(struct catalog_halo)));
  int64_t num_halos = length / sizeof(struct catalog_halo);
  int64_t total_shear[3][2][NUM_BINS];
  int64_t counts[3][2];
  float lgm[3] = {9.8, 10.2, 10.6};
  float masses[3];

  int64_t k, q;
  for (i=0; i<3; i++)
    for (j=0; j<2; j++) {
      for (k=0; k<NUM_BINS; k++)
	total_shear[i][j][k] = 0;
      counts[i][j] = 0;
    }
  for (i=0; i<3; i++) masses[i] = pow(10, lgm[i]);

  for (i=0; i<num_halos; i++) {
    if (ch[i].sm < masses[0]) continue;
    for (j=0; j<3; j++) {
      q = (ch[i].sm*1e-11 > ch[i].obs_sfr) ? 1 : 0;
      int64_t offset = ch[i].id - min_id;
      assert(ch[i].id <= max_id && ch[i].id >= min_id);
      struct shear *s = sh+idx[offset];
      if (ch[i].sm > masses[j]) {
	counts[j][q]++;
	for (k=0; k<NUM_BINS; k++) {
	  total_shear[j][q][k] += s->bins[k];
	  //if (!k && s->bins[k]>1000) {
	    //	    printf("#WTF?? ID: %"PRId64"\n", ch[i].id);
	  //}
	}
      }
    }
  }
  munmap(ch, length);

  double r[NUM_BINS];
  double t_r[NUM_BINS];
  double dens[3][2][NUM_BINS];
  double total_dens[3][2][NUM_BINS];
  gsl_spline *g_dens[3][2];
  gsl_spline *g_tdens[3][2];

  for (j=0; j<3; j++) {
    for (q=0; q<2; q++) {
      printf("#%f %"PRId64" %"PRId64, lgm[j], q, counts[j][q]);
      int64_t total = 0;
      for (k=0; k<NUM_BINS; k++) {
	double r2 = pow(10, MIN_R+k/(double)BPDEX)/0.7;
	double r1 = pow(10, MIN_R+(k-1)/(double)BPDEX)/0.7;
	if (k==0) r1 = 0;
	double a = M_PI*(r2*r2-r1*r1);
	double av_r = sqrt(r1*r1+ 0.5*(r2*r2-r1*r1));
	double total_a = M_PI*r2*r2;
	r[k] = av_r;
	t_r[k] = r2;
	dens[j][q][k] = total_shear[j][q][k]/(double)(counts[j][q]*a);
	total += total_shear[j][q][k];
	total_dens[j][q][k] = total/(double)(counts[j][q]*total_a);
	//printf(" %f(%f) %f(%f)", dens[j][q][k], total_shear[j][q][k]/(double)(counts[j][q]), total_dens[j][q][k], total/(double)(counts[j][q]));
      }
      printf("\n");
      g_dens[j][q] = gsl_spline_alloc(gsl_interp_akima, NUM_BINS);
      g_tdens[j][q] = gsl_spline_alloc(gsl_interp_akima, NUM_BINS);
      gsl_spline_init(g_dens[j][q], r, dens[j][q], NUM_BINS);
      gsl_spline_init(g_tdens[j][q], t_r, total_dens[j][q], NUM_BINS);
    }
  }

  gsl_interp_accel ga = {0};
  printf("#R");
  for (j=0; j<3; j++)
    for (q=0; q<2; q++)
      printf(" SM(%f,%"PRId64")", lgm[j], q);
  printf("\n");

  double norm = 1e-12;
  for (k=1; k<NUM_BINS; k++) {
    double r2 = pow(10, MIN_R+k/(double)BPDEX)/0.7;
    double r1 = pow(10, MIN_R+(k-1)/(double)BPDEX)/0.7;
    if (k==0) r1 = 0;
    double av_r = sqrt(r1*r1+ 0.5*(r2*r2-r1*r1));
    printf("%f", av_r);
    for (j=0; j<3; j++) {
      for (q=0; q<2; q++) {
	double r_dens = gsl_spline_eval(g_dens[j][q], av_r, &ga);
	double enc_dens = gsl_spline_eval(g_tdens[j][q], av_r, &ga);
	printf(" %e", norm*pmass*(enc_dens-r_dens));
      }
    }
    printf("\n");
  }
  return 0;
}
