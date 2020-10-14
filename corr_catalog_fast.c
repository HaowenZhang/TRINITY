#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <stdint.h>
#include <sys/mman.h>
#include "check_syscalls.h"
#define MINIMAL_CORR_STRUCT
#include "corr.h"
#include "stringparse.h"
#include "make_sf_catalog.h"

float BOX_SIZE = 0;
float HALF_BOX_SIZE = 0;
float RMAX = 40;
float MAX_D2 = 0;
//int64_t PERIODIC = 1;

#define FAST3TREE_TYPE struct point
#define FAST3TREE_DIM 3
#include "fast3tree.c"


inline int dist_range(struct tree3_node *n1, struct tree3_node *n2, float *min_dist, float *max_dist) {
  float d[3], e[3];
  int64_t i;
  for (i=0; i<3; i++) {
    float d1 = fabs(fabs(n2->min[i] - n1->max[i])-HALF_BOX_SIZE);
    float d2 = fabs(fabs(n1->min[i] - n2->max[i])-HALF_BOX_SIZE);
    if (d1 < d2) {
      float c = d2;
      d2 = d1;
      d1 = c;
    }
    d[i] = HALF_BOX_SIZE - d1;
    e[i] = HALF_BOX_SIZE - d2;
  }
  if (d[2] > RMAX) return 0;
  *min_dist = d[0]*d[0] + d[1]*d[1];
  *max_dist = e[0]*e[0] + e[1]*e[1];
  if (e[2] > RMAX) return 1;
  return 2;
}

inline int dist_bin(float d) {
  return (ilogb(d*1e6));
}

int calc_div(struct tree3_node *n) {
  int div = 0, divs = 0;
  while (n->parent != n) {
    n = n->parent;
    if (n->div_dim > 1) continue;
    div = div*2 + n->div_dim;
    divs++;
  }
  if (divs != F_LOG2_DIVS) return 0;
  return div+1;
}

void _calc_xcorr(struct tree3_node *n1, int64_t div1, struct tree3_node *n2, int64_t div2, int64_t dd[F_NUM_BINS+1][F_TOTAL_DIV+1], int64_t counts[F_TOTAL_DIV+1]) {
  float min_dist, max_dist;
  int res = dist_range(n1, n2, &min_dist, &max_dist);
  int64_t i;
  if (!res || min_dist > MAX_D2) return;

  if (!div1) {
    div1 = calc_div(n1);
    if (div1) counts[div1] += n1->num_points;
  }
  if (!div2) {
    div2 = calc_div(n2);
    if (div2) counts[div2] += n2->num_points;
  }

  int64_t tp = n2->num_points*n1->num_points;
  //Check if all points in n1 and n2 are within same bin
  if (res==2) {
    int64_t mb = dist_bin(min_dist);
    if (mb==dist_bin(max_dist)) {
      if (mb < 0) return;
      dd[mb][div1] += tp;
      dd[mb][div2] += tp;
      return;
    }
    if (tp < F_NUM_N2) {
      float offsets[2] = {0};
      for (i=0; i<2; i++) {
	if (n1->min[i] - n2->min[i] > HALF_BOX_SIZE) offsets[i] = -BOX_SIZE;
	if (n1->min[i] - n2->min[i] < -HALF_BOX_SIZE) offsets[i] = BOX_SIZE;
      }
      for (i=0; i<n1->num_points; i++) {
	int64_t j = (n1==n2) ? i+1 : 0;
	for (; j<n2->num_points; j++) {
	  float d0 = n1->points[i].pos[0] - n2->points[j].pos[0] + offsets[0];
	  float d1 = n1->points[i].pos[1] - n2->points[j].pos[1] + offsets[1];
	  int64_t mb = dist_bin(d0*d0+d1*d1);
	  if (mb < 0) continue;
	  dd[mb][div1]++;
	  dd[mb][div2]++;
	}
      }
      return;
    }
  }

  //Do n^2 counting
  if (tp < F_NUM_N2 || (n2->div_dim < 0 && n1->div_dim < 0)) {
    float offsets[3] = {0};
    for (i=0; i<3; i++) {
      if (n1->min[i] - n2->min[i] > HALF_BOX_SIZE) offsets[i] = -BOX_SIZE;
      if (n1->min[i] - n2->min[i] < -HALF_BOX_SIZE) offsets[i] = BOX_SIZE;
    }
    for (i=0; i<n1->num_points; i++) {
      int64_t j = (n1==n2) ? i+1 : 0;
      for (; j<n2->num_points; j++) {
	float d0 = n1->points[i].pos[0] - n2->points[j].pos[0] + offsets[0];
	float d1 = n1->points[i].pos[1] - n2->points[j].pos[1] + offsets[1];
	float d2 = n1->points[i].pos[2] - n2->points[j].pos[2] + offsets[2];
	if (fabs(d2) > RMAX) continue;
	int64_t mb = dist_bin(d0*d0+d1*d1);
	if (mb < 0 || mb >= NUM_BINS) continue;
	dd[mb][div1]++;
	dd[mb][div2]++;
      }
    }
    return;
  }

  //If same node, special splitting is necessary
  if (n1 == n2) {
    _calc_xcorr(n1->left, div1, n2->left, div2, dd, counts);
    _calc_xcorr(n1->left, div1, n2->right, div2, dd, counts);
    _calc_xcorr(n1->right, div1, n2->right, div2, dd, counts);
    return;
  }

  //Otherwise, split on the largest node
  if ((((n1->max[0]-n1->min[0])*(n1->max[1]-n1->min[1])) <
       ((n2->max[0]-n2->min[0])*(n2->max[1]-n2->min[1]))) && (n2->div_dim > -1)) {
    struct tree3_node *n3 = n2;
    int64_t div3 = div2;
    n2 = n1;
    div2 = div1;
    n1 = n3;
    div1 = div3;
  }
  
  _calc_xcorr(n1->left, div1, n2, div2, dd, counts);
  _calc_xcorr(n1->right, div1, n2, div2, dd, counts);
}

void calc_corr(struct fast3tree *t, int64_t num_points, char *filename) {
  int64_t dd[F_NUM_BINS+1][F_TOTAL_DIV+1], corr[F_NUM_BINS+1][F_TOTAL_DIV+1];
  int64_t i, j, k, div, div_bin, counts[F_TOTAL_DIV+1];
  double dens, v1, v2, ravg, v;

  FILE *out = check_fopen(filename, "w");
  fprintf(out, "#R Wp (Err)\n");
  fprintf(out, "#Num objects: %"PRId64"\n", num_points);

  dens = num_points / (BOX_SIZE*BOX_SIZE*BOX_SIZE);
  for (j=0; j<F_NUM_BINS+1; j++) 
    for (div=0; div<F_TOTAL_DIV+1; div++) dd[j][div] = counts[div] = 0;
  HALF_BOX_SIZE = BOX_SIZE / 2.0;
  MAX_D2 = F_MAX_DIST*F_MAX_DIST;
  
  _calc_xcorr(t->root, 0, t->root, 0, dd, counts);

  for (j=0; j<F_NUM_BINS; j++)
    for (div=1; div<F_TOTAL_DIV; div++)
      dd[j][0] += dd[j][div];

  for (j=0; j<NUM_BINS-1; j++) {
    float r1 = F_MIN_DIST*pow(2.0,j/2.0);
    float r2 = F_MIN_DIST*pow(2.0,(j+1)/2.0);
    v1 = M_PI * r1*r1;
    v2 = M_PI * r2*r2;
    ravg = sqrt((v1+v2)/(2.0*M_PI));
    v = v2 - v1;
    corr[j][0] = 
      (dd[j][0])/(double)num_points / (dens*v) - 2.0*RMAX;
    /*    for (div=0; div<TOTAL_DIV+1; div++) {
      corr[j][div] = (counts[j][div] && counts[j+1][div]) ? 
	(((dd[j+1][div]/(double)counts[j+1][div]-dd[j][div]/(double)counts[j][div]) / (dens*v))-BOX_SIZE) : 0;
	}
    var = 0;
    for (div=0; div<TOTAL_DIV; div++) {
      dx = corr[j][div]-corr[j][TOTAL_DIV];
      var += dx*dx;
    }
    var *= (double)(TOTAL_DIV - 1)/(double)(TOTAL_DIV);
    if (corr[j][TOTAL_DIV] > -1)
    */
    //fprintf(out, "%f %f %f\n", ravg/0.7, corr[j][TOTAL_DIV]/0.7, sqrt(var)/0.7);
    fprintf(out, "%f %f\n", ravg/0.7, corr[j][0]/0.7);
  }
  fclose(out);
}

int main(int argc, char **argv)
{
  struct point *points = NULL;
  int64_t num_points = 0;
  int64_t i,j,q;
  char buffer[1024];
  struct point p;
  struct fast3tree *tree = NULL;

  double lgm[3] = {9.8, 10.2, 10.6};
  
  tree = fast3tree_init(0,NULL);
  if (argc < 3) {
    fprintf(stderr, "Usage: %s sfr_catalog box_size rmax\n", argv[0]);
    exit(1);
  }

  BOX_SIZE = atof(argv[2]);
  if (argc > 3) RMAX = atof(argv[3]);

  struct catalog_halo *ch = NULL;
  int64_t length = 0, num_halos;
  ch = check_mmap_file(argv[1], 'r', &length);
  assert(!(length % sizeof(struct catalog_halo)));
  num_halos = length / sizeof(struct catalog_halo);
  mlock(ch, length);
  struct reduced_halo *rh = NULL, trh;
  int64_t num_rh=0;
  double mthresh = pow(10, lgm[0]);
  for (i=0; i<num_halos; i++) {
    if (ch[i].sm < mthresh) continue;
    check_realloc_every(rh, sizeof(struct reduced_halo), num_rh, 1000);
    memcpy(trh.pos, ch[i].pos, sizeof(float)*3);
    trh.sm = ch[i].sm;
    trh.nq = (ch[i].obs_sfr < 1e-11 * ch[i].sm) ? 0 : 1;
    rh[num_rh] = trh;
    num_rh++;    
  }
  munlock(ch, length);
  munmap(ch, length);

  char *names[3] = {"sf", "q", "all"};
  for (j=0; j<3; j++) {
    double mthresh = pow(10, lgm[j]);
    for (q=0; q<3; q++) {
      snprintf(buffer, 1024, "wp_sm%g_%s.dat", lgm[j], names[q]);
      num_points = 0;
      for (i=0; i<num_rh; i++) {
	if (rh[i].sm < mthresh) continue;
	if (q==rh[i].nq) continue;
	memcpy(p.pos, rh[i].pos, sizeof(float)*3);
	check_realloc_every(points, sizeof(struct point), num_points, 1000);
	points[num_points] = p;
	num_points++;
      }
      fast3tree_rebuild(tree, num_points, points);
      _fast3tree_set_minmax(tree, 0, BOX_SIZE);
      calc_corr(tree, num_points, buffer);
    }
  }
  return 0;
}

