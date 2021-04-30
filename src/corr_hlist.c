#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <stdint.h>
#include <sys/mman.h>
#include "check_syscalls.h"
#include "corr.h"
#include "stringparse.h"
//#include "make_sf_catalog.h"

float BOX_SIZE = 0;
int64_t PERIODIC = 1;

#define FAST3TREE_TYPE struct point
#define FAST3TREE_DIM 2
#include "fast3tree.c"

void calc_corr(struct fast3tree *t, int64_t num_points, char *filename) {
  double dd[NUM_BINS][TOTAL_DIV+1], corr[NUM_BINS][TOTAL_DIV+1];
  int64_t i, j, k, div, div_bin, counts[NUM_BINS][TOTAL_DIV+1];
  double dens, dist, v1, v2, ravg, v, var=0, dx;
  struct fast3tree_results *res = NULL;

  FILE *out = check_fopen(filename, "w");
  fprintf(out, "#R Wp (Err)\n");
  fprintf(out, "#Num objects: %"PRId64"\n", num_points);

  dens = num_points / (BOX_SIZE*BOX_SIZE*BOX_SIZE);
  for (j=0; j<NUM_BINS; j++) 
    for (div=0; div<TOTAL_DIV+1; div++) dd[j][div] = counts[j][div] = 0;

  res = fast3tree_results_init();
  for (i=0; i<num_points; i++) {
    for (j=0; j<NUM_BINS; j++) {
      dist = pow(10, MIN_DIST + j*INV_DIST_BPDEX);
      int64_t np = 0;
      if (PERIODIC) {
	fast3tree_find_sphere_periodic(t, res, t->root->points[i].pos, dist);
	np = res->num_points;
      }
      else {
	for (k=0; k<2; k++)
	  if (t->root->points[i].pos[k] < dist || t->root->points[i].pos[k] > BOX_SIZE-dist) break;
	if (k!=2) continue;
	fast3tree_find_sphere(t, res, t->root->points[i].pos, dist);
	double d = t->root->points[i].pos[2];
	double d2 = BOX_SIZE - t->root->points[i].pos[2];
	if (d > d2) d = d2;
	for (k=0; k<res->num_points; k++) {
	  double dx = fabs(res->points[k]->pos[2] - t->root->points[i].pos[2]);
	  if (dx*2.0 > BOX_SIZE) continue;
	  if (dx > d) np+=2;
	  else np++;
	}
      }
      for (div = 0; div < TOTAL_DIV; div++) {
	div_bin = (div != t->root->points[i].div) ? div : TOTAL_DIV;
	dd[j][div_bin]+=np;
	counts[j][div_bin]++;
      }
    }
  }

  for (j=0; j<NUM_BINS-1; j++) {
    v1 = M_PI * pow(10, 2*(MIN_DIST + j*INV_DIST_BPDEX));
    v2 = M_PI * pow(10, 2*(MIN_DIST + (j+1)*INV_DIST_BPDEX));
    ravg = sqrt((v1+v2)/(2.0*M_PI));
    v = v2 - v1;
    for (div=0; div<TOTAL_DIV+1; div++) {
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
      fprintf(out, "%f %f %f\n", ravg/1.0, corr[j][TOTAL_DIV]/1.0, sqrt(var)/1.0);
  }

  fprintf(out, "\n");
  fclose(out);
  fast3tree_results_free(res);
}

int sort_by_vmp(const void *a, const void *b) {
  const struct reduced_halo *c = a;
  const struct reduced_halo *d = b;
  if (c->sm > d-> sm) return -1;
  if (c->sm < d-> sm) return 1;
  return 0;
}


int main(int argc, char **argv)
{
  struct point *points = NULL;
  int64_t num_points = 0, i;
  char buffer[1024];
  struct point p;
  struct fast3tree *tree = NULL;

  int64_t ncounts[3] = {300000, 180000, 70000};
  double vmp_limits[3] = {0,0,0};
  
  tree = fast3tree_init(0,NULL);
  if (argc < 3) {
    fprintf(stderr, "Usage: %s hlist.list box_size\n", argv[0]);
    exit(1);
  }

  BOX_SIZE = atof(argv[2]);

  struct reduced_halo *ch = NULL;
  int64_t num_halos=0;
  struct reduced_halo h = {{0}};

  SHORT_PARSETYPE;
  #define NUM_INPUTS 4
  struct parse_format pf[NUM_INPUTS] =     
    {{17, F, &(h.pos[0])}, {18, F, &(h.pos[1])}, {19, F, &(h.pos[2])},
     {71, F, &(h.sm)}};
  
  FILE *input = check_fopen(argv[1], "r");
  while (fgets(buffer, 1024, input)) {
    int64_t n = stringparse_format(buffer, pf, NUM_INPUTS);
    if (n!=NUM_INPUTS) continue;
    if (h.sm < 100) continue;
    check_realloc_every(ch, sizeof(struct reduced_halo), num_halos, 1000);
    ch[num_halos] = h;
    num_halos++;
  }
  fclose(input);
  printf("%"PRId64"\n", num_halos);

  qsort(ch, num_halos, sizeof(struct reduced_halo), sort_by_vmp);
  for (i=0; i<3; i++)
    vmp_limits[i] = ch[ncounts[i]].sm;

  double orig_box_size = BOX_SIZE;
  int64_t bin = 0, j, k;
  check_realloc_s(points, sizeof(struct point), ncounts[0]);
  for (j=0; j<3; j++) {
    PERIODIC = 0;
    BOX_SIZE = orig_box_size / 5.0;
    double vthresh = vmp_limits[j];
    for (bin=0; bin<=125; bin++) {
      int64_t b = bin, bi[3];
      double min[3], max[3];
      for (i=0; i<3; i++) {
	bi[i] = b % 5;
	b /= 5;
	min[i] = bi[i]*BOX_SIZE;
	max[i] = min[i]+BOX_SIZE;
      }
      sprintf(buffer, "corr/corr_vt%.1f_%"PRId64"_%"PRId64"_%"PRId64".dat", vthresh, bi[0], bi[1], bi[2]);

      if (bin == 125) {
	BOX_SIZE = orig_box_size;
	PERIODIC = 1;
	for (i=0; i<3; i++) { min[i] = 0; max[i] = BOX_SIZE; }
	sprintf(buffer, "corr/corr_vt%.1f_all.dat", vthresh);
      }

      num_points = 0;
      for (i=0; i<ncounts[j]; i++) {
	for (k=0; k<3; k++)
	  if ((ch[i].pos[k] < min[k]) || (ch[i].pos[k] >= max[k])) break;
	if (k<3) continue;
	p.pos[0] = ch[i].pos[0] - min[0];
	p.pos[1] = ch[i].pos[1] - min[1];
	p.pos[2] = ch[i].pos[2] - min[2];
	p.div = DIVISIONS*((int64_t)(DIVISIONS*p.pos[0] / BOX_SIZE))
	  + ((int64_t)(DIVISIONS*p.pos[1] / BOX_SIZE));
	points[num_points] = p;
	num_points++;
      }
      //printf("b: %"PRId64" (%.0f, %.0f, %.0f) - (%.0f, %.0f, %.0f) %"PRId64"\n", bin, min[0], min[1], min[2], max[0], max[1], max[2], num_points);
      fast3tree_rebuild(tree, num_points, points);
      calc_corr(tree, num_points, buffer);
    }
  }
  return 0;
}

