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
#include "make_sf_catalog.h"

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
      if (PERIODIC)
	fast3tree_find_sphere_periodic(t, res, t->root->points[i].pos, dist);
      else {
	for (k=0; k<2; k++)
	  if (t->root->points[i].pos[k] < dist || t->root->points[i].pos[k] > BOX_SIZE-dist) break;
	if (k!=2) continue;
	fast3tree_find_sphere(t, res, t->root->points[i].pos, dist);
      }
      for (div = 0; div < TOTAL_DIV; div++) {
	div_bin = (div != t->root->points[i].div) ? div : TOTAL_DIV;
	dd[j][div_bin]+=res->num_points;
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
      fprintf(out, "%f %f %f\n", ravg/0.7, corr[j][TOTAL_DIV]/0.7, sqrt(var)/0.7);
  }

  fprintf(out, "\n");
  fclose(out);
  fast3tree_results_free(res);
}

int main(int argc, char **argv)
{
  struct point *points = NULL;
  int64_t num_points = 0;
  char buffer[1024];
  struct point p;
  struct fast3tree *tree = NULL;

  double lgm[3] = {9.8, 10.2, 10.6};
  
  tree = fast3tree_init(0,NULL);
  if (argc < 3) {
    fprintf(stderr, "Usage: %s sfr_catalog box_size [x,y translation if not periodic]\n", argv[0]);
    exit(1);
  }

  BOX_SIZE = atof(argv[2]);
  double ds[2]={0};
  if (argc>4) {
    PERIODIC = 0;
    ds[0] = atof(argv[3]);
    ds[1] = atof(argv[4]);
  }

  struct catalog_halo *ch = NULL;
  int64_t length = 0, num_halos;
  ch = check_mmap_file(argv[1], 'r', &length);
  assert(!(length % sizeof(struct catalog_halo)));
  num_halos = length / sizeof(struct catalog_halo);
  mlock(ch, length);

  int64_t i,j,q;
  char *names[3] = {"sf", "q", "all"};
  for (j=0; j<3; j++) {
    double mthresh = pow(10, lgm[j]);
    for (q=0; q<3; q++) {
      snprintf(buffer, 1024, "wp_sm%g_%s.dat", lgm[j], names[q]);
      num_points = 0;
      for (i=0; i<num_halos; i++) {
	if (ch[i].sm < mthresh) continue;
	float ssfr = ch[i].obs_sfr / ch[i].sm;
	if (q==0 && ssfr < 1e-11) continue;
	if (q==1 && ssfr > 1e-11) continue;
	check_realloc_every(points, sizeof(struct point), num_points, 1000);
	p.pos[0] = ch[i].pos[0] - ds[0];
	p.pos[1] = ch[i].pos[1] - ds[1];
	p.div = DIVISIONS*((int64_t)(DIVISIONS*p.pos[0] / BOX_SIZE))
	  + ((int64_t)(DIVISIONS*p.pos[1] / BOX_SIZE));
	points[num_points] = p;
	num_points++;
      }
      fast3tree_rebuild(tree, num_points, points);
      calc_corr(tree, num_points, buffer);
    }
  }
  munmap(ch, length);
  return 0;
}

