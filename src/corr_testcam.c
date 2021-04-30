#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <stdint.h>
#include <sys/mman.h>
#include "mt_rand.h"
#include "check_syscalls.h"
#include "corr.h"
#include "stringparse.h"
#include "make_sf_catalog.h"

float BOX_SIZE = 0;
int64_t PERIODIC = 1;
#define NUM_LABELS 72
#define NUM_PROPS (NUM_LABELS+1)
#define LVMP (NUM_PROPS)
#define SAT 6
#define SM (NUM_PROPS+1)
#define SM_PSC (NUM_PROPS+4)
#define VBIN (NUM_PROPS+2)
#define Q (NUM_PROPS+3)
#define VX 20
#define JX 23
#define AX 45
#define AX500 50
#define BTOA 43
#define BTOA500c 48
char *labels[] = {"scale", "id", "desc_scale", "desc_id", "num_prog", "pid", "upid", "desc_pid", "phantom", "sam_mvir", "mvir", "rvir", "rs", "vrms", "mmp", "scale_of_last_MM", "vmax", "x", "y", "z", "-vabs", "vy", "vz", "Jabs", "Jy", "Jz", "Spin", "Breadth_first_ID", "Depth_first_ID", "Tree_root_ID", "Orig_halo_ID", "Snap_num", "Next_coprogenitor_depthfirst_ID", "Last_progenitor_depthfirst_ID", "Rs_Klypin", "Mvir_all", "M200b", "M200c", "M500c", "M2500c", "Xoff", "Voff", "Spin_Bullock", "b_to_a", "c_to_a", "Aabs", "Ay", "Az", "b_to_a500c", "c_to_a500c", "A500c_abs", "Ay500c", "Az500c", "TU", "M_pe_Behroozi", "M_pe_Diemer", "Macc", "Mpeak", "Vacc", "Vpeak", "Halfmass_Scale", "Acc_Rate_Inst", "Acc_Rate_100Myr", "Acc_Rate_1Tdyn", "Acc_Rate_2Tdyn", "Acc_Rate_Mpeak", "Mpeak_Scale", "Acc_Scale", "First_Acc_Scale", "First_Acc_Mvir", "First_Acc_Vmax", "VmaxMpeak", "VpeakVmax"};
int important[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1};


#define VMP 71
#define POS 17
#define VPVM 72
#define VM 16
#define VP 59

int64_t cam_param = 0;

struct halo {
  double properties[NUM_PROPS+5];
};

#define FAST3TREE_TYPE struct point
#define FAST3TREE_DIM 2
#include "fast3tree.c"


#define NORMAL_BINS (1024*8+1)
double p_to_normal[NORMAL_BINS];

struct sms_vmp {
  double sm, vmp;
};

struct sms_vmp *sms_vmps = NULL;
int64_t num_sms_vmp = 0;

double perc_to_normal(double p) {
  double x1 = -14;
  double x2 = 14;
  while (x2-x1 > 1e-7) {
    double half = 0.5*(x1+x2);
    double perc = 0.5*(1.0+erf(half));
    if (perc > p) { x2 = half; }
    else { x1 = half; }    
  }
  return ((x1+x2)*M_SQRT1_2);
}

void populate_pcache(void) {
  int64_t i=0;
  for (i=0; i<NORMAL_BINS; i++) {
    double p = (double)i / ((double)NORMAL_BINS-1);
    p_to_normal[i] = perc_to_normal(p);
  }
}

double pcache(double p) {
  double f = p*(NORMAL_BINS-1.0);
  if (f<0) return p_to_normal[0];
  if (f>=(NORMAL_BINS-1)) return p_to_normal[NORMAL_BINS-1];
  int64_t b = f;
  f -= b;
  return (p_to_normal[b] + f*(p_to_normal[b+1]-p_to_normal[b]));
}

int sort_by_vmp(const void *a, const void *b) {
  const struct sms_vmp *c = a;
  const struct sms_vmp *d = b;
  if (c->vmp < d->vmp) return -1;
  if (c->vmp > d->vmp) return 1;
  return 0;
}

void load_sm_vmp(char *file) {
  char buffer[1024];
  FILE *input;
  struct sms_vmp s;
  int64_t n;
  input = check_fopen(file, "r");
  while (fgets(buffer, 1024, input)) {
    n = sscanf(buffer, "%lf %lf", &s.sm, &s.vmp);
    if (n<2) continue;
    check_realloc_every(sms_vmps, sizeof(struct sms_vmp), num_sms_vmp, 1000);
    sms_vmps[num_sms_vmp] = s;
    num_sms_vmp++;
  }
  fclose(input);
  qsort(sms_vmps, num_sms_vmp, sizeof(struct sms_vmp), sort_by_vmp);
}

double assign_sm(double vmp) {
  double lvmp = log10(vmp);
  int64_t i;
  for (i=0; i<num_sms_vmp; i++)
    if (sms_vmps[i].vmp > lvmp) break;
  if (i == num_sms_vmp) return sms_vmps[num_sms_vmp-1].sm;
  if (i == 0) return sms_vmps[0].sm;
  i--;
  double f = (lvmp - sms_vmps[i].vmp)/(sms_vmps[i+1].vmp - sms_vmps[i].vmp);
  return (sms_vmps[i].sm + f*(sms_vmps[i+1].sm - sms_vmps[i].sm));
}

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

double fq_vmpeak(double lvmp) {
  return (1.0 - 1.0/(1.0+exp(5.1102*(lvmp-2.28082))));
}

//Sorts larger-to-smaller on vmp.
int sort_by_vmax_mpeak(const void *a, const void *b) {
  const struct halo *c = a;
  const struct halo *d = b;
  double vpa = c->properties[VMP];
  double vpb = d->properties[VMP];
  if (vpa > vpb) return -1;
  if (vpa < vpb) return 1;
  return 0;
}

double vec_abs(double *v) {
  int64_t i;
  double s = 0;
  for (i=0; i<3; i++) s+=v[i]*v[i];
  return sqrt(s);
}

int sort_cam(const void *a, const void *b) {
  const struct halo *c = a;
  const struct halo *d = b;
  double vpa = c->properties[cam_param];
  double vpb = d->properties[cam_param];
  if (vpa < vpb) return -1;
  if (vpa > vpb) return 1;
  return 0;
}

int main(int argc, char **argv)
{
  struct point *points = NULL;
  int64_t num_points = 0, i;
  char buffer[1024];
  struct point p;
  struct fast3tree *tree = NULL;
  double scatter = 0;
  double covariance = 0;
  char *out_dir;

  double lgm[3] = {9.8, 10.2, 10.6};
  //double counts[3] = {54969, 82654, 70314};
  populate_pcache();

  tree = fast3tree_init(0,NULL);
  if (argc < 7) {
    fprintf(stderr, "Usage: %s hlist box_size sm_vmp.dat scatter covariance out_dir\n", argv[0]);
    exit(1);
  }

  out_dir = argv[6];
  BOX_SIZE = atof(argv[2]);
  PERIODIC = 1;
  struct halo d;
  enum parsetype types[NUM_LABELS];
  void *data[NUM_LABELS];
  SHORT_PARSETYPE;
  for (i=0; i<NUM_LABELS; i++) {
    types[i] = F64;
    data[i] = d.properties+i;
  }

  load_sm_vmp(argv[3]);
  scatter = atof(argv[4]);
  covariance = atof(argv[5]);
  r250_init(87L);

  struct halo *halos = NULL;
  int64_t num_halos = 0;
  FILE *in = check_fopen(argv[1], "r");
  while (fgets(buffer, 1024, in)) {
    if (buffer[0] == '#') continue;
    int64_t n = stringparse(buffer, data, types, NUM_LABELS);
    if (n<NUM_LABELS) continue;
    d.properties[VX] = -1.0*vec_abs(d.properties+VX);
    d.properties[JX] = vec_abs(d.properties+JX);
    d.properties[AX] = vec_abs(d.properties+AX);
    d.properties[AX500] = vec_abs(d.properties+AX500);
    d.properties[SM] = assign_sm(d.properties[VMP]);
    d.properties[BTOA] = -d.properties[BTOA]; 
    d.properties[BTOA+1] = -d.properties[BTOA+1]; 
    d.properties[BTOA500c] = -d.properties[BTOA500c]; 
    d.properties[BTOA500c+1] = -d.properties[BTOA500c+1]; 
    if (d.properties[SM] < 9.8 - scatter*4.0) continue;
    d.properties[VPVM] = (d.properties[VM]-d.properties[VP])/d.properties[VP];
    check_realloc_every(halos, sizeof(struct halo), num_halos, 1000);
    halos[num_halos] = d;
    num_halos++;
  }
  fclose(in);

  qsort(halos, num_halos, sizeof(struct halo), sort_by_vmax_mpeak);
  double vmin = 1e4;
  double vmax = -1e4;
  for (i=0; i<num_halos; i++) {
    double lgvmp = log10(halos[i].properties[VMP]);
    halos[i].properties[LVMP] = lgvmp;
    if (lgvmp>vmax) vmax = lgvmp;
    if (lgvmp<vmin) vmin = lgvmp;
  }

  int64_t total_bins = ((vmax-vmin) / 0.02)+1;
  double bpdex = 1.0/0.02;
  int64_t *bcounts = NULL, *offsets = NULL;
  double *av_lgvmp = NULL, *qf = NULL;
  check_realloc_s(bcounts, sizeof(int64_t), total_bins);
  check_realloc_s(offsets, sizeof(int64_t), total_bins);
  check_realloc_s(av_lgvmp, sizeof(double), total_bins);
  check_realloc_s(qf, sizeof(double), total_bins);
  memset(bcounts, 0, sizeof(int64_t)*total_bins);
  memset(av_lgvmp, 0, sizeof(double)*total_bins);
  memset(qf, 0, sizeof(double)*total_bins);
  for (i=0; i<total_bins; i++) offsets[i] = -1;

  for (i=0; i<num_halos; i++) {
    int64_t bin = (halos[i].properties[LVMP]-vmin)*bpdex;
    if (offsets[bin] < 0) offsets[bin] = i;
    halos[i].properties[VBIN] = bin;
    bcounts[bin]++;
    av_lgvmp[bin] += halos[i].properties[LVMP];
  }

  for (i=0; i<total_bins; i++) {
    if (!bcounts[i]) continue;
    av_lgvmp[i] /= (double)bcounts[i];
    qf[i] = fq_vmpeak(av_lgvmp[i]);
  }

  int64_t j,k,q;

  double corr_scatter = covariance*scatter;
  double random_scatter = sqrt(1.0-covariance*covariance)*scatter;

  snprintf(buffer, 1024, "%s/sat_stats.dat", out_dir);
  FILE *sat_stats = check_fopen(buffer, "w");
  fprintf(sat_stats, "#Label Count %%Sat %%Sat_Q %%Cen_Q (all for 10 < log M* < 10.5)\n");
  for (k=0; k<NUM_PROPS; k++) {
    if (!important[k]) continue;
    //Decide which halos are q/sf
    cam_param = k;
    for (i=0; i<total_bins; i++) {
      int64_t q = 0;
      if (!bcounts[i]) continue;
      qsort(halos+offsets[i], bcounts[i], sizeof(struct halo), sort_cam);
      double sep = halos[(int64_t)(offsets[i]+qf[i]*bcounts[i])].properties[cam_param];
      for (j=offsets[i]; j<offsets[i]+bcounts[i]; j++) {
	halos[j].properties[Q] = 0;
	if (halos[j].properties[cam_param]<sep) { halos[j].properties[Q] = 1; q++; }
	halos[j].properties[SM_PSC] = halos[j].properties[SM];
	if (corr_scatter) halos[j].properties[SM_PSC] += pcache(((double)(j-offsets[i]))/(double)bcounts[i])*corr_scatter;
	if (random_scatter) halos[j].properties[SM_PSC] += normal_random(0, random_scatter);
      }
      //printf("%f %f %"PRId64" %"PRId64"\n", av_lgvmp[i], qf[i], q, bcounts[i]);
    }

    //Correlate!
    char *names[3] = {"sf", "q", "all"};
    int64_t count_c=0, qc=0, count_s=0, qs=0;
    for (j=0; j<3; j++) { //SM Thresh
      for (q=0; q<3; q++) {
	FILE *catalog = NULL;
	if (q==2 && j==0) {
	  snprintf(buffer, 1024, "%s/catalog_%s.dat", out_dir, labels[k]);
	  catalog = check_fopen(buffer, "w");
	  fprintf(catalog, "#ID Vpeak@Mpeak %s Q? Sat? SM X Y Z VX VY VZ\n", labels[k]);
	}
	snprintf(buffer, 1024, "%s/wp_sm_%s_%g_%s.dat", out_dir, labels[k], lgm[j], names[q]);
	num_points = 0;
	for (i=0; i<num_halos; i++) {
	  if (catalog) {
	    fprintf(catalog, "%g %g %g %g %d %g %g %g %g %g %g %g\n", halos[i].properties[1], halos[i].properties[VMP], halos[i].properties[cam_param], halos[i].properties[Q], (halos[i].properties[SAT] > -1) ? 1 : 0, halos[i].properties[SM_PSC], halos[i].properties[POS], halos[i].properties[POS+1], halos[i].properties[POS+2], halos[i].properties[POS+3], halos[i].properties[POS+4], halos[i].properties[POS+5]);
	  }
	  if (j==0 && q==2 && (halos[i].properties[SM_PSC] > 10.0) && (halos[i].properties[SM_PSC] < 10.5)) {
	    if (halos[i].properties[SAT] < 0) { count_c++; if (halos[i].properties[Q]) qc++; }
	    else { count_s++; if (halos[i].properties[Q]) qs++; }
	  }
	  if (halos[i].properties[SM_PSC]<lgm[j]) continue;
	  if (q<2 && q!=halos[i].properties[Q]) continue;
	  check_realloc_every(points, sizeof(struct point), num_points, 1000);
	  p.pos[0] = halos[i].properties[POS];
	  p.pos[1] = halos[i].properties[POS+1];
	  p.div = DIVISIONS*((int64_t)(DIVISIONS*p.pos[0] / BOX_SIZE))
	    + ((int64_t)(DIVISIONS*p.pos[1] / BOX_SIZE));
	  points[num_points] = p;
	  num_points++;
	}
	if (catalog) fclose(catalog);
	fast3tree_rebuild(tree, num_points, points);
	calc_corr(tree, num_points, buffer);
      }
    }
    fprintf(sat_stats, "%s %"PRId64" %f %f %f\n", labels[k], count_c+count_s, (double)count_s / ((double)(count_c+count_s)), (double)qs / (double)count_s, (double)qc / (double)count_c);
    fflush(sat_stats);
  }
  fclose(sat_stats);
  return 0;
}

