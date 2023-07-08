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
#include "mah.h"
#include "mt_rand.h"

#define H0 0.7

double merger_rate_mpeak(double mass, double ratio, double scale);

struct packed_halo {
  int64_t id, descid;
  float pos[6], m, v, mp, vp, r, sat;
  float acc_q, acc_qne;
};

#define NUM_Z0_HALOS 10000
#define START_MASS 8
#define END_MASS 15

#define Gc (4.30117902e-9)  //Actually, Gc * (Msun / Mpc) in (km/s)^2
#define CRITICAL_DENSITY 2.77519737e11 // 3H^2/8piG in (Msun / h) / (Mpc / h)^3

double vir_density(double a) {
  double Om = 0.27;
  double Ol = 0.73;
  double x = (Om/pow(a,3))/(Om/pow(a,3) + Ol) - 1.0;
  return ((18*M_PI*M_PI + 82.0*x - 39*x*x)/(1.0+x))*Om*CRITICAL_DENSITY;
}

float Dyn_time(float a) {
  double vir = vir_density(a)*pow(a, -3);
  double t = 9.77813106e11/(sqrt((4.0/3.0)*M_PI*Gc*vir)*H0);
  return t;
}

float quenching_timescale(float hm, float a) {
  //  return 8.71e9*exp(-9.68674e-6*pow(10, 0.46794*sm))*Dyn_time(a)/Dyn_time(1);
  return 1.06e10*exp(-5.57145e-5*pow(10, 0.347749*hm))*Dyn_time(a)/Dyn_time(1);
}


int main(int argc, char **argv) {
  int64_t i, n;
  FILE *out, *in;
  char buffer[1024];
  //  struct smf_fit the_smf;

  if (argc < 2) {
    fprintf(stderr, "Usage: %s scales.txt\n", argv[0]);
    exit(1);
  }

  r250_init(87L);
  init_time_table(0.27, H0);

  struct packed_halo *h = NULL;
  check_realloc_s(h, sizeof(struct packed_halo), NUM_Z0_HALOS);
  memset(h, 0, sizeof(struct packed_halo)*NUM_Z0_HALOS);

  for (i=0; i<NUM_Z0_HALOS; i++) {
    h[i].id = i;
    h[i].m = (double)END_MASS + ((double)(START_MASS-END_MASS))*((double)i/(double)NUM_Z0_HALOS);
    h[i].r = normal_random(0, 1);
  }
  int64_t num_halos = NUM_Z0_HALOS;
  int64_t new_halo_id = NUM_Z0_HALOS;

  float scales[1000];
  in = check_fopen(argv[1], "r");
  int64_t num_scales=0;
  while (fgets(buffer, 1024, in)) scales[num_scales++] = atof(buffer);
  fclose(in);

  for (n=num_scales-1; n>=0; n--) {
    for (i=0; i<num_halos; i++) h[i].descid = h[i].id;

    //merger_rate_mpeak(double mass, double ratio, double scale) {
    //Add mergers
    if (n<num_scales-1) {
      int64_t new_halos = 0;
      for (i=0; i<num_halos; i++) {
	float m, mi = log10(h[i].mp/0.7)-0.1;
	for (m=mi; m>8; m-=0.2) {
	  double mr = merger_rate_mpeak(mi, m-mi, scales[n])*0.2*(1/scales[n]-1/scales[n+1]);
	  if (dr250() < mr) {
	    if (!(new_halos%1000)) 
	      check_realloc_s(h, sizeof(struct packed_halo), (num_halos+new_halos+1000));
	    memset(h+num_halos+new_halos, 0, sizeof(struct packed_halo));
	    h[num_halos+new_halos].id = new_halo_id++;
	    h[num_halos+new_halos].descid = h[i].id;
	    float mz0 = m_now(scales[n], m);
	    //printf("Added new halo (%f -> %f), a=%f!!!\n", m, mi, scales[n]);
	    h[num_halos+new_halos].m = mz0;
	    h[num_halos+new_halos].r = 0;
	    new_halos++;
	  }
	}
      }
    }

    for (i=0; i<num_halos; i++) {
      float m_now = m_at_a(h[i].m, scales[n])+m_scatter_at_a(h[i].m, scales[n])*h[i].r;
      h[i].mp = pow(10, m_now)*0.7;
      float qt = quenching_timescale(m_now, scales[n]);
      float scale_at_qt = years_to_scale(scale_to_years(scales[n])-qt);
      float m_at_qt = 0;
      if (scale_at_qt > 0)
	m_at_qt = m_at_a(h[i].m, scale_at_qt)+m_scatter_at_a(h[i].m, scale_at_qt)*h[i].r;
      h[i].acc_q = (h[i].mp - m_at_qt)/qt;
    }

    for (i=0; i<num_halos; i++) {
      if (h[i].mp < 1e8) {
	num_halos--;
	h[i] = h[num_halos];
	i--;
      }
    }

    snprintf(buffer, 1024, "acc_list_%f.bin", scales[n]);
    out = check_fopen(buffer, "w");
    fwrite(h, sizeof(struct packed_halo), num_halos, out);
    fclose(out);
  }
  return 0;
}

